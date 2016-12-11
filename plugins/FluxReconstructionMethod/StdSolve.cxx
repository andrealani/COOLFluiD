// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/MethodCommandProvider.hh"

#include "Framework/CFSide.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"

#include "MathTools/MathFunctions.hh"

#include "FluxReconstructionMethod/StdSolve.hh"
#include "FluxReconstructionMethod/FluxReconstruction.hh"
#include "FluxReconstructionMethod/FluxReconstructionElementData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider< StdSolve,FluxReconstructionSolverData,FluxReconstructionModule >
  stdSolveProvider("StdSolve");
  
//////////////////////////////////////////////////////////////////////////////
  
StdSolve::StdSolve(const std::string& name) :
  FluxReconstructionSolverCom(name),
  socket_gradients("gradients"),
  socket_normals("normals"),
  m_cellBuilder(CFNULL),
  m_faceBuilder(CFNULL),
  m_solPntsLocalCoords(),
  m_iElemType(),
  m_cell(),
  m_cellStates(),
  m_faceNormals(),
  m_faceFlxPntConn(),
  m_nbrEqs(),
  m_dim(),
  m_orient(),
  m_face(),
  m_cells(),
  m_riemannFluxComputer(),
  m_states(),
  m_faceConnPerOrient(),
  m_flxPntRiemannFlux(),
  m_corrFlxFactor(),
  m_contFlx(),
  m_divContFlx()
  {
  }
  
//////////////////////////////////////////////////////////////////////////////

std::vector< Common::SafePtr< BaseDataSocketSink > >
StdSolve::needsSockets()
{
  std::vector< Common::SafePtr< BaseDataSocketSink > > result;
  result.push_back(&socket_gradients);
  result.push_back(&socket_normals);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

void StdSolve::execute()
{
  CFAUTOTRACE;
  
  CFLog(INFO, "StdSolve::execute()\n");
  
  // get the elementTypeData
  SafePtr< vector<ElementTypeData> > elemType = MeshDataStack::getActive()->getElementTypeData();

  // get InnerCells TopologicalRegionSet
  SafePtr<TopologicalRegionSet> cells = MeshDataStack::getActive()->getTrs("InnerCells");

  // get the geodata of the geometric entity builder and set the TRS
  StdTrsGeoBuilder::GeoData& geoData = m_cellBuilder->getDataGE();
  geoData.trs = cells;
  
  // get the local FR data
  vector< FluxReconstructionElementData* >& frLocalData = getMethodData().getFRLocalData();
  
  // get InnerFaces TopologicalRegionSet
  SafePtr<TopologicalRegionSet> faces = MeshDataStack::getActive()->getTrs("InnerFaces");

  // get the face start indexes
  vector< CFuint >& innerFacesStartIdxs = getMethodData().getInnerFacesStartIdxs();

  // get number of face orientations
  const CFuint nbrFaceOrients = innerFacesStartIdxs.size()-1;

  // get the geodata of the face builder and set the TRSs
  FaceToCellGEBuilder::GeoData& geoData2 = m_faceBuilder->getDataGE();
  geoData2.cellsTRS = cells;
  geoData2.facesTRS = faces;
  geoData2.isBoundary = false;
  
  // get number of solution points
  CFuint nbrSolPnt = frLocalData[0]->getNbrOfSolPnts();
    
  // get number of flux points
  CFuint nbrFlxPnt1D = (frLocalData[0]->getFlxPntsLocalCoord1D())->size();
  
  // get solution point local coordinates
  m_solPntsLocalCoords = frLocalData[0]->getSolPntsLocalCoords();
    
  // get flux point local coordinates
  m_flxPntsLocalCoords = frLocalData[0]->getFlxPntsLocalCoords();
    
  // get the face normals
  m_faceNormals = frLocalData[0]->getFaceNormals();
   
  // get the face - flx pnt connectivity per orient
  m_faceFlxPntConn = frLocalData[0]->getFaceFlxPntConnPerOrient();
    
  // get the face connectivity per orientation
  m_faceConnPerOrient = frLocalData[0]->getFaceConnPerOrient();
    
  SafePtr< ConvectiveVarSet > updateVarSet = getMethodData().getUpdateVar();

  // loop over different orientations
  for (m_orient = 0; m_orient < nbrFaceOrients; ++m_orient)
  {
    // start and stop index of the faces with this orientation
    const CFuint faceStartIdx = innerFacesStartIdxs[m_orient  ];
    const CFuint faceStopIdx  = innerFacesStartIdxs[m_orient+1];

    // loop over faces with this orientation
    for (CFuint faceID = faceStartIdx; faceID < faceStopIdx; ++faceID)
    {
      // build the face GeometricEntity
      geoData2.idx = faceID;
      m_face = m_faceBuilder->buildGE();

      // get the neighbouring cells
      m_cells[LEFT ] = m_face->getNeighborGeo(LEFT );
      m_cells[RIGHT] = m_face->getNeighborGeo(RIGHT);

      // get the states in the neighbouring cells
      m_states[LEFT ] = m_cells[LEFT ]->getStates();
      m_states[RIGHT] = m_cells[RIGHT]->getStates();
      
      (m_cellStatesFlxPnt[LEFT])->resize(nbrFlxPnt1D);
      (m_cellStatesFlxPnt[RIGHT])->resize(nbrFlxPnt1D);
      
      // get solution polynomial values at nodes
      vector< vector< CFreal > > solPolyValsAtNodes
	  = frLocalData[0]->getSolPolyValsAtNode(*m_flxPntsLocalCoords);
	  
      RealVector normalL = (*m_faceNormals)[(*m_faceConnPerOrient)[m_orient][LEFT]];
      RealVector normalR = (*m_faceNormals)[(*m_faceConnPerOrient)[m_orient][RIGHT]];
      vector< RealVector > normalListL;
      normalListL.resize(nbrFlxPnt1D);
      for (CFuint i = 0; i < nbrFlxPnt1D; ++i)
      {
	normalListL[i] = normalL;
      }
	  
      // Loop over flux points to extrapolate the states to the flux points and get the 
      // discontinuous normal flux in the flux points and the Riemann flux
      for (CFuint iFlxPnt = 0; iFlxPnt < nbrFlxPnt1D; ++iFlxPnt)
      {
	CFuint flxPntIdxL = (*m_faceFlxPntConn)[m_orient][LEFT][iFlxPnt];
	CFuint flxPntIdxR = (*m_faceFlxPntConn)[m_orient][RIGHT][iFlxPnt];
	
        for (CFuint iSol = 0; iSol < nbrSolPnt; ++iSol)
        {
          *(*(m_cellStatesFlxPnt[LEFT]))[iFlxPnt] += solPolyValsAtNodes[flxPntIdxL][iSol]*(*(*(m_states[LEFT]))[iSol]);
	  *(*(m_cellStatesFlxPnt[RIGHT]))[iFlxPnt] += solPolyValsAtNodes[flxPntIdxR][iSol]*(*(*(m_states[RIGHT]))[iSol]);
        }
        
        RealVector pData;
	// compute the normal flux at the current flux point
	updateVarSet->computePhysicalData(*(*(m_cellStatesFlxPnt[LEFT]))[iFlxPnt], pData);
	(*(m_cellFlx[LEFT]))[iFlxPnt] = updateVarSet->getFlux()(pData,normalL);
	updateVarSet->computePhysicalData(*(*(m_cellStatesFlxPnt[RIGHT]))[iFlxPnt], pData);
	(*(m_cellFlx[RIGHT]))[iFlxPnt] = updateVarSet->getFlux()(pData,normalR);
      }
      
      m_flxPntRiemannFlux = m_riemannFluxComputer->computeFlux(*m_cellStatesFlxPnt[LEFT],*m_cellStatesFlxPnt[RIGHT],
							       normalListL,nbrFlxPnt1D);

      for (CFuint jFlxPnt = 0; jFlxPnt < nbrFlxPnt1D; ++jFlxPnt)
      {
        (*(m_cellFlx[LEFT]))[jFlxPnt] = m_flxPntRiemannFlux[jFlxPnt] - (*(m_cellFlx[LEFT]))[jFlxPnt];
	(*(m_cellFlx[RIGHT]))[jFlxPnt] = -m_flxPntRiemannFlux[jFlxPnt] + (*(m_cellFlx[RIGHT]))[jFlxPnt];
      }
      
      CFuint leftID = (m_cells[LEFT])->getID();
      CFuint rightID = (m_cells[RIGHT])->getID();
      CFuint leftLocalFace = (*m_faceConnPerOrient)[m_orient][LEFT];
      CFuint rightLocalFace = (*m_faceConnPerOrient)[m_orient][RIGHT];
      
      (*(m_corrFlxFactor[leftID]))[leftLocalFace] = (*(m_cellFlx[LEFT]));
      (*(m_corrFlxFactor[rightID]))[rightLocalFace] = (*(m_cellFlx[RIGHT]));

      // release the GeometricEntity
      m_faceBuilder->releaseGE();
    }
  }
  
  // loop over element types
  const CFuint nbrElemTypes = elemType->size();
  for (m_iElemType = 0; m_iElemType < nbrElemTypes; ++m_iElemType)
  {
    // get start and end indexes for this type of element
    const CFuint startIdx = (*elemType)[m_iElemType].getStartIdx();
    const CFuint endIdx   = (*elemType)[m_iElemType].getEndIdx();
    
    // get number of solution points
    CFuint nbrSolPnt = frLocalData[m_iElemType]->getNbrOfSolPnts();
    
    // get number of flux points
    CFuint nbrFlxPnt = frLocalData[m_iElemType]->getNbrOfFlxPnts();
    
    // get number of flux points
    CFuint nbrFlxPnt1D = (frLocalData[0]->getFlxPntsLocalCoord1D())->size();
    
    // get number of faces
    CFuint nbrFaces = frLocalData[m_iElemType]->getNbrCellFaces();

    // get solution point local coordinates
    m_solPntsLocalCoords = frLocalData[m_iElemType]->getSolPntsLocalCoords();
    
    // get flux point local coordinates
    m_flxPntsLocalCoords = frLocalData[m_iElemType]->getFlxPntsLocalCoords();
    
    // get the face normals
    m_faceNormals = frLocalData[m_iElemType]->getFaceNormals();
    
    // get solution derivative polynomial values at solution points
    vector< vector< vector< CFreal > > > solDerivPolyValsAtNodes
	  = frLocalData[m_iElemType]-> getSolPolyDerivsAtNode(*m_solPntsLocalCoords);


    // loop over cells
    for (CFuint elemIdx = startIdx; elemIdx < endIdx; ++elemIdx)
    {
      // build the GeometricEntity
      geoData.idx = elemIdx;
      m_cell = m_cellBuilder->buildGE();

      // get the states in this cell
      m_cellStates = m_cell->getStates();
      
      m_contFlx.resize(nbrSolPnt);
      m_divContFlx.resize(nbrSolPnt);
        
      // Loop over solution points to calculate the discontinuous flux.
      for (CFuint iSolPnt = 0; iSolPnt < nbrSolPnt; ++iSolPnt)
      {
	
	RealVector pData;
	// compute the discontinuous flux at the current solution point
	updateVarSet->computePhysicalData(*(*m_cellStates)[iSolPnt], pData);
	m_contFlx[iSolPnt] = updateVarSet->getFlux()(pData);
      }
      
      // Loop over solution pnts to calculate the divergence of the discontinuous flux
      for (CFuint iSolPnt = 0; iSolPnt < nbrSolPnt; ++iSolPnt)
      {
	// Loop over solution pnt to count factor of all sol pnt polys
	for (CFuint jSolPnt = 0; jSolPnt < nbrSolPnt; ++jSolPnt)
	{
	  // Loop over deriv directions and sum them to compute divergence
	  for (CFuint iDir = 0; iDir < m_dim; ++iDir)
	  {
	    // Loop over conservative fluxes 
	    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
	    {
	      m_divContFlx[iSolPnt][iEq] -= solDerivPolyValsAtNodes[iSolPnt][iDir][jSolPnt]*(m_contFlx[iSolPnt])(iEq,iDir);
	    }
	  }
	}
      }
      for (CFuint iSolPnt = 0; iSolPnt < nbrSolPnt; ++iSolPnt)
      {
	
	RealVector Temp;
	Temp[0] = 1;
	Temp[1] = 1;
	Temp[2] = 1;
	CFreal TempDiv = 0;
	
	CFuint iFlxPnt = 0;
	for (CFuint iFace = 0; iFace < nbrFaces; ++iFace)
	{
	  for (CFuint iFlxPnt1D = 0; iFlxPnt1D < nbrFlxPnt1D; ++iFlxPnt1D, ++iFlxPnt)
	  {
	    RealVector currentCorrFactor = (*(m_corrFlxFactor[elemIdx]))[iFace][iFlxPnt1D];
	    cf_assert(currentCorrFactor.size() == m_nbrEqs);
	    
	    // Fill in the matrix
	    for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
	    {
	      for (CFuint iDim = 0; iDim < m_dim; ++iDim)
	      {
		(m_contFlx[iSolPnt])(iVar,iDim) += currentCorrFactor[iVar] * Temp[iDim];
	      }
	      m_divContFlx[iSolPnt][iVar] -= currentCorrFactor[iVar] * TempDiv; 
	    }
	  }
	}
      }


      //release the GeometricEntity
      m_cellBuilder->releaseGE();
    }
  }
  

  

  
  
}

//////////////////////////////////////////////////////////////////////////////

void StdSolve::setup()
{
  CFAUTOTRACE;
  FluxReconstructionSolverCom::setup();
  
  // get face builder
  m_faceBuilder = getMethodData().getFaceBuilder();
  
  // get cell builder
  m_cellBuilder = getMethodData().getStdTrsGeoBuilder();
  
  // get the Riemann flux
  m_riemannFluxComputer = getMethodData().getRiemannFlux();

//   // get and setup the face term computer
//   m_faceTermComputer = getMethodData().getFaceTermComputer();

  // dimensionality and number of equations
  m_dim   = PhysicalModelStack::getActive()->getDim();
  m_nbrEqs = PhysicalModelStack::getActive()->getNbEq();


  
}

//////////////////////////////////////////////////////////////////////////////

void StdSolve::unsetup()
{
  CFAUTOTRACE;
  FluxReconstructionSolverCom::unsetup();
}


//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD

