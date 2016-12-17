// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/MethodCommandProvider.hh"

#include "Framework/CFSide.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"
#include "Framework/BaseTerm.hh"

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
  socket_rhs("rhs"),
  socket_updateCoeff("updateCoeff"),
  socket_faceJacobVecSizeFaceFlxPnts("faceJacobVecSizeFaceFlxPnts"),
  m_updateVarSet(CFNULL),
  m_cellBuilder(CFNULL),
  m_faceBuilder(CFNULL),
  m_solPntsLocalCoords(CFNULL),
  m_faceIntegrationCoefs(CFNULL),
  m_faceMappedCoordDir(CFNULL),
  m_iElemType(),
  m_cell(),
  m_cellStates(),
  m_faceNormals(CFNULL),
  m_faceFlxPntConn(CFNULL),
  m_faceFlxPntConnPerOrient(CFNULL),
  m_nbrEqs(),
  m_dim(),
  m_orient(),
  m_face(),
  m_cells(),
  m_riemannFluxComputer(CFNULL),
  m_corrFctComputer(CFNULL),
  m_states(),
  m_faceConnPerOrient(CFNULL),
  m_flxPntRiemannFlux(CFNULL),
  m_corrFlxFactor(),
  m_contFlx(),
  m_cellFlx(),
  m_divContFlx(),
  m_corrFct(),
  m_corrFctDiv(),
  m_cellStatesFlxPnt(),
  m_faceJacobVecAbsSizeFlxPnts(),
  m_faceJacobVecSizeFlxPnts(),
  m_unitNormalFlxPnts(),
  m_cellFluxProjVects(),
  m_flxPntCoords(),
  m_waveSpeedUpd()
  {
    addConfigOptionsTo(this);
  }
  
  
//////////////////////////////////////////////////////////////////////////////

void StdSolve::defineConfigOptions(Config::OptionList& options)
{
}

//////////////////////////////////////////////////////////////////////////////

void StdSolve::configure ( Config::ConfigArgs& args )
{
  FluxReconstructionSolverCom::configure(args);
}  
  
//////////////////////////////////////////////////////////////////////////////

std::vector< Common::SafePtr< BaseDataSocketSink > >
StdSolve::needsSockets()
{
  std::vector< Common::SafePtr< BaseDataSocketSink > > result;
  result.push_back(&socket_gradients);
  result.push_back(&socket_normals);
  result.push_back(&socket_rhs);
  result.push_back(&socket_updateCoeff);
  result.push_back(&socket_faceJacobVecSizeFaceFlxPnts);
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
  const CFuint nbrSolPnt = frLocalData[0]->getNbrOfSolPnts();
    
  // get number of flux points
  const CFuint nbrFlxPnt1D = (frLocalData[0]->getFlxPntsLocalCoord1D())->size();
  
  // get solution point local coordinates
  m_solPntsLocalCoords = frLocalData[0]->getSolPntsLocalCoords();
    
  // get flux point local coordinates
  m_flxPntsLocalCoords = frLocalData[0]->getFlxPntsLocalCoords();
    
  // get the face normals
  m_faceNormals = frLocalData[0]->getFaceNormals();
   
  // get the face - flx pnt connectivity per orient
  m_faceFlxPntConnPerOrient = frLocalData[0]->getFaceFlxPntConnPerOrient();
    
  // get the face connectivity per orientation
  m_faceConnPerOrient = frLocalData[0]->getFaceConnPerOrient();
  
  // get the face integration coefficient
  m_faceIntegrationCoefs = frLocalData[0]->getFaceIntegrationCoefs();
  
  // get flux point mapped coordinate directions
  m_faceMappedCoordDir = frLocalData[0]->getFaceMappedCoordDirPerOrient();
  
  // get flux point local coordinates
  m_flxPntsLocalCoords1D = frLocalData[0]->getFlxPntsLocalCoord1D();
  
  // Resize vectors
  m_cells.resize(0);
  m_cells.resize(2);
  m_states.resize(0);
  m_states.resize(2);
  m_cellStatesFlxPnt.resize(0);
  m_cellStatesFlxPnt.resize(2);
  m_cellFlx.resize(0);
  m_cellFlx.resize(2);
  m_corrFlxFactor.resize(0);
  m_corrFlxFactor.resize(cells->getNbStatesInTrs()/nbrSolPnt);
  m_faceJacobVecAbsSizeFlxPnts.resize(nbrFlxPnt1D);
  
  for (CFuint iElem = 0; iElem < cells->getLocalNbGeoEnts(); ++iElem)
  {
    m_corrFlxFactor[iElem].resize(m_flxPntsLocalCoords->size()/nbrFlxPnt1D);
    for (CFuint iFace = 0; iFace < m_flxPntsLocalCoords->size()/nbrFlxPnt1D; ++iFace)
    {
      m_corrFlxFactor[iElem][iFace].resize(nbrFlxPnt1D);
      for (CFuint jFlxPnt = 0; jFlxPnt < nbrFlxPnt1D; ++jFlxPnt)
      {
	m_corrFlxFactor[iElem][iFace][jFlxPnt].resize(m_nbrEqs);
      }
    }
  }
  

  // loop over different orientations
  for (m_orient = 0; m_orient < nbrFaceOrients; ++m_orient)
  {
    //CFLog(VERBOSE, "Orient = " << m_orient << "\n");
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

      // Resize vectors
      m_cellStatesFlxPnt[LEFT].resize(0);
      m_cellStatesFlxPnt[RIGHT].resize(0);
      m_cellStatesFlxPnt[LEFT].resize(nbrFlxPnt1D);
      m_cellStatesFlxPnt[RIGHT].resize(nbrFlxPnt1D);
      m_cellFlx[LEFT].resize(0);
      m_cellFlx[RIGHT].resize(0);
      m_cellFlx[LEFT].resize(nbrFlxPnt1D);
      m_cellFlx[RIGHT].resize(nbrFlxPnt1D);
      
      // get solution polynomial values at nodes
      vector< vector< CFreal > > solPolyValsAtNodes
	  = frLocalData[0]->getSolPolyValsAtNode(*m_flxPntsLocalCoords);
	  
      RealVector normalL = (*m_faceNormals)[(*m_faceConnPerOrient)[m_orient][LEFT]];
      RealVector normalR = (*m_faceNormals)[(*m_faceConnPerOrient)[m_orient][RIGHT]];
      for (CFuint iDim = 0; iDim < m_dim; ++iDim)
      {
	cf_assert(normalL[iDim] == -normalR[iDim]);
      }
      
      // compute flux point coordinates
      //vector< RealVector > normalListL;
      //normalListL.resize(nbrFlxPnt1D);
      vector<RealVector> flxCoords1D;
      flxCoords1D.resize(nbrFlxPnt1D);
      for (CFuint iFlxPnt = 0; iFlxPnt < nbrFlxPnt1D; ++iFlxPnt)
      {
	flxCoords1D[iFlxPnt].resize(1);
	flxCoords1D[iFlxPnt][0] = (*(m_flxPntsLocalCoords1D))[iFlxPnt];
	m_flxPntCoords[iFlxPnt] = m_face->computeCoordFromMappedCoord(flxCoords1D[iFlxPnt]);
	//normalListL[iFlxPnt].resize(normalL.size());
	//normalListL[iFlxPnt] = normalL;
      }

      // compute face Jacobian vectors
      vector< RealVector > faceJacobVecs = m_face->computeFaceJacobDetVectorAtMappedCoords(flxCoords1D);
	  
      // Loop over flux points to extrapolate the states to the flux points and get the 
      // discontinuous normal flux in the flux points and the Riemann flux
      for (CFuint iFlxPnt = 0; iFlxPnt < nbrFlxPnt1D; ++iFlxPnt)
      {
	// get face Jacobian vector sizes in the flux points
	DataHandle< vector< CFreal > >
	  faceJacobVecSizeFaceFlxPnts = socket_faceJacobVecSizeFaceFlxPnts.getDataHandle();
        // get face Jacobian vector size
        m_faceJacobVecAbsSizeFlxPnts[iFlxPnt] = faceJacobVecSizeFaceFlxPnts[faceID][iFlxPnt];
	
	// set face Jacobian vector size with sign depending on mapped coordinate direction
	m_faceJacobVecSizeFlxPnts[iFlxPnt].resize(2);
        m_faceJacobVecSizeFlxPnts[iFlxPnt][LEFT] = m_faceJacobVecAbsSizeFlxPnts[iFlxPnt]*(*m_faceMappedCoordDir)[m_orient][LEFT];
	m_faceJacobVecSizeFlxPnts[iFlxPnt][RIGHT] = m_faceJacobVecAbsSizeFlxPnts[iFlxPnt]*(*m_faceMappedCoordDir)[m_orient][RIGHT];

	m_unitNormalFlxPnts[iFlxPnt].resize(m_dim);
	
	// set unit normal vector
        m_unitNormalFlxPnts[iFlxPnt] = faceJacobVecs[iFlxPnt]/m_faceJacobVecAbsSizeFlxPnts[iFlxPnt];
	
	CFuint flxPntIdxL = (*m_faceFlxPntConnPerOrient)[m_orient][LEFT][iFlxPnt];
	CFuint flxPntIdxR = (*m_faceFlxPntConnPerOrient)[m_orient][RIGHT][iFlxPnt];
	
	//m_cellStatesFlxPnt[LEFT][iFlxPnt].resize(m_nbrEqs);
	//m_cellStatesFlxPnt[RIGHT][iFlxPnt].resize(m_nbrEqs);
	m_cellFlx[LEFT][iFlxPnt].resize(0);
	m_cellFlx[RIGHT][iFlxPnt].resize(0);
	m_cellFlx[LEFT][iFlxPnt].resize(m_nbrEqs);
	m_cellFlx[RIGHT][iFlxPnt].resize(m_nbrEqs);
	
	State* tempVectorL = new State(solPolyValsAtNodes[flxPntIdxL][0]*(*((*(m_states[LEFT]))[0]->getData())));
	State* tempVectorR = new State(solPolyValsAtNodes[flxPntIdxR][0]*(*((*(m_states[RIGHT]))[0]->getData())));
	
	m_cellStatesFlxPnt[LEFT][iFlxPnt] = tempVectorL;
	m_cellStatesFlxPnt[RIGHT][iFlxPnt] = tempVectorR;

        for (CFuint iSol = 1; iSol < nbrSolPnt; ++iSol)
        {
          *(m_cellStatesFlxPnt[LEFT][iFlxPnt]) = (State) (*(m_cellStatesFlxPnt[LEFT][iFlxPnt]->getData()) + 
							  solPolyValsAtNodes[flxPntIdxL][iSol]*(*((*(m_states[LEFT]))[iSol]->getData())));
	  *(m_cellStatesFlxPnt[RIGHT][iFlxPnt]) = (State) (*(m_cellStatesFlxPnt[RIGHT][iFlxPnt]->getData()) +
							   solPolyValsAtNodes[flxPntIdxR][iSol]*(*((*(m_states[RIGHT]))[iSol]->getData())));
        }
        //CFLog(VERBOSE, "state in flx pnt = " << *(m_cellStatesFlxPnt[LEFT][iFlxPnt]->getData()) << "\n");
        
	// compute the normal flux at the current flux point
	m_updateVarSet->computePhysicalData(*(m_cellStatesFlxPnt[LEFT][iFlxPnt]), m_pData);
	m_cellFlx[LEFT][iFlxPnt] = m_updateVarSet->getFlux()(m_pData,m_unitNormalFlxPnts[iFlxPnt]);
	//CFLog(VERBOSE, "Flux Left = " << m_cellFlx[LEFT][iFlxPnt] << "\n");
	m_updateVarSet->computePhysicalData(*(m_cellStatesFlxPnt[RIGHT][iFlxPnt]), m_pData);
	m_cellFlx[RIGHT][iFlxPnt] = m_updateVarSet->getFlux()(m_pData,m_unitNormalFlxPnts[iFlxPnt]);
	//CFLog(VERBOSE, "Flux Right = " << m_cellFlx[RIGHT][iFlxPnt] << "\n");
	
      }
      m_flxPntRiemannFlux.resize(nbrFlxPnt1D);

      for (CFuint iFlxPnt = 0; iFlxPnt < nbrFlxPnt1D; ++iFlxPnt)
      {
	m_flxPntRiemannFlux[iFlxPnt].resize(m_nbrEqs);
	m_flxPntRiemannFlux[iFlxPnt] = m_riemannFluxComputer->computeFlux(*(m_cellStatesFlxPnt[LEFT][iFlxPnt]),
									  *(m_cellStatesFlxPnt[RIGHT][iFlxPnt]),
									  m_unitNormalFlxPnts[iFlxPnt]);
	//CFLog(VERBOSE, "RiemannFlux = " << m_flxPntRiemannFlux[iFlxPnt] << "\n");
      }
      
      for (CFuint jFlxPnt = 0; jFlxPnt < nbrFlxPnt1D; ++jFlxPnt)
      {
        m_cellFlx[LEFT][jFlxPnt] = (m_flxPntRiemannFlux[jFlxPnt] - m_cellFlx[LEFT][jFlxPnt])*m_faceJacobVecSizeFlxPnts[jFlxPnt][LEFT];
	m_cellFlx[RIGHT][jFlxPnt] = (m_flxPntRiemannFlux[jFlxPnt] - m_cellFlx[RIGHT][jFlxPnt])*m_faceJacobVecSizeFlxPnts[jFlxPnt][RIGHT];
	//CFLog(VERBOSE,"delta fluxL = " << m_cellFlx[LEFT][jFlxPnt] << "\n");
	//CFLog(VERBOSE,"delta fluxR = " << m_cellFlx[RIGHT][jFlxPnt] << "\n");
      }
      
      CFuint leftID = (m_cells[LEFT])->getID();
      CFuint rightID = (m_cells[RIGHT])->getID();
      CFuint leftLocalFace = (*m_faceConnPerOrient)[m_orient][LEFT];
      CFuint rightLocalFace = (*m_faceConnPerOrient)[m_orient][RIGHT];
      
      m_corrFlxFactor[leftID][leftLocalFace] = m_cellFlx[LEFT];
      m_corrFlxFactor[rightID][rightLocalFace] = m_cellFlx[RIGHT];

      computeWaveSpeedUpdates(m_waveSpeedUpd);

      // update the wave speeds
      updateWaveSpeed();

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
    
    // get the face - flx pnt connectivity per orient
    m_faceFlxPntConn = frLocalData[m_iElemType]->getFaceFlxPntConn();
    
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
      
      m_contFlx.resize(0);
      m_contFlx.resize(nbrSolPnt);
      m_divContFlx.resize(0);
      m_divContFlx.resize(nbrSolPnt);
      m_corrFctDiv.resize(0);
      m_corrFctDiv.resize(nbrSolPnt);
      m_cellFluxProjVects.resize(0);
      m_cellFluxProjVects.resize(m_dim);
      
      for (CFuint iDim = 0; iDim < m_dim; ++iDim)
      {
	m_cellFluxProjVects[iDim].resize(nbrSolPnt);
	vector<CFuint> dimList;
	dimList.resize(nbrSolPnt);
	for (CFuint iSolPnt = 0; iSolPnt < nbrSolPnt; ++iSolPnt)
	{
	  dimList[iSolPnt] = iDim;
	  m_cellFluxProjVects[iDim][iSolPnt].resize(m_dim);
	}
	
	m_cellFluxProjVects[iDim] = m_cell->computeMappedCoordPlaneNormalAtMappedCoords(dimList,
                                                                            *m_solPntsLocalCoords);
      }
      
      // Loop over solution points to calculate the discontinuous flux.
      for (CFuint iSolPnt = 0; iSolPnt < nbrSolPnt; ++iSolPnt)
      {
	m_contFlx[iSolPnt].resize(0);
        m_contFlx[iSolPnt].resize(m_dim);
	for (CFuint iDim = 0; iDim < m_dim; ++iDim)
	{
	  m_contFlx[iSolPnt][iDim].resize(m_nbrEqs);
	  // compute the discontinuous flux at the current solution point
	  m_updateVarSet->computePhysicalData(*(*m_cellStates)[iSolPnt], m_pData);
	  m_contFlx[iSolPnt][iDim] = m_updateVarSet->getFlux()(m_pData,m_cellFluxProjVects[iDim][iSolPnt]);
	  m_divContFlx[iSolPnt].resize(m_nbrEqs);
	  m_corrFctDiv[iSolPnt].resize(nbrFlxPnt);
	}
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
	      m_divContFlx[iSolPnt][iEq] -= solDerivPolyValsAtNodes[iSolPnt][iDir][jSolPnt]*(m_contFlx[iSolPnt][iDir][iEq]);
	      if (fabs(m_divContFlx[iSolPnt][iEq]) < MathTools::MathConsts::CFrealEps())
	      {
		m_divContFlx[iSolPnt][iEq] = 0;
	      }
	    }
	  }
	}
	//CFLog(VERBOSE, "Before: -div FD = " << m_divContFlx[iSolPnt] << "\n");
      }
      
      //m_corrFctComputer->computeCorrectionFunction(frLocalData[m_iElemType],m_corrFct);
      m_corrFctComputer->computeDivCorrectionFunction(frLocalData[m_iElemType],m_corrFctDiv);
      
      
      
      for (CFuint iSolPnt = 0; iSolPnt < nbrSolPnt; ++iSolPnt)
      {
	//CFLog(VERBOSE, "After: -div FD = " << m_divContFlx[iSolPnt] << "\n");
	for (CFuint iFace = 0; iFace < nbrFaces; ++iFace)
	{
	  for (CFuint iFlxPnt1D = 0; iFlxPnt1D < nbrFlxPnt1D; ++iFlxPnt1D)
	  {
	    RealVector currentCorrFactor = m_corrFlxFactor[elemIdx][iFace][iFlxPnt1D];
	    cf_assert(currentCorrFactor.size() == m_nbrEqs);
	    
	    // Fill in the matrix
	    for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
	    {
// 	      for (CFuint iDim = 0; iDim < m_dim; ++iDim)
// 	      {
// 		(m_contFlx[iSolPnt])(iVar,iDim) += currentCorrFactor[iVar] * m_corrFct[iSolPnt][iFlxPnt][iDim];
// 	      }
	      m_divContFlx[iSolPnt][iVar] -= currentCorrFactor[iVar] * m_corrFctDiv[iSolPnt][(*m_faceFlxPntConn)[iFace][iFlxPnt1D]]; 
	    }
	  }
	}
	if(m_cell->getID() == 323 || m_cell->getID() == 72)
	      {
		//CFLog(VERBOSE, "FI-FD = " << currentCorrFactor << "\n");
		//CFLog(VERBOSE, "div h = " << m_corrFctDiv[iSolPnt][iFlxPnt] << "\n");
		//CFLog(VERBOSE, "-div FC = " << m_divContFlx[iSolPnt] << "\n");
	      }
      }
      
      updateRHS();
      
      divideByJacobDet();
      
      
      if(m_cell->getID() == 323 || m_cell->getID() == 72)
	      {
		CFLog(VERBOSE, "ID  = " << (*m_cellStates)[0]->getLocalID() << "\n");
      CFLog(VERBOSE, "Update = \n");
      // get the datahandle of the rhs
      DataHandle< CFreal > rhs = socket_rhs.getDataHandle();
      for (CFuint iState = 0; iState < nbrSolPnt; ++iState)
      {
        CFuint resID = m_nbrEqs*( (*m_cellStates)[iState]->getLocalID() );
        for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
        {
          CFLog(VERBOSE, "" << rhs[resID+iVar] << " ");
    }
    CFLog(VERBOSE,"\n");

  }
	      }
	      
	      


      //release the GeometricEntity
      m_cellBuilder->releaseGE();
      

    }
  }
  

  

  
  
}


//////////////////////////////////////////////////////////////////////////////

void StdSolve::updateRHS()
{
  // get the datahandle of the rhs
  DataHandle< CFreal > rhs = socket_rhs.getDataHandle();

  // get residual factor
  const CFreal resFactor = getMethodData().getResFactor();
  
  const CFuint nbrStates = m_cellStates->size();

  // update rhs
  for (CFuint iState = 0; iState < nbrStates; ++iState)
  {
    CFuint resID = m_nbrEqs*( (*m_cellStates)[iState]->getLocalID() );
    for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
    {
      //CFLog(VERBOSE, "Update = " << m_divContFlx[iState][iVar] << "\n");
      rhs[resID+iVar] += resFactor*m_divContFlx[iState][iVar];
      //CFLog(VERBOSE,"-div F = " << m_divContFlx[iState][iVar] << "\n");
    }
//     if (std::abs(rhs[resID]) > 1e-6)
//     {
//       CF_DEBUG_OBJ(rhs[resID]);
//     }
  }
//   CF_DEBUG_POINT;
//   CF_DEBUG_EXIT;
}


//////////////////////////////////////////////////////////////////////////////

void StdSolve::updateWaveSpeed()
{
  // get the datahandle of the update coefficients
  DataHandle<CFreal> updateCoeff = socket_updateCoeff.getDataHandle();

  for (CFuint iSide = 0; iSide < 2; ++iSide)
  {
    const CFuint nbrSolPnts = m_states[iSide]->size();
    for (CFuint iSol = 0; iSol < nbrSolPnts; ++iSol)
    {
      const CFuint solID = (*m_states[iSide])[iSol]->getLocalID();
      updateCoeff[solID] += m_waveSpeedUpd[iSide];
      //CFLog(VERBOSE, "updateCoeff = " << updateCoeff[solID] << "\n");
    }
  }
  
}

//////////////////////////////////////////////////////////////////////////////

void StdSolve::computeWaveSpeedUpdates(vector< CFreal >& waveSpeedUpd)
{
  // compute the wave speed updates for the neighbouring cells
  cf_assert(waveSpeedUpd.size() == 2);
  for (CFuint iSide = 0; iSide < 2; ++iSide)
  {
    waveSpeedUpd[iSide] = 0.0;
    for (CFuint iFlx = 0; iFlx < m_cellFlx[iSide].size(); ++iFlx)
    {
      const CFreal jacobXIntCoef = m_faceJacobVecAbsSizeFlxPnts[iFlx]*
                                   (*m_faceIntegrationCoefs)[iFlx];
      // transform update states to physical data to calculate eigenvalues

// cout << m_flxPntSol[iSide][iFlx]->getLocalID() << " " << flush;

      m_updateVarSet->computePhysicalData(*(m_cellStatesFlxPnt[iSide][iFlx]), m_pData);
      waveSpeedUpd[iSide] += jacobXIntCoef*
          m_updateVarSet->getMaxAbsEigenValue(m_pData,m_unitNormalFlxPnts[iFlx]);
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void StdSolve::divideByJacobDet()
{
  // this is achieved by multiplying the update coefs with the Jacobian determinant
  // (and dividing by the cell volume)

  // get the updateCoeff
  DataHandle< CFreal > updateCoeff = socket_updateCoeff.getDataHandle();

//   // get the cell volume
//   const CFreal invCellVolume = 1.0/m_cell->computeVolume();

  // get jacobian determinants at solution points
  const std::valarray<CFreal> jacobDet =
      m_cell->computeGeometricShapeFunctionJacobianDeterminant(*m_solPntsLocalCoords);

  // get number of solution points
  const CFuint nbrSolPnts = m_cellStates->size();

  // loop over residuals
  for (CFuint iSol = 0; iSol < nbrSolPnts; ++iSol)
  {
    // get solution point ID
    const CFuint solID = (*m_cellStates)[iSol]->getLocalID();

    // divide update coeff by volume
    updateCoeff[solID] *= jacobDet[iSol];//*invCellVolume
  }
}


//////////////////////////////////////////////////////////////////////////////

void StdSolve::setup()
{
  CFAUTOTRACE;
  FluxReconstructionSolverCom::setup();
  
  // get the update varset
  m_updateVarSet = getMethodData().getUpdateVar();
  
  // get face builder
  m_faceBuilder = getMethodData().getFaceBuilder();
  
  // get cell builder
  m_cellBuilder = getMethodData().getStdTrsGeoBuilder();
  
  // get the Riemann flux
  m_riemannFluxComputer = getMethodData().getRiemannFlux();
  
  // get the correction function computer
  m_corrFctComputer = getMethodData().getCorrectionFunction();
  
  m_waveSpeedUpd.resize(2);
  
  // get the local spectral FD data
  vector< FluxReconstructionElementData* >& frLocalData = getMethodData().getFRLocalData();
  cf_assert(frLocalData.size() > 0);
  
  // number of flux points
  const CFuint nbrFlxPnts = frLocalData[0]->getFlxPntsLocalCoord1D()->size();
  
  // resize m_faceJacobVecSizeFlxPnts
  m_faceJacobVecSizeFlxPnts.resize(nbrFlxPnts);
  
  // resize m_unitNormalFlxPnts
  m_unitNormalFlxPnts.resize(nbrFlxPnts);
  

//   // get and setup the face term computer
//   m_faceTermComputer = getMethodData().getFaceTermComputer();

  // dimensionality and number of equations
  m_dim   = PhysicalModelStack::getActive()->getDim();
  m_nbrEqs = PhysicalModelStack::getActive()->getNbEq();
  
  // resize the physical data temporary vector
  SafePtr<BaseTerm> convTerm = PhysicalModelStack::getActive()->getImplementor()->getConvectiveTerm(); 
  convTerm->resizePhysicalData(m_pData);
  
  m_flxPntCoords.resize(nbrFlxPnts);
  for (CFuint iFlx = 0; iFlx < nbrFlxPnts; ++iFlx)
  {
    m_flxPntCoords[iFlx].resize(m_dim);
    
  }


  
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

