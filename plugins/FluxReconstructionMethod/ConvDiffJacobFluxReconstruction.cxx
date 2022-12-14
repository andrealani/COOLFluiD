// Copyright (C) 2016 KU Leuven, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Common/PE.hh"

#include "Framework/MethodCommandProvider.hh"
#include "Framework/BlockAccumulator.hh"
#include "Framework/LSSMatrix.hh"

#include "Framework/CFSide.hh"
#include "Framework/MeshData.hh"
#include "Framework/BaseTerm.hh"

#include "MathTools/MathFunctions.hh"

#include "FluxReconstructionMethod/ConvDiffJacobFluxReconstruction.hh"
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

MethodCommandProvider< ConvDiffJacobFluxReconstruction,
		       FluxReconstructionSolverData,
		       FluxReconstructionModule >
convDiffRHSJacobFluxReconstructionProvider("ConvDiffRHSJacob");

//////////////////////////////////////////////////////////////////////////////
  
ConvDiffJacobFluxReconstruction::ConvDiffJacobFluxReconstruction(const std::string& name) :
  DiffRHSJacobFluxReconstruction(name),
  m_updateVarSet(CFNULL),
  m_elemIdx(),
  m_facesCell(),
  m_jacob(),
  socket_solPntNormals("solPntNormals"),
  socket_flxPntNormals("flxPntNormals"),
  socket_cellVolumes("cellVolumes"),
  socket_volumes("volumes"),
  m_tempSolPntVec(),
  m_tempSolPntVec2(),
  m_fluxJacobian(),
  m_riemannFluxJacobian(),
  m_flxPntRiemannFluxDiff(),
  m_flxPntRiemannFluxPert(),
  m_tempFlux(),
  m_gradientFluxJacobian(),
  m_gradVarsToStateJacobian(),
  m_gradientStateJacobian(),
  m_riemannFluxGradJacobian(),
  m_contFlxBackupDiff()
  {
    addConfigOptionsTo(this);
  }

//////////////////////////////////////////////////////////////////////////////

void ConvDiffJacobFluxReconstruction::configure ( Config::ConfigArgs& args )
{
  DiffRHSJacobFluxReconstruction::configure(args);
}  

//////////////////////////////////////////////////////////////////////////////

void ConvDiffJacobFluxReconstruction::defineConfigOptions(Config::OptionList& options)
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector< Common::SafePtr< BaseDataSocketSource > >
  ConvDiffJacobFluxReconstruction::providesSockets()
{
  std::vector< Common::SafePtr< BaseDataSocketSource > > result = DiffRHSJacobFluxReconstruction::providesSockets();
  return result;
}

//////////////////////////////////////////////////////////////////////////////

std::vector< Common::SafePtr< BaseDataSocketSink > >
ConvDiffJacobFluxReconstruction::needsSockets()
{
  std::vector< Common::SafePtr< BaseDataSocketSink > > result = DiffRHSJacobFluxReconstruction::needsSockets();
  result.push_back(&socket_solPntNormals);
  result.push_back(&socket_flxPntNormals);
  result.push_back(&socket_cellVolumes);
  result.push_back(&socket_volumes);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffJacobFluxReconstruction::execute()
{
  CFAUTOTRACE;
  
  CFLog(VERBOSE, "ConvDiffJacobFluxReconstruction::execute()\n");
  
  ////////////////////INITIALIZATION/////////////////////////
  
  // get the elementTypeData
  SafePtr< vector<ElementTypeData> > elemType = MeshDataStack::getActive()->getElementTypeData();

  // get InnerCells TopologicalRegionSet
  SafePtr<TopologicalRegionSet> cells = MeshDataStack::getActive()->getTrs("InnerCells");

  // get the geodata of the geometric entity builder and set the TRS
  CellToFaceGEBuilder::GeoData& geoDataCell = m_cellBuilder->getDataGE();
  geoDataCell.trs = cells;
  
  // reset cell flags, these tell whether the cell flux part of the cell has been added yet
  for (CFuint iCell = 0; iCell < m_cellFlags.size(); ++iCell)
  {
    m_cellFlags[iCell] = false;
  }
  
  // get InnerFaces TopologicalRegionSet
  SafePtr<TopologicalRegionSet> faces = MeshDataStack::getActive()->getTrs("InnerFaces");

  // get the face start indexes
  vector< CFuint >& innerFacesStartIdxs = getMethodData().getInnerFacesStartIdxs();

  // get number of face orientations
  const CFuint nbrFaceOrients = innerFacesStartIdxs.size()-1;

  // get the geodata of the face builder and set the TRSs
  FaceToCellGEBuilder::GeoData& geoDataFace = m_faceBuilder->getDataGE();
  geoDataFace.cellsTRS = cells;
  geoDataFace.facesTRS = faces;
  geoDataFace.isBoundary = false;
  
  // get the geodata of the cell builders and set the TRS
  CellToFaceGEBuilder::GeoData& geoDataCBL = m_cellBuilders[LEFT]->getDataGE();
  geoDataCBL.trs = cells;
  CellToFaceGEBuilder::GeoData& geoDataCBR = m_cellBuilders[RIGHT]->getDataGE();
  geoDataCBR.trs = cells;
  
  ////////////////////COMPUTE GRADIENTS/////////////////////////
  
  //// Loop over faces to add the face part to the gradients
  
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
      geoDataFace.idx = faceID;
      m_face = m_faceBuilder->buildGE();

      // get the neighbouring cells
      m_cells[LEFT ] = m_face->getNeighborGeo(LEFT );
      m_cells[RIGHT] = m_face->getNeighborGeo(RIGHT);

      // get the states in the neighbouring cells
      m_states[LEFT ] = m_cells[LEFT ]->getStates();
      m_states[RIGHT] = m_cells[RIGHT]->getStates();
 
      // set the needed data
      setFaceDataForGradients(m_face->getID());
          
      // compute the face correction term of the corrected gradients
      computeGradientFaceCorrections();

      // release the cells
      m_cellBuilders[LEFT ]->releaseGE();
      m_cellBuilders[RIGHT]->releaseGE();
      
      // release the GeometricEntity
      m_faceBuilder->releaseGE();
    }
  }
  
  //// Loop over the elements to compute the cell
  //// part of gradients
  
  // loop over element types, for the moment there should only be one
  const CFuint nbrElemTypes = elemType->size();
  cf_assert(nbrElemTypes == 1);
  for (m_iElemType = 0; m_iElemType < nbrElemTypes; ++m_iElemType)
  {
    // get start and end indexes for this type of element
    const CFuint startIdx = (*elemType)[m_iElemType].getStartIdx();
    const CFuint endIdx   = (*elemType)[m_iElemType].getEndIdx();

    // loop over cells
    for (CFuint elemIdx = startIdx; elemIdx < endIdx; ++elemIdx)
    {
      // build the GeometricEntity
      geoDataCell.idx = elemIdx;
      m_elemIdx = elemIdx;
      m_cell = m_cellBuilder->buildGE();

      // get the states in this cell
      m_cellStates = m_cell->getStates();
      
      // add the cell part to the gradients
      computeGradients();
      
      //release the GeometricEntity
      m_cellBuilder->releaseGE();
    }
  }
  
  ////////////////////COMPUTE RHS////////////////////////
  
  // get the cell volumes
  DataHandle< CFreal > cellVolumes = socket_cellVolumes.getDataHandle();
  
  //// Loop over faces to compute the RHS
  
  // loop over different orientations
  for (m_orient = 0; m_orient < nbrFaceOrients; ++m_orient)
  {
    CFLog(VERBOSE, "Orient = " << m_orient << "\n");
    // start and stop index of the faces with this orientation
    const CFuint faceStartIdx = innerFacesStartIdxs[m_orient  ];
    const CFuint faceStopIdx  = innerFacesStartIdxs[m_orient+1];

    // loop over faces with this orientation
    for (CFuint faceID = faceStartIdx; faceID < faceStopIdx; ++faceID)
    {
      // build the face GeometricEntity
      geoDataFace.idx = faceID;
      m_face = m_faceBuilder->buildGE();

      // get the neighbouring cells
      m_cells[LEFT ] = m_face->getNeighborGeo(LEFT );
      m_cells[RIGHT] = m_face->getNeighborGeo(RIGHT);

      // get the states in the neighbouring cells
      m_states[LEFT ] = m_cells[LEFT ]->getStates();
      m_states[RIGHT] = m_cells[RIGHT]->getStates();
      
      // compute volume
      m_cellVolume[LEFT] = cellVolumes[m_cells[LEFT]->getID()];
      m_cellVolume[RIGHT] = cellVolumes[m_cells[RIGHT]->getID()];
      
      cf_assert(m_cellVolume[LEFT] > 0.0);
      cf_assert(m_cellVolume[RIGHT] > 0.0);
      
      // if one of the neighbouring cells is parallel updatable, compute the correction flux
      if ((*m_states[LEFT ])[0]->isParUpdatable() || (*m_states[RIGHT])[0]->isParUpdatable())
      {
	// build the neighbouring cells
        const CFuint cellIDL = m_face->getNeighborGeo(LEFT)->getID();
        geoDataCBL.idx = cellIDL;
        m_cells[LEFT] = m_cellBuilders[LEFT ]->buildGE();
        const CFuint cellIDR = m_face->getNeighborGeo(RIGHT)->getID();
        geoDataCBR.idx = cellIDR;
        m_cells[RIGHT] = m_cellBuilders[RIGHT]->buildGE();

	// set the face data
	setFaceData(m_face->getID());//faceID

	// compute the left and right states and gradients in the flx pnts
	computeFlxPntStatesAndGrads();

	// compute FI
	computeInterfaceFlxCorrection();

	// compute the wave speed updates
        computeWaveSpeedUpdates(m_waveSpeedUpd);

        // update the wave speed
        updateWaveSpeed();

	// compute the correction for the left neighbour
	computeCorrection(LEFT, m_divContFlxL);
	m_divContFlx = m_divContFlxL;

	// update RHS
	updateRHS();
	
	// compute the correction for the right neighbour
	computeCorrection(RIGHT, m_divContFlxR);
	m_divContFlx = m_divContFlxR;
	
	// update RHS
	updateRHS();
        
        // compute needed cell contributions: what used to be cell loop is incorporated here!!
        if (!m_cellFlags[cellIDL] && (*m_states[LEFT ])[0]->isParUpdatable())
        {
          // compute cell contribution 
          computeUnpertCellDiffResiduals(LEFT);
          
	  // update RHS
	  updateRHSUnpertCell(LEFT);
        }
        if (!m_cellFlags[cellIDR] && (*m_states[RIGHT])[0]->isParUpdatable())
        {
          // compute cell contribution 
          computeUnpertCellDiffResiduals(RIGHT);

	  // update RHS
	  updateRHSUnpertCell(RIGHT);
        }

	// get all the faces neighbouring the cells
        m_faces[LEFT ] = m_cells[LEFT ]->getNeighborGeos();
        m_faces[RIGHT] = m_cells[RIGHT]->getNeighborGeos();

        // set the local indexes of the other faces than the current faces
        setOtherFacesLocalIdxs();

        // get the sol pnt jacobians
        DataHandle< CFreal > volumes = socket_volumes.getDataHandle();
        
        // compute solution points Jacobian determinants
	for (CFuint iSide = 0; iSide < 2; ++iSide)
        {
          for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
          {
            m_solJacobDet[iSide][iSol] = volumes[(*(m_states[iSide]))[iSol]->getLocalID()];
          }
	}
        
        const CFuint iter = SubSystemStatusStack::getActive()->getNbIter();
    
        const CFuint iterFreeze = getMethodData().getFreezeJacobIter();
    
        const CFuint interval = iter - iterFreeze;
        
        if (!getMethodData().freezeJacob() || iter < iterFreeze || interval % getMethodData().getFreezeJacobInterval() == 0)
        {

        // compute the diffusive face term contribution to the jacobian
        if ((*m_states[LEFT])[0]->isParUpdatable() && (*m_states[RIGHT])[0]->isParUpdatable())
        {
          computeBothJacobsDiffFaceTerm();
        }
        else if ((*m_states[LEFT])[0]->isParUpdatable())
        {
          computeOneJacobDiffFaceTerm(LEFT );
        }
        else if ((*m_states[RIGHT])[0]->isParUpdatable())
        {
          computeOneJacobDiffFaceTerm(RIGHT);
        }
        
        }

        // release the cells
        m_cellBuilders[LEFT ]->releaseGE();
        m_cellBuilders[RIGHT]->releaseGE();
        
        m_cellFlags[cellIDL] = true;
        m_cellFlags[cellIDR] = true;
      }
      
      // release the GeometricEntity
      m_faceBuilder->releaseGE();
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffJacobFluxReconstruction::setFaceDataForGradients(CFuint faceID)
{   
  // get the face flux point normals
  DataHandle< CFreal > flxPntNormals = socket_flxPntNormals.getDataHandle();
  
  // Loop over flux points to set the normal vectors
  for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
  { 
    for (CFuint iDim = 0; iDim < m_dim; ++iDim)
    {
      m_faceJacobVecs[iFlxPnt][iDim] = flxPntNormals[m_face->getID()*m_nbrFaceFlxPnts*m_dim+iFlxPnt*m_dim+iDim];
    }
  }

  // Loop over flux points to set the normal vectors
  for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
  {
    // get face Jacobian vector sizes in the flux points
    DataHandle< vector< CFreal > > faceJacobVecSizeFaceFlxPnts = socket_faceJacobVecSizeFaceFlxPnts.getDataHandle();
  
    // get face Jacobian vector size
    m_faceJacobVecAbsSizeFlxPnts[iFlxPnt] = faceJacobVecSizeFaceFlxPnts[faceID][iFlxPnt];

    // set face Jacobian vector size with sign depending on mapped coordinate direction
    m_faceJacobVecSizeFlxPnts[iFlxPnt][LEFT] = m_faceJacobVecAbsSizeFlxPnts[iFlxPnt]*(*m_faceMappedCoordDir)[m_orient][LEFT];
    m_faceJacobVecSizeFlxPnts[iFlxPnt][RIGHT] = m_faceJacobVecAbsSizeFlxPnts[iFlxPnt]*(*m_faceMappedCoordDir)[m_orient][RIGHT];

    // set unit normal vector
    m_unitNormalFlxPnts[iFlxPnt] = m_faceJacobVecs[iFlxPnt]/m_faceJacobVecAbsSizeFlxPnts[iFlxPnt];
  }

  // loop over flx pnts to extrapolate the states to the flux points
  for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
  {     
    // local flux point indices in the left and right cell
    const CFuint flxPntIdxL = (*m_faceFlxPntConnPerOrient)[m_orient][LEFT][iFlxPnt];
    const CFuint flxPntIdxR = (*m_faceFlxPntConnPerOrient)[m_orient][RIGHT][iFlxPnt];
    
    // reset states in flx pnt
    *(m_cellStatesFlxPnt[LEFT][iFlxPnt]) = 0.0;
    *(m_cellStatesFlxPnt[RIGHT][iFlxPnt]) = 0.0;

    // extrapolate the left and right states to the flx pnts
    for (CFuint iSol = 0; iSol < m_nbrSolDep; ++iSol)
    {
      const CFuint solIdxL = (*m_flxSolDep)[flxPntIdxL][iSol];
      const CFuint solIdxR = (*m_flxSolDep)[flxPntIdxR][iSol];
 
      // add the contributions of the current sol pnt
      *(m_cellStatesFlxPnt[LEFT][iFlxPnt]) += (*m_solPolyValsAtFlxPnts)[flxPntIdxL][solIdxL]*(*((*(m_states[LEFT]))[solIdxL]));
      *(m_cellStatesFlxPnt[RIGHT][iFlxPnt]) += (*m_solPolyValsAtFlxPnts)[flxPntIdxR][solIdxR]*(*((*(m_states[RIGHT]))[solIdxR]));
    }
  }
}


//////////////////////////////////////////////////////////////////////////////

void ConvDiffJacobFluxReconstruction::setFaceData(CFuint faceID)
{   
  // get the face flux point normals
  DataHandle< CFreal > flxPntNormals = socket_flxPntNormals.getDataHandle();
  
  // Loop over flux points to set the normal vectors
  for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
  { 
    for (CFuint iDim = 0; iDim < m_dim; ++iDim)
    {
      m_faceJacobVecs[iFlxPnt][iDim] = flxPntNormals[m_face->getID()*m_nbrFaceFlxPnts*m_dim+iFlxPnt*m_dim+iDim];
    }
  }

  // Loop over flux points to set the normal vectors
  for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
  {
    // get face Jacobian vector sizes in the flux points
    DataHandle< vector< CFreal > > faceJacobVecSizeFaceFlxPnts = socket_faceJacobVecSizeFaceFlxPnts.getDataHandle();
  
    // get face Jacobian vector size
    m_faceJacobVecAbsSizeFlxPnts[iFlxPnt] = faceJacobVecSizeFaceFlxPnts[faceID][iFlxPnt];

    // set face Jacobian vector size with sign depending on mapped coordinate direction
    m_faceJacobVecSizeFlxPnts[iFlxPnt][LEFT] = m_faceJacobVecAbsSizeFlxPnts[iFlxPnt]*(*m_faceMappedCoordDir)[m_orient][LEFT];
    m_faceJacobVecSizeFlxPnts[iFlxPnt][RIGHT] = m_faceJacobVecAbsSizeFlxPnts[iFlxPnt]*(*m_faceMappedCoordDir)[m_orient][RIGHT];

    // set unit normal vector
    m_unitNormalFlxPnts[iFlxPnt] = m_faceJacobVecs[iFlxPnt]/m_faceJacobVecAbsSizeFlxPnts[iFlxPnt];
  }

  // compute inverse characteristic lengths
  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    m_faceInvCharLengths[iFlx] = 2.*m_faceJacobVecAbsSizeFlxPnts[iFlx]/(m_cellVolume[LEFT] + m_cellVolume[RIGHT]);
  }
  
  // get the gradients datahandle
  DataHandle< vector< RealVector > > gradients = socket_gradients.getDataHandle();

  for (CFuint iSide = 0; iSide < 2; ++iSide)
  {
    for (CFuint iState = 0; iState < m_nbrSolPnts; ++iState)
    {
      const CFuint stateID = (*(m_states[iSide]))[iState]->getLocalID();
      m_cellGrads[iSide][iState] = &gradients[stateID];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffJacobFluxReconstruction::computeFlxPntStatesAndGrads()
{
  // loop over flx pnts to extrapolate the states to the flux points
  for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
  {     
    // local flux point indices in the left and right cell
    const CFuint flxPntIdxL = (*m_faceFlxPntConnPerOrient)[m_orient][LEFT][iFlxPnt];
    const CFuint flxPntIdxR = (*m_faceFlxPntConnPerOrient)[m_orient][RIGHT][iFlxPnt];
    
    // reset states in flx pnt
    *(m_cellStatesFlxPnt[LEFT][iFlxPnt]) = 0.0;
    *(m_cellStatesFlxPnt[RIGHT][iFlxPnt]) = 0.0;

    // reset the grads in the flx pnts
    for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
    {
      *(m_cellGradFlxPnt[LEFT][iFlxPnt][iVar]) = 0.0;
      *(m_cellGradFlxPnt[RIGHT][iFlxPnt][iVar]) = 0.0;
    }

    // extrapolate the left and right states to the flx pnts
    for (CFuint iSol = 0; iSol < m_nbrSolDep; ++iSol)
    {
      const CFuint solIdxL = (*m_flxSolDep)[flxPntIdxL][iSol];
      const CFuint solIdxR = (*m_flxSolDep)[flxPntIdxR][iSol];
 
      // add the contributions of the current sol pnt
      *(m_cellStatesFlxPnt[LEFT][iFlxPnt]) += (*m_solPolyValsAtFlxPnts)[flxPntIdxL][solIdxL]*(*((*(m_states[LEFT]))[solIdxL]));
      *(m_cellStatesFlxPnt[RIGHT][iFlxPnt]) += (*m_solPolyValsAtFlxPnts)[flxPntIdxR][solIdxR]*(*((*(m_states[RIGHT]))[solIdxR]));

      for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
      {
        *(m_cellGradFlxPnt[LEFT][iFlxPnt][iVar]) += (*m_solPolyValsAtFlxPnts)[flxPntIdxL][solIdxL]*((*(m_cellGrads[LEFT][solIdxL]))[iVar]);
	*(m_cellGradFlxPnt[RIGHT][iFlxPnt][iVar]) += (*m_solPolyValsAtFlxPnts)[flxPntIdxR][solIdxR]*((*(m_cellGrads[RIGHT][solIdxR]))[iVar]);
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffJacobFluxReconstruction::computeInterfaceFlxCorrection()
{
  // Loop over the flux points to calculate FI
  for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
  { 
    // compute the average sol and grad to use the BR2 scheme
    for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
    {
      *(m_avgGrad[iVar]) = (*(m_cellGradFlxPnt[LEFT][iFlxPnt][iVar]) + *(m_cellGradFlxPnt[RIGHT][iFlxPnt][iVar]))/2.0;
                   
      m_avgSol[iVar] = ((*(m_cellStatesFlxPnt[LEFT][iFlxPnt]))[iVar] + (*(m_cellStatesFlxPnt[RIGHT][iFlxPnt]))[iVar])/2.0; 
    }
    
    prepareFluxComputation();
     
    // compute diffusive flux
    computeFlux(m_avgSol,m_avgGrad,m_unitNormalFlxPnts[iFlxPnt],0,m_flxPntRiemannFluxDiff[iFlxPnt]);
    
    m_flxPntRiemannFlux[iFlxPnt] = m_flxPntRiemannFluxDiff[iFlxPnt];
    
    // compute the convective riemann flux
    m_flxPntRiemannFlux[iFlxPnt] -= m_riemannFluxComputer->computeFlux(*(m_cellStatesFlxPnt[LEFT][iFlxPnt]),
									    *(m_cellStatesFlxPnt[RIGHT][iFlxPnt]),
									    m_unitNormalFlxPnts[iFlxPnt]);
     
    // compute FI in the mapped coord frame
    m_cellFlx[LEFT][iFlxPnt] = (m_flxPntRiemannFlux[iFlxPnt])*m_faceJacobVecSizeFlxPnts[iFlxPnt][LEFT];
    m_cellFlx[RIGHT][iFlxPnt] = (m_flxPntRiemannFlux[iFlxPnt])*m_faceJacobVecSizeFlxPnts[iFlxPnt][RIGHT];
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffJacobFluxReconstruction::computeWaveSpeedUpdates(vector< CFreal >& waveSpeedUpd)
{
  // compute the wave speed updates for the neighbouring cells
  cf_assert(waveSpeedUpd.size() == 2);
  
  // add convective part
  for (CFuint iSide = 0; iSide < 2; ++iSide)
  {
    waveSpeedUpd[iSide] = 0.0;
  
    for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
    {
      const CFreal jacobXIntCoef = m_faceJacobVecAbsSizeFlxPnts[iFlx]*(*m_faceIntegrationCoefs)[iFlx];
      
      // transform update states to physical data to calculate eigenvalues
      m_updateVarSet->computePhysicalData(*(m_cellStatesFlxPnt[iSide][iFlx]), m_pData);
      
      waveSpeedUpd[iSide] += jacobXIntCoef * m_updateVarSet->getMaxAbsEigenValue(m_pData,m_unitNormalFlxPnts[iFlx]);
    }
  }
  
  //////Diffusive part is added in NS!!!!
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffJacobFluxReconstruction::initJacobianComputation()
{
  CFLog(VERBOSE, "initJacobianComputation\n");
    
  // get the sol pnt normals
  DataHandle< CFreal > solPntNormals = socket_solPntNormals.getDataHandle();
  
  // set block row and column indices, proj vectors and make a backup of discontinuous fluxes
  for (m_pertSide = 0; m_pertSide < 2; ++m_pertSide)
  {
    // store the sol pnt normals
    for (CFuint iState = 0; iState < m_nbrSolPnts; ++iState)
    {
      const CFuint solID = (*(m_states[m_pertSide]))[iState]->getLocalID();
      
      for (CFuint iDim = 0; iDim < m_dim+m_ndimplus; ++iDim)
      {
        for (CFuint jDim = 0; jDim < m_dim; ++jDim)
        {
          m_neighbCellFluxProjVects[m_pertSide][iDim][iState][jDim] = solPntNormals[solID*(m_dim+m_ndimplus)*m_dim+iDim*m_dim+jDim];
        }
      }
    }
    
    // Loop over solution points to calculate the discontinuous flux.
    for (m_pertSol = 0; m_pertSol < m_nbrSolPnts; ++m_pertSol)
    { 
      for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
      {
        *(m_tempGrad[iVar]) = (*(m_cellGrads[m_pertSide][m_pertSol]))[iVar];
      }

      m_avgSol = *((*(m_states[m_pertSide]))[m_pertSol]->getData());
      
      m_updateVarSet->computePhysicalData(*((*(m_states[m_pertSide]))[m_pertSol]), m_pData); 

      prepareFluxComputation();

      // calculate the discontinuous flux projected on x, y, z-directions
      for (CFuint iDim = 0; iDim < m_dim+m_ndimplus; ++iDim)
      {
        // diffusive part 
        computeFlux(m_avgSol,m_tempGrad,m_neighbCellFluxProjVects[m_pertSide][iDim][m_pertSol],0,m_contFlxBackupDiff[m_pertSide][m_pertSol][iDim]);
        
        m_contFlxBackup[m_pertSide][m_pertSol][iDim] = m_contFlxBackupDiff[m_pertSide][m_pertSol][iDim];
        
        // convective part
        m_contFlxBackup[m_pertSide][m_pertSol][iDim] -= m_updateVarSet->getFlux()(m_pData,m_neighbCellFluxProjVects[m_pertSide][iDim][m_pertSol]);
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffJacobFluxReconstruction::computeCellFluxJacobianNum(const CFreal resFactor)
{
  CFLog(VERBOSE, "computeCellFluxJacobianNum\n");
    
  // loop over left and right cell to compute flux jacobian
  for (m_pertSide = 0; m_pertSide < 2; ++m_pertSide)
  {
    // variable for the other side
    const CFuint iOtherSide = m_pertSide == LEFT ? RIGHT : LEFT;
    
    // cell ID of the cell at the non-perturbed side
    const CFuint otherCellID = m_cells[iOtherSide]->getID();

    // loop over the states to perturb the states
    for (m_pertSol = 0; m_pertSol < m_nbrSolPnts; ++m_pertSol)
    {
      // dereference state
      State& pertState = *(*m_states[m_pertSide])[m_pertSol];
      
      for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
      {
        *(m_tempGrad[iVar]) = (*(m_cellGrads[m_pertSide][m_pertSol]))[iVar];
      }

      // loop over the variables in the state
      for (m_pertVar = 0; m_pertVar < m_nbrEqs; ++m_pertVar)
      {
        // perturb physical variable in state
        m_numJacob->perturb(m_pertVar,pertState[m_pertVar]);
        
        // compute perturbed physical data
        m_updateVarSet->computePhysicalData(*((*(m_states[m_pertSide]))[m_pertSol]), m_pData); 

        m_avgSol = *((*(m_states[m_pertSide]))[m_pertSol]->getData());
        
        // compute perturbed fluxes
        prepareFluxComputation();

        // calculate the discontinuous flux projected on x, y, z-directions
        for (CFuint iDim = 0; iDim < m_dim+m_ndimplus; ++iDim)
        {       
          // diffusive part
          computeFlux(m_avgSol,m_tempGrad,m_neighbCellFluxProjVects[m_pertSide][iDim][m_pertSol],0,m_contFlxNeighb[m_pertSide][m_pertSol][iDim]);
                    
          // convective part
          m_contFlxNeighb[m_pertSide][m_pertSol][iDim] -= m_updateVarSet->getFlux()(m_pData,m_neighbCellFluxProjVects[m_pertSide][iDim][m_pertSol]);

          // compute the flux current jacobian term
          // compute the finite difference derivative of the face term (note an implicit minus sign is added here, by the ordering of arguments)
          m_numJacob->computeDerivative(m_contFlxNeighb[m_pertSide][m_pertSol][iDim],m_contFlxBackup[m_pertSide][m_pertSol][iDim],m_tempFlux);

          // multiply residual update derivatives with residual factor so it is taken into the final jacobian
          m_tempFlux *= resFactor;
        
          // store the flux jacobian
          m_fluxJacobian[m_pertSide][m_pertSol][m_pertVar][iDim] = m_tempFlux;
        }

        // restore physical variable in state
        m_numJacob->restore(pertState[m_pertVar]);
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffJacobFluxReconstruction::computeRiemannFluxJacobianNum(const CFreal resFactor)
{
  CFLog(VERBOSE, "computeRiemannFluxJacobianNum\n");
    
  // loop over the face flux points to compute the Riemann flux jacobian
  for (m_pertSide = 0; m_pertSide < 2; ++m_pertSide)
  {
    // variable for the other side
    const CFuint iOtherSide = m_pertSide == LEFT ? RIGHT : LEFT;
    
    for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
    {
      // dereference state
      State& pertState = *(m_cellStatesFlxPnt[m_pertSide][iFlxPnt]);
      
      // compute the average grad to use the BR2 scheme
      for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
      {
        *(m_avgGrad[iVar]) = (*(m_cellGradFlxPnt[LEFT][iFlxPnt][iVar]) + *(m_cellGradFlxPnt[RIGHT][iFlxPnt][iVar]))/2.0;            
      }
                
      // loop over the variables in the state
      for (m_pertVar = 0; m_pertVar < m_nbrEqs; ++m_pertVar)
      {
        // perturb physical variable in state
        m_numJacob->perturb(m_pertVar,pertState[m_pertVar]);
        
        // compute the average sol
        for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
        {        
          m_avgSol[iVar] = (pertState[iVar] + (*(m_cellStatesFlxPnt[iOtherSide][iFlxPnt]))[iVar])/2.0; 
        }
    
        prepareFluxComputation();
     
        // compute diffusive flux
        computeFlux(m_avgSol,m_avgGrad,m_unitNormalFlxPnts[iFlxPnt],0,m_flxPntRiemannFluxPert[iFlxPnt]);
        
        if (m_pertSide == LEFT)
        {
          // compute the convective riemann flux
          m_flxPntRiemannFluxPert[iFlxPnt] -= m_riemannFluxComputer->computeFlux(pertState,
									              *(m_cellStatesFlxPnt[RIGHT][iFlxPnt]),
									              m_unitNormalFlxPnts[iFlxPnt]);
        }
        else
        {
          // compute the convective riemann flux
          m_flxPntRiemannFluxPert[iFlxPnt] -= m_riemannFluxComputer->computeFlux(*(m_cellStatesFlxPnt[LEFT][iFlxPnt]),
                                                                                      pertState,
									              m_unitNormalFlxPnts[iFlxPnt]);
        }
       
        // compute the flux current jacobian term
        // compute the finite difference derivative of the face term (note an implicit minus sign is added here, by the ordering of arguments)
        m_numJacob->computeDerivative(m_flxPntRiemannFluxPert[iFlxPnt],m_flxPntRiemannFlux[iFlxPnt],m_riemannFluxJacobian[m_pertSide][iFlxPnt][m_pertVar]);

        // multiply residual update derivatives with residual factor so it is taken into the final jacobian
        m_riemannFluxJacobian[m_pertSide][iFlxPnt][m_pertVar] *= resFactor;

        // restore physical variable in state
        m_numJacob->restore(pertState[m_pertVar]);
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffJacobFluxReconstruction::computeFluxToGradJacobianNum(const CFreal resFactor)
{
  CFLog(VERBOSE, "computeFluxToGradJacobianNum\n");
    
  // loop over left and right cell to compute flux jacobian
  for (m_pertSide = 0; m_pertSide < 2; ++m_pertSide)
  {
    // variable for the other side
    const CFuint iOtherSide = m_pertSide == LEFT ? RIGHT : LEFT;
    
    // cell ID of the cell at the non-perturbed side
    const CFuint otherCellID = m_cells[iOtherSide]->getID();

    // loop over the states to perturb the states
    for (m_pertSol = 0; m_pertSol < m_nbrSolPnts; ++m_pertSol)
    {
      // dereference state
      m_avgSol = *((*(m_states[m_pertSide]))[m_pertSol]->getData());
      
      for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
      {
        *(m_tempGrad[iVar]) = (*(m_cellGrads[m_pertSide][m_pertSol]))[iVar];
      }

      // loop over the variables in the state
      for (m_pertVar = 0; m_pertVar < m_nbrEqs; ++m_pertVar)
      {
        for (CFuint pertDim = 0; pertDim < m_dim; ++pertDim)
        {
          // perturb physical variable in state
          m_numJacob->perturb(m_pertVar,(*(m_tempGrad[m_pertVar]))[pertDim]);

          // compute perturbed fluxes
          prepareFluxComputation();

          // calculate the discontinuous flux projected on x, y, z-directions
          for (CFuint iDim = 0; iDim < m_dim+m_ndimplus; ++iDim)
          {       
            // diffusive part
            computeFlux(m_avgSol,m_tempGrad,m_neighbCellFluxProjVects[m_pertSide][iDim][m_pertSol],0,m_contFlxNeighb[m_pertSide][m_pertSol][iDim]);

            // compute the flux current jacobian term
            // compute the finite difference derivative of the face term (note an implicit minus sign is added here, by the ordering of arguments)
            m_numJacob->computeDerivative(m_contFlxNeighb[m_pertSide][m_pertSol][iDim],m_contFlxBackupDiff[m_pertSide][m_pertSol][iDim],m_tempFlux);

            // multiply residual update derivatives with residual factor so it is taken into the final jacobian
            m_tempFlux *= resFactor;

            // store the flux jacobian
            m_gradientFluxJacobian[m_pertSide][m_pertSol][m_pertVar][pertDim][iDim] = m_tempFlux;
          }

          // restore physical variable in state
          m_numJacob->restore((*(m_tempGrad[m_pertVar]))[pertDim]);
        }
      }
    }
  }  
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffJacobFluxReconstruction::computeRiemannFluxToGradJacobianNum(const CFreal resFactor)
{  
  CFLog(VERBOSE, "computeRiemannFluxToGradJacobianNum\n");
    
  for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
  {      
    // compute the average sol
    for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
    {        
      m_avgSol[iVar] = ((*(m_cellStatesFlxPnt[LEFT][iFlxPnt]))[iVar] + (*(m_cellStatesFlxPnt[RIGHT][iFlxPnt]))[iVar])/2.0; 
    }
    
    // compute the average grad to use the BR2 scheme
    for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
    {
      *(m_avgGrad[iVar]) = (*(m_cellGradFlxPnt[LEFT][iFlxPnt][iVar]) + *(m_cellGradFlxPnt[RIGHT][iFlxPnt][iVar]))/2.0;            
    }
           
    // loop over the variables in the state
    for (m_pertVar = 0; m_pertVar < m_nbrEqs; ++m_pertVar)
    {
      for (CFuint pertDir = 0; pertDir < m_dim; ++pertDir)
      {
        // perturb physical variable in state
        m_numJacob->perturb(m_pertVar,(*(m_avgGrad[m_pertVar]))[pertDir]);

        prepareFluxComputation();
     
        // compute diffusive flux
        computeFlux(m_avgSol,m_avgGrad,m_unitNormalFlxPnts[iFlxPnt],0,m_flxPntRiemannFluxPert[iFlxPnt]);
     
        // compute the flux current jacobian term
        // compute the finite difference derivative of the face term (note an implicit minus sign is added here, by the ordering of arguments)
        m_numJacob->computeDerivative(m_flxPntRiemannFluxPert[iFlxPnt],m_flxPntRiemannFluxDiff[iFlxPnt],m_riemannFluxGradJacobian[iFlxPnt][m_pertVar][pertDir]);

        // multiply residual update derivatives with residual factor so it is taken into the final jacobian
        m_riemannFluxGradJacobian[iFlxPnt][m_pertVar][pertDir] *= resFactor;

        // restore physical variable in state
        m_numJacob->restore((*(m_avgGrad[m_pertVar]))[pertDir]);
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffJacobFluxReconstruction::computeGradToStateJacobianAna()
{
  CFLog(VERBOSE, "computeGradToStateJacobianAna\n");
  
  // reset the jacobian
  for (m_pertSide = 0; m_pertSide < 2; ++m_pertSide)
  { 
    for (m_pertSol = 0; m_pertSol < m_nbrSolPnts; ++m_pertSol)
    {
      for (CFuint jSol = 0; jSol < m_nbrSolPnts; ++jSol)
      {
        for (CFuint iDim = 0; iDim < m_dim; ++iDim)
        { 
          m_gradientStateJacobian[m_pertSide][m_pertSol][LEFT][jSol][iDim] = 0.0;
          m_gradientStateJacobian[m_pertSide][m_pertSol][RIGHT][jSol][iDim] = 0.0;
        }
      }
    }
  }
  
  // get the face flux point normals
  DataHandle< CFreal > flxPntNormals = socket_flxPntNormals.getDataHandle();
    
  for (m_pertSide = 0; m_pertSide < 2; ++m_pertSide)
  {
    // variable for the other side
    const CFuint iOtherSide = m_pertSide == LEFT ? RIGHT : LEFT;
    
    // loop over the sol pnts of which the gradients will be derived
    for (m_pertSol = 0; m_pertSol < m_nbrSolPnts; ++m_pertSol)
    {
      // loop over the sol pnts to which to derive
      for (CFuint jSol = 0; jSol < m_nbrSolSolDep; ++jSol)
      {
        const CFuint jSolIdx = (*m_solSolDep)[m_pertSol][jSol];  
        
        for (CFuint iDim = 0; iDim < m_dim; ++iDim)
        {                    
          for (CFuint jDim = 0; jDim < m_dim; ++jDim)
          {
            m_gradientStateJacobian[m_pertSide][jSolIdx][m_pertSide][m_pertSol][iDim] += (m_neighbCellFluxProjVects[m_pertSide][jDim+m_ndimplus][m_pertSol][iDim]) * (*m_solPolyDerivAtSolPnts)[jSolIdx][jDim][m_pertSol];
          }
        }
      }
    }
    
    for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
    {
      // local flux point indices in the left and right cell
      const CFuint flxPntIdxThis = (*m_faceFlxPntConnPerOrient)[m_orient][m_pertSide][iFlx];
      const CFuint flxPntIdxOther = (*m_faceFlxPntConnPerOrient)[m_orient][iOtherSide][iFlx];
      
      for (m_pertSol = 0; m_pertSol < m_nbrSolDep; ++m_pertSol)
      {
        const CFuint pertSolIdx = (*m_flxSolDep)[flxPntIdxThis][m_pertSol]; 
    
        for (CFuint jSol = 0; jSol < m_nbrSolDep; ++jSol)
        {
          const CFuint jSolIdx = (*m_flxSolDep)[flxPntIdxThis][jSol];  
          const CFuint jSolIdxOther = (*m_flxSolDep)[flxPntIdxOther][jSol];  
        
          for (CFuint iDim = 0; iDim < m_dim; ++iDim)
          {
            m_gradientStateJacobian[m_pertSide][jSolIdx][m_pertSide][pertSolIdx][iDim] -= 0.5 * m_corrFctDiv[jSolIdx][flxPntIdxThis] * (*m_solPolyValsAtFlxPnts)[flxPntIdxThis][pertSolIdx] *
                    (*m_faceMappedCoordDir)[m_orient][m_pertSide]*m_faceJacobVecs[iFlx][iDim];
         
            m_gradientStateJacobian[iOtherSide][jSolIdxOther][m_pertSide][pertSolIdx][iDim] += 0.5 * m_corrFctDiv[jSolIdxOther][flxPntIdxOther] * (*m_solPolyValsAtFlxPnts)[flxPntIdxThis][pertSolIdx] *
                    (*m_faceMappedCoordDir)[m_orient][iOtherSide]*m_faceJacobVecs[iFlx][iDim];
          }
        }   
      }
    }
  
    // Add the contribution of the correction to the gradients for each face
    // compute other face contributions to the gradients
    const CFuint nbrOtherFaces = m_otherFaceLocalIdxs[m_pertSide].size();
  
    for (CFuint iFace = 0; iFace < nbrOtherFaces; ++iFace)
    {
      // get local face index
      const CFuint faceIdx = m_otherFaceLocalIdxs[m_pertSide][iFace];
      
      const CFuint faceID = (*m_faces[m_pertSide])[faceIdx]->getID();
      
      if ((*m_isFaceOnBoundary[m_pertSide])[faceIdx])
      {  
        for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
        {
          const CFuint currFlxIdx = (*m_faceFlxPntConn)[faceIdx][iFlxPnt];
          
          for (m_pertSol = 0; m_pertSol < m_nbrSolDep; ++m_pertSol)
          {
            const CFuint pertSolIdx = (*m_flxSolDep)[currFlxIdx][m_pertSol];        
        
            for (CFuint jSol = 0; jSol < m_nbrSolDep; ++jSol)
            {
              const CFuint jSolIdx = (*m_flxSolDep)[currFlxIdx][jSol]; 
    
              for (CFuint iDim = 0; iDim < m_dim; ++iDim)
              {
                m_gradientStateJacobian[m_pertSide][jSolIdx][m_pertSide][pertSolIdx][iDim] -= 0.5 * m_corrFctDiv[jSolIdx][currFlxIdx] * (*m_solPolyValsAtFlxPnts)[currFlxIdx][pertSolIdx] *
                          (*m_faceLocalDir)[faceIdx]*flxPntNormals[faceID*m_nbrFaceFlxPnts*m_dim+iFlxPnt*m_dim+iDim];
              }
            }
          }
        }
      }
      else
      {
        // Get orientation of face
        const CFuint orient = (*m_faceOrients[m_pertSide])[faceIdx];
        
        // cell side with respect to this face
        const CFuint cellSide = (*m_currCellSide[m_pertSide])[faceIdx];
      
        for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
        {
          const CFuint currFlxIdx = (*m_faceFlxPntConnPerOrient)[orient][cellSide][iFlxPnt];
          
          for (m_pertSol = 0; m_pertSol < m_nbrSolDep; ++m_pertSol)
          {
            const CFuint pertSolIdx = (*m_flxSolDep)[currFlxIdx][m_pertSol];        
        
            for (CFuint jSol = 0; jSol < m_nbrSolDep; ++jSol)
            {
              const CFuint jSolIdx = (*m_flxSolDep)[currFlxIdx][jSol]; 
    
              for (CFuint iDim = 0; iDim < m_dim; ++iDim)
              {
                m_gradientStateJacobian[m_pertSide][jSolIdx][m_pertSide][pertSolIdx][iDim] -= 0.5 * m_corrFctDiv[jSolIdx][currFlxIdx] * (*m_solPolyValsAtFlxPnts)[currFlxIdx][pertSolIdx] *
                          ((*m_faceMappedCoordDir)[orient][cellSide])*flxPntNormals[faceID*m_nbrFaceFlxPnts*m_dim+iFlxPnt*m_dim+iDim];
              }
            }
          }
        } 
      }
    }
  }
  
  for (m_pertSide = 0; m_pertSide < 2; ++m_pertSide)
  { 
    for (m_pertSol = 0; m_pertSol < m_nbrSolPnts; ++m_pertSol)
    {
      
      
      for (CFuint jSol = 0; jSol < m_nbrSolPnts; ++jSol)
      {
        const CFreal invJacobL = 1./m_solJacobDet[LEFT][jSol];
        const CFreal invJacobR = 1./m_solJacobDet[RIGHT][jSol];
          
        for (CFuint iDim = 0; iDim < m_dim; ++iDim)
        { 
          m_gradientStateJacobian[LEFT][jSol][m_pertSide][m_pertSol][iDim] *= invJacobL;
          m_gradientStateJacobian[RIGHT][jSol][m_pertSide][m_pertSol][iDim] *= invJacobR;
        }
//        if (m_cells[LEFT]->getID()==1) CFLog(INFO,"side: " << m_pertSide << ", sol: " << m_pertSol << ", to side: 0, sol: " << jSol << ": " << m_gradientStateJacobian[LEFT][jSol][m_pertSide][m_pertSol] <<"\n");
//      if (m_cells[LEFT]->getID()==1) CFLog(INFO,"side: " << m_pertSide << ", sol: " << m_pertSol << ", to side: 1, sol: " << jSol << ": " << m_gradientStateJacobian[RIGHT][jSol][m_pertSide][m_pertSol] <<"\n");
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffJacobFluxReconstruction::computeGradVarsToStateJacobianNum()
{
  CFLog(VERBOSE, "computeGradVarsToStateJacobianNum\n");

  for (m_pertSide = 0; m_pertSide < 2; ++m_pertSide)
  { 
    for (m_pertSol = 0; m_pertSol < m_nbrSolPnts; ++m_pertSol)
    {
      for (m_pertVar = 0; m_pertVar < m_nbrEqs; ++m_pertVar)
      {
        for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
        {
          m_gradVarsToStateJacobian[m_pertSide][m_pertSol][m_pertVar][iEq] = 0.0;
        }
        m_gradVarsToStateJacobian[m_pertSide][m_pertSol][m_pertVar][m_pertVar] = 1.0;
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffJacobFluxReconstruction::computeBothJacobsDiffFaceTerm()
{
  CFLog(VERBOSE, "computeBothJacobsDiffFaceTerm\n");
    
  // get residual factor
  const CFreal resFactor = getMethodData().getResFactor();

  // dereference accumulator
  BlockAccumulator& acc = *m_acc;
  
  CFuint solIdx = 0;
  for (m_pertSide = 0; m_pertSide < 2; ++m_pertSide)
  {
    for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol, ++solIdx)
    {
      acc.setRowColIndex(solIdx,(*m_states[m_pertSide])[iSol]->getLocalID());
    }
  }

  //// compute the needed flux jacobians
  
  initJacobianComputation();
  
  computeCellFluxJacobianNum(resFactor);
  
  computeRiemannFluxJacobianNum(resFactor);
  
  computeFluxToGradJacobianNum(resFactor);
  
  if (m_addRiemannToGradJacob || m_addRiemannToGradCrossCellJacob) computeRiemannFluxToGradJacobianNum(resFactor);
  
  computeGradToStateJacobianAna();
  
  computeGradVarsToStateJacobianNum();
  
  //// add the total jacobians to the system jacobian
  
  // loop over left and right cell to add the discontinuous (cell) part to the jacobian
  for (m_pertSide = 0; m_pertSide < 2; ++m_pertSide)
  {
    // make sure this is only done once per cell
    if (!m_cellFlags[m_cells[m_pertSide]->getID()]) 
    {
      // variable for the other side
      const CFuint iOtherSide = m_pertSide == LEFT ? RIGHT : LEFT;

      // cell ID of the cell at the non-perturbed side
      const CFuint otherCellID = m_cells[iOtherSide]->getID();
    
      // term depending on iSide
      const CFuint pertSideTerm = m_pertSide*m_nbrSolPnts;

      // term depending on iOtherSide
      const CFuint otherSideTerm = iOtherSide*m_nbrSolPnts;

      // loop over the states to which to derive (l)
      for (m_pertSol = 0; m_pertSol < m_nbrSolPnts; ++m_pertSol)
      {
        // loop over the variables in the state (k)
        for (m_pertVar = 0; m_pertVar < m_nbrEqs; ++m_pertVar)
        { 
          // add the discontinuous part of the jacobian related to the sol pnt (i)
          for (CFuint jSolPnt = 0; jSolPnt < m_nbrSolSolDep; ++jSolPnt)
          {
            const CFuint jSolIdx = (*m_solSolDep)[m_pertSol][jSolPnt];
            
            m_tempFlux = 0.;

            // (d)
            for (CFuint iDim = 0; iDim < m_dim; ++iDim)
            {
              const CFreal polyCoef = (*m_solPolyDerivAtSolPnts)[jSolIdx][iDim][m_pertSol]; 
          
              m_tempFlux += (m_fluxJacobian[m_pertSide][m_pertSol][m_pertVar][iDim+m_ndimplus]) * polyCoef;
            }
            
            acc.addValues(jSolIdx+pertSideTerm,m_pertSol+pertSideTerm,m_pertVar,&m_tempFlux[0]);
          }
            //if (m_cells[LEFT]->getID() == 1 || m_cells[RIGHT]->getID() == 1) CFLog(INFO, "before: " << m_tempFlux << "\n");

          // (i)
          for (CFuint jSolIdx = 0; jSolIdx < m_nbrSolPnts; ++jSolIdx)
          {
            m_tempFlux = 0.0;
            
            // add the discontinuous gradient part of the jacobian (m) ////actually only need to loop over the depending which are also altered
            for (CFuint kSolPnt = 0; kSolPnt < m_nbrSolSolDep; ++kSolPnt)
            {
              const CFuint kSolIdx = (*m_solSolDep)[jSolIdx][kSolPnt];
              
              // (p)
              for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
              {
                const CFreal dgradVar_du = m_gradVarsToStateJacobian[m_pertSide][m_pertSol][m_pertVar][iEq];
                
                // (d)
                for (CFuint iDim = 0; iDim < m_dim; ++iDim)
                {
                  const CFreal dgradVar_du_dl = (*m_solPolyDerivAtSolPnts)[jSolIdx][iDim][kSolIdx] * dgradVar_du;
                
                  // (b)
                  for (CFuint jDim = 0; jDim < m_dim; ++jDim)
                  {
                    m_tempFlux += (m_gradientFluxJacobian[m_pertSide][kSolIdx][iEq][jDim][iDim+m_ndimplus]) * m_gradientStateJacobian[m_pertSide][kSolIdx][m_pertSide][m_pertSol][jDim] * dgradVar_du_dl;
//if (m_cells[LEFT]->getID() == 1 && jSolIdx+pertSideTerm==0 && m_pertSol+pertSideTerm==0 && m_pertVar==3) CFLog(INFO, "j: " << m_tempFlux << ", l: " << (*m_solPolyDerivAtSolPnts)[jSolIdx][iDim][kSolIdx] << ", dFdq: " << m_gradientFluxJacobian[m_pertSide][kSolIdx][iEq][jDim][iDim] 
//        << ", dqdu: " << m_gradientStateJacobian[m_pertSide][m_pertSol][m_pertSide][m_pertSol][0] << ", dudu: " << m_gradVarsToStateJacobian[m_pertSide][m_pertSol][m_pertVar][iEq] << "\n");
                  }
                }
              }
            }

            acc.addValues(jSolIdx+pertSideTerm,m_pertSol+pertSideTerm,m_pertVar,&m_tempFlux[0]);
//            if (m_cells[LEFT]->getID() == 1 && jSolIdx+pertSideTerm==0 && m_pertSol+pertSideTerm == 0 && m_pertVar==3) CFLog(INFO, "adding: " << m_tempFlux << "\n");
          }

          // add the discontinuous part of the jacobian related to the flx pnt (f)
          for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFlxDep; ++iFlxPnt)
          {
            const CFuint flxIdx = (*m_solFlxDep)[m_pertSol][iFlxPnt];
            
            // (df)
            const CFuint dim = (*m_flxPntFlxDim)[flxIdx];
            
            const RealVector dFduL = m_fluxJacobian[m_pertSide][m_pertSol][m_pertVar][dim] * (*m_solPolyValsAtFlxPnts)[flxIdx][m_pertSol];
            
            // add the second part of the discontinuous part of the jacobian (i)
            for (CFuint jSolPnt = 0; jSolPnt < m_nbrSolDep; ++jSolPnt)
            {
              const CFuint jSolIdx = (*m_flxSolDep)[flxIdx][jSolPnt];

              // get the divergence of the correction function
              const CFreal divh = m_corrFctDiv[jSolIdx][flxIdx];
                           
              m_tempFlux = -dFduL * divh;
            
              acc.addValues(jSolIdx+pertSideTerm,m_pertSol+pertSideTerm,m_pertVar,&m_tempFlux[0]);
            }
          }
              
          for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrTotalFlxPnts; ++iFlxPnt)
          {   
            // (df)
            const CFuint dim = (*m_flxPntFlxDim)[iFlxPnt];
              
            for (CFuint jSolPnt = 0; jSolPnt < m_nbrSolDep; ++jSolPnt)
            {
              const CFuint jSolIdx = (*m_flxSolDep)[iFlxPnt][jSolPnt];

              // get the divergence of the correction function
              const CFreal divh = m_corrFctDiv[jSolIdx][iFlxPnt];
              
              m_tempFlux = 0.0;
              
              // add the second part of the discontinuous gradient part of the jacobian (m) ////actually only need to loop over the depending which are also altered
              for (CFuint kSolPnt = 0; kSolPnt < m_nbrSolDep; ++kSolPnt)
              {
                const CFuint kSolIdx = (*m_flxSolDep)[iFlxPnt][kSolPnt];
                
                const CFreal divh_l = -divh * (*m_solPolyValsAtFlxPnts)[iFlxPnt][kSolIdx];
              
                // (b)
                for (CFuint jDim = 0; jDim < m_dim; ++jDim)
                {
                  const CFreal divh_l_dqdu = divh_l * m_gradientStateJacobian[m_pertSide][kSolIdx][m_pertSide][m_pertSol][jDim];
                  
                  // (p)
                  for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
                  {
                    m_tempFlux += divh_l_dqdu * (m_gradientFluxJacobian[m_pertSide][kSolIdx][iEq][jDim][dim]) * m_gradVarsToStateJacobian[m_pertSide][m_pertSol][m_pertVar][iEq];
                  }
                }
              }
              
              acc.addValues(jSolIdx+pertSideTerm,m_pertSol+pertSideTerm,m_pertVar,&m_tempFlux[0]);
            }
          }
        }
      }
    }
  }
  
  // loop over left and right cell to add the riemann flux (face) part to the jacobian
  for (m_pertSide = 0; m_pertSide < 2; ++m_pertSide)
  {
    // variable for the other side
    const CFuint iOtherSide = m_pertSide == LEFT ? RIGHT : LEFT;
    
    // cell ID of the cell at the non-perturbed side
    const CFuint otherCellID = m_cells[iOtherSide]->getID();
    
    // term depending on iSide
    const CFuint pertSideTerm = m_pertSide*m_nbrSolPnts;

    // term depending on iOtherSide
    const CFuint otherSideTerm = iOtherSide*m_nbrSolPnts;
    
    // loop over the variables in the state (k)
    for (m_pertVar = 0; m_pertVar < m_nbrEqs; ++m_pertVar)
    { 
      // loop over face flx pnts (f)
      for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
      {
        // local flux point indices in the left and right cell
        const CFuint flxPntIdxThis = (*m_faceFlxPntConnPerOrient)[m_orient][m_pertSide][iFlxPnt];
        const CFuint flxPntIdxOther = (*m_faceFlxPntConnPerOrient)[m_orient][iOtherSide][iFlxPnt];
        
        const RealVector dFIduJ = m_riemannFluxJacobian[m_pertSide][iFlxPnt][m_pertVar]*m_faceJacobVecSizeFlxPnts[iFlxPnt][m_pertSide];
        const RealVector dFIduJOther = m_riemannFluxJacobian[m_pertSide][iFlxPnt][m_pertVar]*m_faceJacobVecSizeFlxPnts[iFlxPnt][iOtherSide] ;
        
        const CFreal halfFaceJacob = 0.5 * m_faceJacobVecSizeFlxPnts[iFlxPnt][m_pertSide];
        
        // loop over the states to perturb the states (l)
        for (m_pertSol = 0; m_pertSol < m_nbrSolDep; ++m_pertSol)
        {
          const CFuint pertSolIdx = (*m_flxSolDep)[flxPntIdxThis][m_pertSol];
          
          const RealVector dFIduJL = dFIduJ * (*m_solPolyValsAtFlxPnts)[flxPntIdxThis][pertSolIdx];
          const RealVector dFIduJLOther = dFIduJOther * (*m_solPolyValsAtFlxPnts)[flxPntIdxThis][pertSolIdx];
            
          // add the second part of the discontinuous part of the jacobian (i)
          for (CFuint jSolPnt = 0; jSolPnt < m_nbrSolDep; ++jSolPnt)
          {
            const CFuint jSolIdxThis = (*m_flxSolDep)[flxPntIdxThis][jSolPnt];
            const CFuint jSolIdxOther = (*m_flxSolDep)[flxPntIdxOther][jSolPnt];

            // get the divergence of the correction function on this side
            CFreal divh = m_corrFctDiv[jSolIdxThis][flxPntIdxThis];
                          
            // add part on this side of face
            m_tempFlux = dFIduJL * divh;
              
            acc.addValues(jSolIdxThis+pertSideTerm,pertSolIdx+pertSideTerm,m_pertVar,&m_tempFlux[0]);
          
            // get the divergence of the correction function on other side
            divh = m_corrFctDiv[jSolIdxOther][flxPntIdxOther];
                        
            // add cross-cell part 
            m_tempFlux = dFIduJLOther * divh;   
              
            acc.addValues(jSolIdxOther+otherSideTerm,pertSolIdx+pertSideTerm,m_pertVar,&m_tempFlux[0]);
          }
        }
            
        // (i)
        if (m_addRiemannToGradJacob || m_addRiemannToGradCrossCellJacob)
        {
        for (CFuint jSolPnt = 0; jSolPnt < m_nbrSolDep; ++jSolPnt)
        {
          const CFuint jSolIdxThis = (*m_flxSolDep)[flxPntIdxThis][jSolPnt];
          const CFuint jSolIdxOther = (*m_flxSolDep)[flxPntIdxOther][jSolPnt];

          // get the divergence of the correction function on this side
          CFreal divh = m_corrFctDiv[jSolIdxThis][flxPntIdxThis];
            
          const CFreal divh_halfFaceJacob = divh * halfFaceJacob;
          
          // loop over the states to perturb the states (l)
          for (CFuint pertSolIdx = 0; pertSolIdx < m_nbrSolPnts; ++pertSolIdx)
          {
            if (m_addRiemannToGradJacob)
            {
            // add part on this side of face
            m_tempFlux = 0.0;

            // (m)
            for (CFuint kSolPnt = 0; kSolPnt < m_nbrSolDep; ++kSolPnt)
            {
              const CFuint kSolIdx = (*m_flxSolDep)[flxPntIdxThis][kSolPnt];
              const CFuint kSolIdxOther = (*m_flxSolDep)[flxPntIdxOther][kSolPnt];
              
              const CFreal divh_halfFaceJacob_l = divh_halfFaceJacob * (*m_solPolyValsAtFlxPnts)[flxPntIdxThis][kSolIdx];
              const CFreal divh_halfFaceJacob_lOther = divh_halfFaceJacob * (*m_solPolyValsAtFlxPnts)[flxPntIdxOther][kSolIdxOther];
                
              // (b)
              for (CFuint iDim = 0; iDim < m_dim; ++iDim)
              {
                const CFreal divh_halfFaceJacob_l_dqdu = divh_halfFaceJacob_l * m_gradientStateJacobian[m_pertSide][kSolIdx][m_pertSide][pertSolIdx][iDim];
                const CFreal divh_halfFaceJacob_l_dqduOther = divh_halfFaceJacob_lOther * m_gradientStateJacobian[iOtherSide][kSolIdxOther][m_pertSide][pertSolIdx][iDim];
                
                // (p)
                for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
                {
                  m_tempFlux += m_riemannFluxGradJacobian[iFlxPnt][iEq][iDim] * m_gradVarsToStateJacobian[m_pertSide][pertSolIdx][m_pertVar][iEq] * divh_halfFaceJacob_l_dqdu;
                    
                  m_tempFlux += m_riemannFluxGradJacobian[iFlxPnt][iEq][iDim] * m_gradVarsToStateJacobian[m_pertSide][pertSolIdx][m_pertVar][iEq] * divh_halfFaceJacob_l_dqduOther;
                }
              }         
            }
              
            acc.addValues(jSolIdxThis+pertSideTerm,pertSolIdx+pertSideTerm,m_pertVar,&m_tempFlux[0]);
            }
          
            if (m_addRiemannToGradCrossCellJacob)
            {
            // get the divergence of the correction function on other side
            divh = m_corrFctDiv[jSolIdxOther][flxPntIdxOther];
            
            const CFreal divh_halfFaceJacobOther = 0.5 * divh * m_faceJacobVecSizeFlxPnts[iFlxPnt][iOtherSide];

            // add cross-cell part 
            m_tempFlux = 0.0;

            // (m)
            for (CFuint kSolPnt = 0; kSolPnt < m_nbrSolDep; ++kSolPnt)
            {
              const CFuint kSolIdx = (*m_flxSolDep)[flxPntIdxThis][kSolPnt];
              const CFuint kSolIdxOther = (*m_flxSolDep)[flxPntIdxOther][kSolPnt];
              
              const CFreal divh_halfFaceJacobOther_lThis = divh_halfFaceJacobOther * (*m_solPolyValsAtFlxPnts)[flxPntIdxThis][kSolIdx];
              const CFreal divh_halfFaceJacobOther_lOther = divh_halfFaceJacobOther * (*m_solPolyValsAtFlxPnts)[flxPntIdxOther][kSolIdxOther];
                
              // (b)
              for (CFuint iDim = 0; iDim < m_dim; ++iDim)
              {
                const CFreal divh_halfFaceJacobOther_lThis_dqduThis = divh_halfFaceJacobOther_lThis * m_gradientStateJacobian[m_pertSide][kSolIdx][m_pertSide][pertSolIdx][iDim];
                const CFreal divh_halfFaceJacobOther_lOther_dqduOther = divh_halfFaceJacobOther_lOther * m_gradientStateJacobian[iOtherSide][kSolIdxOther][m_pertSide][pertSolIdx][iDim];
              
                // (p)
                for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
                {
                  m_tempFlux += m_riemannFluxGradJacobian[iFlxPnt][iEq][iDim] * m_gradVarsToStateJacobian[m_pertSide][pertSolIdx][m_pertVar][iEq] * divh_halfFaceJacobOther_lThis_dqduThis;
                    
                  m_tempFlux += m_riemannFluxGradJacobian[iFlxPnt][iEq][iDim] * m_gradVarsToStateJacobian[m_pertSide][pertSolIdx][m_pertVar][iEq] * divh_halfFaceJacobOther_lOther_dqduOther;
                }
              }         
            }
            
            acc.addValues(jSolIdxOther+otherSideTerm,pertSolIdx+pertSideTerm,m_pertVar,&m_tempFlux[0]);
            }
          }
        }
        }
        
        if (m_addFluxToGradCrossCellJacob)
        {
        //// add the cross-element gradient part
        for (CFuint jSolPnt = 0; jSolPnt < m_nbrSolPnts; ++jSolPnt)
        { 
          // loop over the states to perturb the states (l)
          for (m_pertSol = 0; m_pertSol < m_nbrSolDep; ++m_pertSol)
          {
            const CFuint pertSolIdx = (*m_flxSolDep)[flxPntIdxThis][m_pertSol];
              
            m_tempFlux = 0.0;
              
            // add the first and second part of the discontinuous gradient part of the jacobian
            for (CFuint iInfluencedFlx = 0; iInfluencedFlx < m_nbrFlxDep; ++iInfluencedFlx)
            {
              const CFuint iInfluencedFlxIdx = (*m_solFlxDep)[jSolPnt][iInfluencedFlx];
            
              // get the divergence of the correction function on this side
              CFreal divh = m_corrFctDiv[jSolPnt][iInfluencedFlxIdx]; 
              
              const CFuint dimOther = (*m_flxPntFlxDim)[iInfluencedFlxIdx];
                
              for (CFuint kSolPnt = 0; kSolPnt < m_nbrSolDep; ++kSolPnt)
              {
                const CFuint kSolIdxOther = (*m_flxSolDep)[iInfluencedFlxIdx][kSolPnt];

                const CFreal divh_lOther = -divh * (*m_solPolyValsAtFlxPnts)[iInfluencedFlxIdx][kSolIdxOther];
              
                for (CFuint jDim = 0; jDim < m_dim; ++jDim)
                {
                  const CFreal divh_l_dqduOther = divh_lOther * m_gradientStateJacobian[iOtherSide][kSolIdxOther][m_pertSide][pertSolIdx][jDim];
                  
                  for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
                  {
                    m_tempFlux += (m_gradientFluxJacobian[iOtherSide][kSolIdxOther][iEq][jDim][dimOther]) * m_gradVarsToStateJacobian[m_pertSide][pertSolIdx][m_pertVar][iEq] * divh_l_dqduOther;
                  }
                }
              }
            }
             
            for (CFuint kSolPnt = 0; kSolPnt < m_nbrSolSolDep; ++kSolPnt)
            {
              const CFuint kSolIdxOther = (*m_solSolDep)[jSolPnt][kSolPnt];
                
              for (CFuint iDim = 0; iDim < m_dim; ++iDim)
              {
                const CFreal lOther = (*m_solPolyDerivAtSolPnts)[jSolPnt][iDim][kSolIdxOther];
                  
                for (CFuint jDim = 0; jDim < m_dim; ++jDim)
                {
                  const CFreal l_dqduOther = lOther * m_gradientStateJacobian[iOtherSide][kSolIdxOther][m_pertSide][pertSolIdx][jDim];
                  
                  for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
                  {
                    m_tempFlux += (m_gradientFluxJacobian[iOtherSide][kSolIdxOther][iEq][jDim][iDim+m_ndimplus]) * m_gradVarsToStateJacobian[m_pertSide][pertSolIdx][m_pertVar][iEq] * l_dqduOther;
                  }
                }
              }
            }
            
            acc.addValues(jSolPnt+otherSideTerm,pertSolIdx+pertSideTerm,m_pertVar,&m_tempFlux[0]);
          }
        }
        }
      }
    }
  }
  
  if (m_cells[LEFT]->getID() == 1) 
  {
      //CFLog(INFO, "ACC: " << acc.getValue(0,4,3,3) << "\n");
      //acc.printToScreen();
  }
  if (m_cells[RIGHT]->getID() == 1) 
  {
      //CFLog(INFO, "ACC: " << acc.getValue(4,4,3,3) << "\n");
      //acc.printToScreen();
  }

  if (getMethodData().doComputeJacobian())
  {
    // add the values to the jacobian matrix
    m_lss->getMatrix()->addValues(acc);
  }

  // reset to zero the entries in the block accumulator
  acc.reset();
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffJacobFluxReconstruction::computeOneJacobDiffFaceTerm(const CFuint side)
{
  CFLog(VERBOSE, "computeOneJacobDiffFaceTerm\n");
    
  // get residual factor
  const CFreal resFactor = getMethodData().getResFactor();

  // dereference accumulator
  BlockAccumulator& acc = *m_acc;
  
  CFuint solIdx = 0;
  for (m_pertSide = 0; m_pertSide < 2; ++m_pertSide)
  {
    for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol, ++solIdx)
    {
      acc.setRowColIndex(solIdx,(*m_states[m_pertSide])[iSol]->getLocalID());
    }
  }

  //// compute the needed flux jacobians
  
  initJacobianComputation();
  
  computeCellFluxJacobianNum(resFactor);
  
  computeRiemannFluxJacobianNum(resFactor);
  
  computeFluxToGradJacobianNum(resFactor);
  
  if (m_addRiemannToGradJacob || m_addRiemannToGradCrossCellJacob) computeRiemannFluxToGradJacobianNum(resFactor);
  
  computeGradToStateJacobianAna();
  
  computeGradVarsToStateJacobianNum();
  
  //// add the total jacobians to the system jacobian
  
  // loop over left and right cell to add the discontinuous (cell) part to the jacobian
  for (m_pertSide = 0; m_pertSide < 2; ++m_pertSide)
  {
    // make sure this is only done once per cell
    if (!m_cellFlags[m_cells[m_pertSide]->getID()]) 
    {
      // variable for the other side
      const CFuint iOtherSide = m_pertSide == LEFT ? RIGHT : LEFT;

      // cell ID of the cell at the non-perturbed side
      const CFuint otherCellID = m_cells[iOtherSide]->getID();
    
      // term depending on iSide
      const CFuint pertSideTerm = m_pertSide*m_nbrSolPnts;

      // term depending on iOtherSide
      const CFuint otherSideTerm = iOtherSide*m_nbrSolPnts;

      // loop over the states to which to derive (l)
      for (m_pertSol = 0; m_pertSol < m_nbrSolPnts; ++m_pertSol)
      {
        // loop over the variables in the state (k)
        for (m_pertVar = 0; m_pertVar < m_nbrEqs; ++m_pertVar)
        { 
          // add the discontinuous part of the jacobian related to the sol pnt (i)
          for (CFuint jSolPnt = 0; jSolPnt < m_nbrSolSolDep; ++jSolPnt)
          {
            const CFuint jSolIdx = (*m_solSolDep)[m_pertSol][jSolPnt];
            
            m_tempFlux = 0.;

            // (d)
            for (CFuint iDim = 0; iDim < m_dim; ++iDim)
            {
              const CFreal polyCoef = (*m_solPolyDerivAtSolPnts)[jSolIdx][iDim][m_pertSol]; 
          
              m_tempFlux += (m_fluxJacobian[m_pertSide][m_pertSol][m_pertVar][iDim+m_ndimplus]) * polyCoef;
            }
            
            acc.addValues(jSolIdx+pertSideTerm,m_pertSol+pertSideTerm,m_pertVar,&m_tempFlux[0]);
          }
            //if (m_cells[LEFT]->getID() == 1 || m_cells[RIGHT]->getID() == 1) CFLog(INFO, "before: " << m_tempFlux << "\n");
          
          // (i)
          for (CFuint jSolIdx = 0; jSolIdx < m_nbrSolPnts; ++jSolIdx)
          {
            m_tempFlux = 0.0;
            
            // add the discontinuous gradient part of the jacobian (m) ////actually only need to loop over the depending which are also altered
            for (CFuint kSolPnt = 0; kSolPnt < m_nbrSolSolDep; ++kSolPnt)
            {
              const CFuint kSolIdx = (*m_solSolDep)[jSolIdx][kSolPnt];
              
              // (p)
              for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
              {
                const CFreal dgradVar_du = m_gradVarsToStateJacobian[m_pertSide][m_pertSol][m_pertVar][iEq];
                
                // (d)
                for (CFuint iDim = 0; iDim < m_dim; ++iDim)
                {
                  const CFreal dgradVar_du_dl = (*m_solPolyDerivAtSolPnts)[jSolIdx][iDim][kSolIdx] * dgradVar_du;
                
                  // (b)
                  for (CFuint jDim = 0; jDim < m_dim; ++jDim)
                  {
                    m_tempFlux += (m_gradientFluxJacobian[m_pertSide][kSolIdx][iEq][jDim][iDim+m_ndimplus]) * m_gradientStateJacobian[m_pertSide][kSolIdx][m_pertSide][m_pertSol][jDim] * dgradVar_du_dl;
//if (m_cells[LEFT]->getID() == 1 && jSolIdx+pertSideTerm==0 && m_pertSol+pertSideTerm==0 && m_pertVar==3) CFLog(INFO, "j: " << m_tempFlux << ", l: " << (*m_solPolyDerivAtSolPnts)[jSolIdx][iDim][kSolIdx] << ", dFdq: " << m_gradientFluxJacobian[m_pertSide][kSolIdx][iEq][jDim][iDim] 
//        << ", dqdu: " << m_gradientStateJacobian[m_pertSide][m_pertSol][m_pertSide][m_pertSol][0] << ", dudu: " << m_gradVarsToStateJacobian[m_pertSide][m_pertSol][m_pertVar][iEq] << "\n");
                  }
                }
              }
            }
            
            acc.addValues(jSolIdx+pertSideTerm,m_pertSol+pertSideTerm,m_pertVar,&m_tempFlux[0]);
//            if (m_cells[LEFT]->getID() == 1 && jSolIdx+pertSideTerm==0 && m_pertSol+pertSideTerm == 0 && m_pertVar==3) CFLog(INFO, "adding: " << m_tempFlux << "\n");
          }
          
          // add the discontinuous part of the jacobian related to the flx pnt (f)
          for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFlxDep; ++iFlxPnt)
          {
            const CFuint flxIdx = (*m_solFlxDep)[m_pertSol][iFlxPnt];
            
            // (df)
            const CFuint dim = (*m_flxPntFlxDim)[flxIdx];
            
            const RealVector dFduL = m_fluxJacobian[m_pertSide][m_pertSol][m_pertVar][dim] * (*m_solPolyValsAtFlxPnts)[flxIdx][m_pertSol];
            
            // add the second part of the discontinuous part of the jacobian (i)
            for (CFuint jSolPnt = 0; jSolPnt < m_nbrSolDep; ++jSolPnt)
            {
              const CFuint jSolIdx = (*m_flxSolDep)[flxIdx][jSolPnt];

              // get the divergence of the correction function
              const CFreal divh = m_corrFctDiv[jSolIdx][flxIdx];
                           
              m_tempFlux = -dFduL * divh;
            
              acc.addValues(jSolIdx+pertSideTerm,m_pertSol+pertSideTerm,m_pertVar,&m_tempFlux[0]);
            }
          }
              
          for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrTotalFlxPnts; ++iFlxPnt)
          {   
            // (df)
            const CFuint dim = (*m_flxPntFlxDim)[iFlxPnt];
              
            for (CFuint jSolPnt = 0; jSolPnt < m_nbrSolDep; ++jSolPnt)
            {
              const CFuint jSolIdx = (*m_flxSolDep)[iFlxPnt][jSolPnt];

              // get the divergence of the correction function
              const CFreal divh = m_corrFctDiv[jSolIdx][iFlxPnt];
              
              m_tempFlux = 0.0;
              
              // add the second part of the discontinuous gradient part of the jacobian (m) ////actually only need to loop over the depending which are also altered
              for (CFuint kSolPnt = 0; kSolPnt < m_nbrSolDep; ++kSolPnt)
              {
                const CFuint kSolIdx = (*m_flxSolDep)[iFlxPnt][kSolPnt];
                
                const CFreal divh_l = -divh * (*m_solPolyValsAtFlxPnts)[iFlxPnt][kSolIdx];
              
                // (b)
                for (CFuint jDim = 0; jDim < m_dim; ++jDim)
                {
                  const CFreal divh_l_dqdu = divh_l * m_gradientStateJacobian[m_pertSide][kSolIdx][m_pertSide][m_pertSol][jDim];
                  
                  // (p)
                  for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
                  {
                    m_tempFlux += divh_l_dqdu * (m_gradientFluxJacobian[m_pertSide][kSolIdx][iEq][jDim][dim]) * m_gradVarsToStateJacobian[m_pertSide][m_pertSol][m_pertVar][iEq];
                  }
                }
              }
              
              acc.addValues(jSolIdx+pertSideTerm,m_pertSol+pertSideTerm,m_pertVar,&m_tempFlux[0]);
            }
          }
        }
      }
    }
  }
  
  // loop over left and right cell to add the riemann flux (face) part to the jacobian
  for (m_pertSide = 0; m_pertSide < 2; ++m_pertSide)
  {
    // variable for the other side
    const CFuint iOtherSide = m_pertSide == LEFT ? RIGHT : LEFT;
    
    // cell ID of the cell at the non-perturbed side
    const CFuint otherCellID = m_cells[iOtherSide]->getID();
    
    // term depending on iSide
    const CFuint pertSideTerm = m_pertSide*m_nbrSolPnts;

    // term depending on iOtherSide
    const CFuint otherSideTerm = iOtherSide*m_nbrSolPnts;
    
    // loop over the variables in the state (k)
    for (m_pertVar = 0; m_pertVar < m_nbrEqs; ++m_pertVar)
    { 
      // loop over face flx pnts (f)
      for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
      {
        // local flux point indices in the left and right cell
        const CFuint flxPntIdxThis = (*m_faceFlxPntConnPerOrient)[m_orient][m_pertSide][iFlxPnt];
        const CFuint flxPntIdxOther = (*m_faceFlxPntConnPerOrient)[m_orient][iOtherSide][iFlxPnt];
        
        const RealVector dFIduJ = m_riemannFluxJacobian[m_pertSide][iFlxPnt][m_pertVar]*m_faceJacobVecSizeFlxPnts[iFlxPnt][m_pertSide];
        const RealVector dFIduJOther = m_riemannFluxJacobian[m_pertSide][iFlxPnt][m_pertVar]*m_faceJacobVecSizeFlxPnts[iFlxPnt][iOtherSide] ;
        
        const CFreal halfFaceJacob = 0.5 * m_faceJacobVecSizeFlxPnts[iFlxPnt][m_pertSide];
        
        // loop over the states to perturb the states (l)
        for (m_pertSol = 0; m_pertSol < m_nbrSolDep; ++m_pertSol)
        {
          const CFuint pertSolIdx = (*m_flxSolDep)[flxPntIdxThis][m_pertSol];
          
          const RealVector dFIduJL = dFIduJ * (*m_solPolyValsAtFlxPnts)[flxPntIdxThis][pertSolIdx];
          const RealVector dFIduJLOther = dFIduJOther * (*m_solPolyValsAtFlxPnts)[flxPntIdxThis][pertSolIdx];
            
          // add the second part of the discontinuous part of the jacobian (i)
          for (CFuint jSolPnt = 0; jSolPnt < m_nbrSolDep; ++jSolPnt)
          {
            const CFuint jSolIdxThis = (*m_flxSolDep)[flxPntIdxThis][jSolPnt];
            const CFuint jSolIdxOther = (*m_flxSolDep)[flxPntIdxOther][jSolPnt];

            // get the divergence of the correction function on this side
            CFreal divh = m_corrFctDiv[jSolIdxThis][flxPntIdxThis];
                          
            // add part on this side of face
            m_tempFlux = dFIduJL * divh;
              
            acc.addValues(jSolIdxThis+pertSideTerm,pertSolIdx+pertSideTerm,m_pertVar,&m_tempFlux[0]);
          
            // get the divergence of the correction function on other side
            divh = m_corrFctDiv[jSolIdxOther][flxPntIdxOther];
                        
            // add cross-cell part 
            m_tempFlux = dFIduJLOther * divh;   
              
            acc.addValues(jSolIdxOther+otherSideTerm,pertSolIdx+pertSideTerm,m_pertVar,&m_tempFlux[0]);
          }
        }
            
        // (i)
        if (m_addRiemannToGradJacob || m_addRiemannToGradCrossCellJacob)
        {
        for (CFuint jSolPnt = 0; jSolPnt < m_nbrSolDep; ++jSolPnt)
        {
          const CFuint jSolIdxThis = (*m_flxSolDep)[flxPntIdxThis][jSolPnt];
          const CFuint jSolIdxOther = (*m_flxSolDep)[flxPntIdxOther][jSolPnt];

          // get the divergence of the correction function on this side
          CFreal divh = m_corrFctDiv[jSolIdxThis][flxPntIdxThis];
            
          const CFreal divh_halfFaceJacob = divh * halfFaceJacob;
          
          // loop over the states to perturb the states (l)
          for (CFuint pertSolIdx = 0; pertSolIdx < m_nbrSolPnts; ++pertSolIdx)
          {
            if (m_addRiemannToGradJacob)
            {
            // add part on this side of face
            m_tempFlux = 0.0;

            // (m)
            for (CFuint kSolPnt = 0; kSolPnt < m_nbrSolDep; ++kSolPnt)
            {
              const CFuint kSolIdx = (*m_flxSolDep)[flxPntIdxThis][kSolPnt];
              const CFuint kSolIdxOther = (*m_flxSolDep)[flxPntIdxOther][kSolPnt];
              
              const CFreal divh_halfFaceJacob_l = divh_halfFaceJacob * (*m_solPolyValsAtFlxPnts)[flxPntIdxThis][kSolIdx];
              const CFreal divh_halfFaceJacob_lOther = divh_halfFaceJacob * (*m_solPolyValsAtFlxPnts)[flxPntIdxOther][kSolIdxOther];
                
              // (b)
              for (CFuint iDim = 0; iDim < m_dim; ++iDim)
              {
                const CFreal divh_halfFaceJacob_l_dqdu = divh_halfFaceJacob_l * m_gradientStateJacobian[m_pertSide][kSolIdx][m_pertSide][pertSolIdx][iDim];
                const CFreal divh_halfFaceJacob_l_dqduOther = divh_halfFaceJacob_lOther * m_gradientStateJacobian[iOtherSide][kSolIdxOther][m_pertSide][pertSolIdx][iDim];
                
                // (p)
                for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
                {
                  m_tempFlux += m_riemannFluxGradJacobian[iFlxPnt][iEq][iDim] * m_gradVarsToStateJacobian[m_pertSide][pertSolIdx][m_pertVar][iEq] * divh_halfFaceJacob_l_dqdu;
                    
                  m_tempFlux += m_riemannFluxGradJacobian[iFlxPnt][iEq][iDim] * m_gradVarsToStateJacobian[m_pertSide][pertSolIdx][m_pertVar][iEq] * divh_halfFaceJacob_l_dqduOther;
                }
              }         
            }
              
            acc.addValues(jSolIdxThis+pertSideTerm,pertSolIdx+pertSideTerm,m_pertVar,&m_tempFlux[0]);
            }
          
            if (m_addRiemannToGradCrossCellJacob)
            {
            // get the divergence of the correction function on other side
            divh = m_corrFctDiv[jSolIdxOther][flxPntIdxOther];
            
            const CFreal divh_halfFaceJacobOther = 0.5 * divh * m_faceJacobVecSizeFlxPnts[iFlxPnt][iOtherSide];

            // add cross-cell part 
            m_tempFlux = 0.0;

            // (m)
            for (CFuint kSolPnt = 0; kSolPnt < m_nbrSolDep; ++kSolPnt)
            {
              const CFuint kSolIdx = (*m_flxSolDep)[flxPntIdxThis][kSolPnt];
              const CFuint kSolIdxOther = (*m_flxSolDep)[flxPntIdxOther][kSolPnt];
              
              const CFreal divh_halfFaceJacobOther_lThis = divh_halfFaceJacobOther * (*m_solPolyValsAtFlxPnts)[flxPntIdxThis][kSolIdx];
              const CFreal divh_halfFaceJacobOther_lOther = divh_halfFaceJacobOther * (*m_solPolyValsAtFlxPnts)[flxPntIdxOther][kSolIdxOther];
                
              // (b)
              for (CFuint iDim = 0; iDim < m_dim; ++iDim)
              {
                const CFreal divh_halfFaceJacobOther_lThis_dqduThis = divh_halfFaceJacobOther_lThis * m_gradientStateJacobian[m_pertSide][kSolIdx][m_pertSide][pertSolIdx][iDim];
                const CFreal divh_halfFaceJacobOther_lOther_dqduOther = divh_halfFaceJacobOther_lOther * m_gradientStateJacobian[iOtherSide][kSolIdxOther][m_pertSide][pertSolIdx][iDim];
              
                // (p)
                for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
                {
                  m_tempFlux += m_riemannFluxGradJacobian[iFlxPnt][iEq][iDim] * m_gradVarsToStateJacobian[m_pertSide][pertSolIdx][m_pertVar][iEq] * divh_halfFaceJacobOther_lThis_dqduThis;
                    
                  m_tempFlux += m_riemannFluxGradJacobian[iFlxPnt][iEq][iDim] * m_gradVarsToStateJacobian[m_pertSide][pertSolIdx][m_pertVar][iEq] * divh_halfFaceJacobOther_lOther_dqduOther;
                }
              }         
            }
            
            acc.addValues(jSolIdxOther+otherSideTerm,pertSolIdx+pertSideTerm,m_pertVar,&m_tempFlux[0]);
            }
          }
        }
        }
        
        if (m_addFluxToGradCrossCellJacob)
        {
        //// add the cross-element gradient part
        for (CFuint jSolPnt = 0; jSolPnt < m_nbrSolPnts; ++jSolPnt)
        { 
          // loop over the states to perturb the states (l)
          for (m_pertSol = 0; m_pertSol < m_nbrSolDep; ++m_pertSol)
          {
            const CFuint pertSolIdx = (*m_flxSolDep)[flxPntIdxThis][m_pertSol];
              
            m_tempFlux = 0.0;
              
            // add the first and second part of the discontinuous gradient part of the jacobian
            for (CFuint iInfluencedFlx = 0; iInfluencedFlx < m_nbrFlxDep; ++iInfluencedFlx)
            {
              const CFuint iInfluencedFlxIdx = (*m_solFlxDep)[jSolPnt][iInfluencedFlx];
            
              // get the divergence of the correction function on this side
              CFreal divh = m_corrFctDiv[jSolPnt][iInfluencedFlxIdx]; 
              
              const CFuint dimOther = (*m_flxPntFlxDim)[iInfluencedFlxIdx];
                
              for (CFuint kSolPnt = 0; kSolPnt < m_nbrSolDep; ++kSolPnt)
              {
                const CFuint kSolIdxOther = (*m_flxSolDep)[iInfluencedFlxIdx][kSolPnt];

                const CFreal divh_lOther = -divh * (*m_solPolyValsAtFlxPnts)[iInfluencedFlxIdx][kSolIdxOther];
              
                for (CFuint jDim = 0; jDim < m_dim; ++jDim)
                {
                  const CFreal divh_l_dqduOther = divh_lOther * m_gradientStateJacobian[iOtherSide][kSolIdxOther][m_pertSide][pertSolIdx][jDim];
                  
                  for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
                  {
                    m_tempFlux += (m_gradientFluxJacobian[iOtherSide][kSolIdxOther][iEq][jDim][dimOther]) * m_gradVarsToStateJacobian[m_pertSide][pertSolIdx][m_pertVar][iEq] * divh_l_dqduOther;
                  }
                }
              }
            }
             
            for (CFuint kSolPnt = 0; kSolPnt < m_nbrSolSolDep; ++kSolPnt)
            {
              const CFuint kSolIdxOther = (*m_solSolDep)[jSolPnt][kSolPnt];
                
              for (CFuint iDim = 0; iDim < m_dim; ++iDim)
              {
                const CFreal lOther = (*m_solPolyDerivAtSolPnts)[jSolPnt][iDim][kSolIdxOther];
                  
                for (CFuint jDim = 0; jDim < m_dim; ++jDim)
                {
                  const CFreal l_dqduOther = lOther * m_gradientStateJacobian[iOtherSide][kSolIdxOther][m_pertSide][pertSolIdx][jDim];
                  
                  for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
                  {
                    m_tempFlux += (m_gradientFluxJacobian[iOtherSide][kSolIdxOther][iEq][jDim][iDim+m_ndimplus]) * m_gradVarsToStateJacobian[m_pertSide][pertSolIdx][m_pertVar][iEq] * l_dqduOther;
                  }
                }
              }
            }
            
            acc.addValues(jSolPnt+otherSideTerm,pertSolIdx+pertSideTerm,m_pertVar,&m_tempFlux[0]);
          }
        }
        }
      }
    }
  }
  
  if (m_cells[LEFT]->getID() == 1) 
  {
      //CFLog(INFO, "ACC: " << acc.getValue(0,4,3,3) << "\n");
      //acc.printToScreen();
  }
  if (m_cells[RIGHT]->getID() == 1) 
  {
      //CFLog(INFO, "ACC: " << acc.getValue(4,4,3,3) << "\n");
      //acc.printToScreen();
  }

  if (getMethodData().doComputeJacobian())
  {
    // add the values to the jacobian matrix
    m_lss->getMatrix()->addValues(acc);
  }

  // reset to zero the entries in the block accumulator
  acc.reset();
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffJacobFluxReconstruction::computeUnpertCellDiffResiduals(const CFuint side)
{
  // get the sol pnt normals
  DataHandle< CFreal > solPntNormals = socket_solPntNormals.getDataHandle();
  
  // store the sol pnt normals
  for (CFuint iState = 0; iState < m_nbrSolPnts; ++iState)
  {
    const CFuint solID = (*(m_states[side]))[iState]->getLocalID();
      
    for (CFuint iDim = 0; iDim < m_dim+m_ndimplus; ++iDim)
    {
      for (CFuint jDim = 0; jDim < m_dim; ++jDim)
      {
        m_cellFluxProjVects[iDim][iState][jDim] = solPntNormals[solID*(m_dim+m_ndimplus)*m_dim+iDim*m_dim+jDim];
      }
    }
  }
  
  // reset the extrapolated fluxes
  for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrTotalFlxPnts; ++iFlxPnt)
  {
    m_extrapolatedFluxes[iFlxPnt] = 0.0;
  }

  // Loop over solution points to calculate the discontinuous flux.
  for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolPnts; ++iSolPnt)
  { 
    for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
    {
      *(m_tempGrad[iVar]) = (*(m_cellGrads[side][iSolPnt]))[iVar];
    }
    
    m_updateVarSet->computePhysicalData(*((*(m_states[side]))[iSolPnt]), m_pData); 

    m_avgSol = *((*(m_states[side]))[iSolPnt]->getData());

    prepareFluxComputation();

    // calculate the discontinuous flux projected on x, y, z-directions
    for (CFuint iDim = 0; iDim < m_dim+m_ndimplus; ++iDim)
    {
      // add diffusive part 
      computeFlux(m_avgSol,m_tempGrad,m_cellFluxProjVects[iDim][iSolPnt],0,m_contFlx[iSolPnt][iDim]);
      
      // add convective part
      m_contFlx[iSolPnt][iDim] -= m_updateVarSet->getFlux()(m_pData,m_cellFluxProjVects[iDim][iSolPnt]);
    }

    for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFlxDep; ++iFlxPnt)
    {
      const CFuint flxIdx = (*m_solFlxDep)[iSolPnt][iFlxPnt];
      const CFuint dim = (*m_flxPntFlxDim)[flxIdx];

      m_extrapolatedFluxes[flxIdx] += (*m_solPolyValsAtFlxPnts)[flxIdx][iSolPnt]*(m_contFlx[iSolPnt][dim]);
    }
  }

  // Loop over solution pnts to calculate the divergence of the discontinuous flux
  for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolPnts; ++iSolPnt)
  {
    // reset the divergence of FC
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      m_unpertCellDiffRes[side][m_nbrEqs*iSolPnt+iEq] = 0.0;
    }

    // Loop over solution pnt to count factor of all sol pnt polys
    for (CFuint jSolPnt = 0; jSolPnt < m_nbrSolSolDep; ++jSolPnt)
    {
      const CFuint jSolIdx = (*m_solSolDep)[iSolPnt][jSolPnt];

      // Loop over deriv directions and sum them to compute divergence
      for (CFuint iDir = 0; iDir < m_dim; ++iDir)
      {
        const CFreal polyCoef = (*m_solPolyDerivAtSolPnts)[iSolPnt][iDir][jSolIdx]; 

        // Loop over conservative fluxes 
        for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
        {
          // Store divFD in the vector that will be divFC
          m_unpertCellDiffRes[side][m_nbrEqs*iSolPnt+iEq] += polyCoef*(m_contFlx[jSolIdx][iDir+m_ndimplus][iEq]);
	}
      }
    }

    for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFlxDep; ++iFlxPnt)
    {
      const CFuint flxIdx = (*m_solFlxDep)[iSolPnt][iFlxPnt];

      // get the divergence of the correction function
      const CFreal divh = m_corrFctDiv[iSolPnt][flxIdx];
  
      // Fill in the corrections
      for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
      {
        m_unpertCellDiffRes[side][m_nbrEqs*iSolPnt+iVar] += -m_extrapolatedFluxes[flxIdx][iVar] * divh; 
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffJacobFluxReconstruction::computeGradientFaceCorrections()
{ 
  // get the face flux point normals
  DataHandle< CFreal > flxPntNormals = socket_flxPntNormals.getDataHandle();
  
  // Loop over flux points to set the normal vectors
  for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
  { 
    for (CFuint iDim = 0; iDim < m_dim; ++iDim)
    {
      m_faceJacobVecs[iFlxPnt][iDim] = flxPntNormals[m_face->getID()*m_nbrFaceFlxPnts*m_dim+iFlxPnt*m_dim+iDim];
    }
  }

  // Loop over solution pnts to reset the grad updates
  for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolPnts; ++iSolPnt)
  {
    // Loop over  variables
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      //set the grad updates to 0 
      m_gradUpdates[LEFT][iSolPnt][iEq] = 0.0;
      m_gradUpdates[RIGHT][iSolPnt][iEq] = 0.0;
    }
  }
  
  // compute the face corrections to the gradients
  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    const CFuint flxIdxL = (*m_faceFlxPntConnPerOrient)[m_orient][LEFT][iFlx];
    const CFuint flxIdxR = (*m_faceFlxPntConnPerOrient)[m_orient][RIGHT][iFlx];
    
    // reset states in flx pnt
    *(m_cellStatesFlxPnt[LEFT][iFlx]) = 0.0;
    *(m_cellStatesFlxPnt[RIGHT][iFlx]) = 0.0;

    // extrapolate the left and right states to the flx pnts
    for (CFuint iSol = 0; iSol < m_nbrSolDep; ++iSol)
    {
      const CFuint solIdxL = (*m_flxSolDep)[flxIdxL][iSol];
      const CFuint solIdxR = (*m_flxSolDep)[flxIdxR][iSol];
 
      // add the contributions of the current sol pnt
      *(m_cellStatesFlxPnt[LEFT][iFlx]) += (*m_solPolyValsAtFlxPnts)[flxIdxL][solIdxL]*(*((*(m_states[LEFT]))[solIdxL]));
      *(m_cellStatesFlxPnt[RIGHT][iFlx]) += (*m_solPolyValsAtFlxPnts)[flxIdxR][solIdxR]*(*((*(m_states[RIGHT]))[solIdxR]));
    }

    // Loop over  variables
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      const CFreal avgSol = ((*m_cellStatesFlxPnt[LEFT][iFlx])[iEq]+(*m_cellStatesFlxPnt[RIGHT][iFlx])[iEq])/2.0;
      m_projectedCorrL = (avgSol-(*m_cellStatesFlxPnt[LEFT][iFlx])[iEq])*(*m_faceMappedCoordDir)[m_orient][LEFT]*m_faceJacobVecs[iFlx];
      m_projectedCorrR = (avgSol-(*m_cellStatesFlxPnt[RIGHT][iFlx])[iEq])*(*m_faceMappedCoordDir)[m_orient][RIGHT]*m_faceJacobVecs[iFlx];

      // Loop over solution pnts to calculate the grad updates
      for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolDep; ++iSolPnt)
      {
        const CFuint iSolIdxL = (*m_flxSolDep)[flxIdxL][iSolPnt];
        const CFuint iSolIdxR = (*m_flxSolDep)[flxIdxR][iSolPnt];

	/// @todo Check if this is also OK for triangles!!
	m_gradUpdates[LEFT][iSolIdxL][iEq] += m_projectedCorrL*m_corrFctDiv[iSolIdxL][flxIdxL];
	m_gradUpdates[RIGHT][iSolIdxR][iEq] += m_projectedCorrR*m_corrFctDiv[iSolIdxR][flxIdxR];
      }
    }
  }
  
  // get the gradients
  DataHandle< vector< RealVector > > gradients = socket_gradients.getDataHandle();

  for (CFuint iSide = 0; iSide < 2; ++iSide)
  {
    for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
    {
      // get state ID
      const CFuint solID = (*m_states[iSide])[iSol]->getLocalID();

      // update gradients
      for (CFuint iGrad = 0; iGrad < m_nbrEqs; ++iGrad)
      {
        gradients[solID][iGrad] += m_gradUpdates[iSide][iSol][iGrad];
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffJacobFluxReconstruction::computeGradients()
{
  // get the sol pnt normals
  DataHandle< CFreal > solPntNormals = socket_solPntNormals.getDataHandle();
  
  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      //set the grad updates to 0 
      m_gradUpdates[0][iSol][iEq] = 0.0;
    }
  }
  
  // Loop over solution pnts to calculate the grad updates
  for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolPnts; ++iSolPnt)
  {
    const CFuint solID = (*m_cellStates)[iSolPnt]->getLocalID();
    
    // Loop over  variables
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      // Loop over gradient directions
      for (CFuint iDir = 0; iDir < m_dim; ++iDir)
      {
        for (CFuint jDir = 0; jDir < m_dim; ++jDir)
        {
	  // project the state on a normal and reuse a RealVector variable of the class to store
	  m_projectedCorrL[jDir] = ((*(*m_cellStates)[iSolPnt])[iEq]) * solPntNormals[solID*(m_dim+m_ndimplus)*m_dim+(iDir+m_ndimplus)*m_dim+jDir];
        }

        // Loop over solution pnts to count factor of all sol pnt polys
        for (CFuint jSolPnt = 0; jSolPnt < m_nbrSolSolDep; ++jSolPnt)
        { 
          const CFuint jSolIdx = (*m_solSolDep)[iSolPnt][jSolPnt];
          // compute the grad updates
          m_gradUpdates[0][jSolIdx][iEq] += (*m_solPolyDerivAtSolPnts)[jSolIdx][iDir][iSolPnt]*m_projectedCorrL;
	}
      }
    }
  }
  
  // get the gradients
  DataHandle< vector< RealVector > > gradients = socket_gradients.getDataHandle();
  
  // get the volumes
  DataHandle<CFreal> volumes = socket_volumes.getDataHandle();

//  // get jacobian determinants at solution points
//  m_jacobDet = m_cell->computeGeometricShapeFunctionJacobianDeterminant(*m_solPntsLocalCoords);

  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {
    // get state ID
    const CFuint solID = (*m_cellStates)[iSol]->getLocalID();

    // inverse Jacobian determinant
    const CFreal invJacobDet = 1.0/volumes[solID];
    
    // update gradients
    for (CFuint iGrad = 0; iGrad < m_nbrEqs; ++iGrad)
    {
      gradients[solID][iGrad] += m_gradUpdates[0][iSol][iGrad];
      gradients[solID][iGrad] *= invJacobDet;
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffJacobFluxReconstruction::computeFlux(const RealVector& values, const std::vector< RealVector* >& gradients, const RealVector& normal, const CFreal& radius, RealVector& flux)
{
  flux = m_diffusiveVarSet->getFlux(values,gradients,normal,radius);
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffJacobFluxReconstruction::setup()
{
  CFAUTOTRACE;

  // setup parent class
  DiffRHSJacobFluxReconstruction::setup();
  
  m_unpertAllCellDiffRes.resize(0);
  
  // get the update varset
  m_updateVarSet = getMethodData().getUpdateVar();
  
  // resize the physical data temporary vector
  SafePtr<BaseTerm> convTerm = PhysicalModelStack::getActive()->getImplementor()->getConvectiveTerm(); 
  convTerm->resizePhysicalData(m_pData);
  
  // get the local spectral FD data
  vector< FluxReconstructionElementData* >& frLocalData = getMethodData().getFRLocalData();
  cf_assert(frLocalData.size() > 0);
  // for now, there should be only one type of element
  cf_assert(frLocalData.size() == 1);
  
  // get the elementTypeData
  SafePtr< vector<ElementTypeData> > elemType = MeshDataStack::getActive()->getElementTypeData();
  
  // get the number of cells in the mesh
  const CFuint nbrCells = (*elemType)[0].getEndIdx();
  
  const CFuint nbStates = nbrCells*m_nbrSolPnts;
  
  m_tempSolPntVec.resize(m_nbrSolPnts);
  m_tempSolPntVec2.resize(m_nbrSolPnts);
  m_fluxJacobian.resize(2);
  m_fluxJacobian[LEFT].resize(m_nbrSolPnts);
  m_fluxJacobian[RIGHT].resize(m_nbrSolPnts);
  m_riemannFluxJacobian.resize(2);
  m_riemannFluxJacobian[LEFT].resize(m_nbrFaceFlxPnts);
  m_riemannFluxJacobian[RIGHT].resize(m_nbrFaceFlxPnts);
  m_flxPntRiemannFluxDiff.resize(m_nbrFaceFlxPnts);
  m_flxPntRiemannFluxPert.resize(m_nbrFaceFlxPnts);
  m_tempFlux.resize(m_nbrEqs);
  m_gradientFluxJacobian.resize(2);
  m_gradientFluxJacobian[LEFT].resize(m_nbrSolPnts);
  m_gradientFluxJacobian[RIGHT].resize(m_nbrSolPnts);
  m_gradientStateJacobian.resize(2);
  m_gradientStateJacobian[LEFT].resize(m_nbrSolPnts);
  m_gradientStateJacobian[RIGHT].resize(m_nbrSolPnts);
  m_riemannFluxGradJacobian.resize(m_nbrFaceFlxPnts);
  m_gradVarsToStateJacobian.resize(2);
  m_gradVarsToStateJacobian[LEFT].resize(m_nbrSolPnts);
  m_gradVarsToStateJacobian[RIGHT].resize(m_nbrSolPnts);
  m_contFlxBackupDiff.resize(2);
  m_contFlxBackupDiff[LEFT].resize(m_nbrSolPnts);
  m_contFlxBackupDiff[RIGHT].resize(m_nbrSolPnts);
  
  for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolPnts; ++iSolPnt)
  {
    m_contFlxBackupDiff[LEFT][iSolPnt].resize(m_dim+m_ndimplus);
    m_contFlxBackupDiff[RIGHT][iSolPnt].resize(m_dim+m_ndimplus);
    
    for (CFuint iDim = 0; iDim < m_dim+m_ndimplus; ++iDim)
    {
      m_contFlxBackupDiff[LEFT][iSolPnt][iDim].resize(m_nbrEqs);
      m_contFlxBackupDiff[RIGHT][iSolPnt][iDim].resize(m_nbrEqs);
    }
  }
  
  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  { 
    m_fluxJacobian[LEFT][iSol].resize(m_nbrEqs);
    m_fluxJacobian[RIGHT][iSol].resize(m_nbrEqs);
    m_gradientFluxJacobian[LEFT][iSol].resize(m_nbrEqs);
    m_gradientFluxJacobian[RIGHT][iSol].resize(m_nbrEqs);
    m_gradientStateJacobian[LEFT][iSol].resize(2);
    m_gradientStateJacobian[RIGHT][iSol].resize(2);
    m_gradVarsToStateJacobian[LEFT][iSol].resize(m_nbrEqs);
    m_gradVarsToStateJacobian[RIGHT][iSol].resize(m_nbrEqs);
    
    for (CFuint iSide = 0; iSide < 2; ++iSide)
    {
      m_gradientStateJacobian[LEFT][iSol][iSide].resize(m_nbrSolPnts);
      m_gradientStateJacobian[RIGHT][iSol][iSide].resize(m_nbrSolPnts);
      
      for (CFuint jSol = 0; jSol < m_nbrSolPnts; ++jSol)
      {
        m_gradientStateJacobian[LEFT][iSol][iSide][jSol].resize(m_dim);
        m_gradientStateJacobian[RIGHT][iSol][iSide][jSol].resize(m_dim);
      }
    }
  
    for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
    {
      m_fluxJacobian[LEFT][iSol][iVar].resize(m_dim+m_ndimplus);
      m_fluxJacobian[RIGHT][iSol][iVar].resize(m_dim+m_ndimplus);
      m_gradientFluxJacobian[LEFT][iSol][iVar].resize(m_dim);
      m_gradientFluxJacobian[RIGHT][iSol][iVar].resize(m_dim);
      m_gradVarsToStateJacobian[LEFT][iSol][iVar].resize(m_nbrEqs);
      m_gradVarsToStateJacobian[RIGHT][iSol][iVar].resize(m_nbrEqs);

      for (CFuint iDim = 0; iDim < m_dim+m_ndimplus; ++iDim)
      {
        m_fluxJacobian[LEFT][iSol][iVar][iDim].resize(m_nbrEqs);
        m_fluxJacobian[RIGHT][iSol][iVar][iDim].resize(m_nbrEqs);
      }

      for (CFuint iDim = 0; iDim < m_dim; ++iDim)
      {
        m_gradientFluxJacobian[LEFT][iSol][iVar][iDim].resize(m_dim+m_ndimplus);
        m_gradientFluxJacobian[RIGHT][iSol][iVar][iDim].resize(m_dim+m_ndimplus);
        
        for (CFuint jDim = 0; jDim < m_dim+m_ndimplus; ++jDim)
        {
          m_gradientFluxJacobian[LEFT][iSol][iVar][iDim][jDim].resize(m_nbrEqs);
          m_gradientFluxJacobian[RIGHT][iSol][iVar][iDim][jDim].resize(m_nbrEqs); 
        }
      }
    }
  }
  
  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    m_riemannFluxJacobian[LEFT][iFlx].resize(m_nbrEqs);
    m_riemannFluxJacobian[RIGHT][iFlx].resize(m_nbrEqs);
    m_riemannFluxGradJacobian[iFlx].resize(m_nbrEqs);
    m_flxPntRiemannFluxDiff[iFlx].resize(m_nbrEqs);
    m_flxPntRiemannFluxPert[iFlx].resize(m_nbrEqs);
      
    for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
    {
      m_riemannFluxJacobian[LEFT][iFlx][iVar].resize(m_nbrEqs);
      m_riemannFluxJacobian[RIGHT][iFlx][iVar].resize(m_nbrEqs);
      m_riemannFluxGradJacobian[iFlx][iVar].resize(m_dim);
      
      for (CFuint iDim = 0; iDim < m_dim; ++iDim)
      {
        m_riemannFluxGradJacobian[iFlx][iVar][iDim].resize(m_nbrEqs);
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffJacobFluxReconstruction::unsetup()
{
  CFAUTOTRACE;
  
  // unsetup parent class
  DiffRHSJacobFluxReconstruction::unsetup();
}


//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD
