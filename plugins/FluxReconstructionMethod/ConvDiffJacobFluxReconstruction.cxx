// Copyright (C) 2016 KU Leuven, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Common/BadValueException.hh"
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
  CombinedJacobFluxReconstruction(name),
  m_elemIdx(),
  m_facesCell(),
  m_jacob(),
  socket_solPntNormals("solPntNormals"),
  socket_flxPntNormals("flxPntNormals"),
  socket_cellVolumes("cellVolumes"),
  socket_volumes("volumes"),
  m_tempSolPntVec(),
  m_tempSolPntVec2(),
  m_flxPntRiemannFluxPert(),
  m_tempFlux(),
  m_gradVarsToStateJacobian(),
  m_gradientStateJacobian(),
  m_contFlxBackupDiff()
  {
    addConfigOptionsTo(this);
  }

//////////////////////////////////////////////////////////////////////////////

void ConvDiffJacobFluxReconstruction::configure ( Config::ConfigArgs& args )
{
  CombinedJacobFluxReconstruction::configure(args);
}  

//////////////////////////////////////////////////////////////////////////////

void ConvDiffJacobFluxReconstruction::defineConfigOptions(Config::OptionList& options)
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector< Common::SafePtr< BaseDataSocketSource > >
  ConvDiffJacobFluxReconstruction::providesSockets()
{
  std::vector< Common::SafePtr< BaseDataSocketSource > > result = CombinedJacobFluxReconstruction::providesSockets();
  return result;
}

//////////////////////////////////////////////////////////////////////////////

std::vector< Common::SafePtr< BaseDataSocketSink > >
ConvDiffJacobFluxReconstruction::needsSockets()
{
  std::vector< Common::SafePtr< BaseDataSocketSink > > result = CombinedJacobFluxReconstruction::needsSockets();
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

      // Reset the value of m_nbrFaceFlxPnts in case it is not the same for all faces (Prism)
      m_nbrFaceFlxPnts = (*m_faceFlxPntConnPerOrient)[m_orient][0].size();

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

      // Reset the value of m_nbrFaceFlxPnts in case it is not the same for all faces (Prism)
      m_nbrFaceFlxPnts = (*m_faceFlxPntConnPerOrient)[m_orient][0].size();

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

	// fill the per-side cell metrics the compact face gradient needs
        prepareFaceCellMetrics();

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
        
        if (getMethodData().doComputeJacobian() && (!getMethodData().freezeJacob() || iter < iterFreeze || interval % getMethodData().getFreezeJacobInterval() == 0))
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

  // residual and Jacobian of the cells whose faces are all boundary faces
  computeCellsWithoutInnerFace();
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
      m_faceJacobVecs[iFlxPnt][iDim] = flxPntNormals[m_face->getID()*m_nbFaceFlxPntsMax*m_dim+iFlxPnt*m_dim+iDim];
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
    m_unitNormalFlxPnts[iFlxPnt] = m_mappedFaceNormalDir*m_faceJacobVecs[iFlxPnt]/m_faceJacobVecAbsSizeFlxPnts[iFlxPnt];
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
    m_nbrSolDep = ((*m_flxSolDep)[flxPntIdxL]).size();
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
      m_faceJacobVecs[iFlxPnt][iDim] = flxPntNormals[m_face->getID()*m_nbFaceFlxPntsMax*m_dim+iFlxPnt*m_dim+iDim];
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
    m_unitNormalFlxPnts[iFlxPnt] = m_mappedFaceNormalDir*m_faceJacobVecs[iFlxPnt]/m_faceJacobVecAbsSizeFlxPnts[iFlxPnt];
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
    m_nbrSolDep = ((*m_flxSolDep)[flxPntIdxL]).size();
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
  // common diffusive flux with the compact face gradients
  DiffRHSFluxReconstruction::computeInterfaceFlxCorrection();

  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    m_flxPntRiemannFluxDiff[iFlx] = m_flxPntRiemannFlux[iFlx];

    // subtract the convective Riemann flux
    m_flxPntRiemannFlux[iFlx] -= m_riemannFluxComputer->computeFlux(*m_cellStatesFlxPnt[LEFT][iFlx],*m_cellStatesFlxPnt[RIGHT][iFlx],m_unitNormalFlxPnts[iFlx]);

    // compute FI in the mapped coord frame
    for (CFuint iSide = 0; iSide < 2; ++iSide)
    {
      m_cellFlx[iSide][iFlx] = m_flxPntRiemannFlux[iFlx]*m_faceJacobVecSizeFlxPnts[iFlx][iSide];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffJacobFluxReconstruction::computeWaveSpeedUpdates(vector< CFreal >& waveSpeedUpd)
{
  // compute the wave speed updates for the neighbouring cells
  cf_assert(waveSpeedUpd.size() == 2);

  // get the correct m_faceIntegrationCoefs depending on the face type (only applicable for Prism for now) @todo should be updated for hybrid grid
  if (m_dim>2)
  {
    // get face geo
    const CFGeoShape::Type geo = m_face->getShape(); 

    if (geo == CFGeoShape::TRIAG) // triag face
    {
      //(*m_faceIntegrationCoefs).resize(m_nbrFaceFlxPnts);
      (m_faceIntegrationCoefs) = &(*m_faceIntegrationCoefsPerType)[0];
    }
    else  // quad face
    {
      //(*m_faceIntegrationCoefs).resize(m_nbrFaceFlxPnts);
      (m_faceIntegrationCoefs) = &(*m_faceIntegrationCoefsPerType)[1];
    } 
  }
  
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

      prepareSolPntFluxComputation((*(m_states[m_pertSide]))[m_pertSol]->getLocalID());

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
        prepareSolPntFluxComputation((*(m_states[m_pertSide]))[m_pertSol]->getLocalID());

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
    
        prepareFlxPntFluxComputation(iFlxPnt);
     
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
          prepareSolPntFluxComputation((*(m_states[m_pertSide]))[m_pertSol]->getLocalID());

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

        prepareFlxPntFluxComputation(iFlxPnt);
     
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
            m_gradientStateJacobian[m_pertSide][jSolIdx][m_pertSide][m_pertSol][iDim] += (m_neighbCellFluxProjVects[m_pertSide][jDim][m_pertSol][iDim]) * (*m_solPolyDerivAtSolPnts)[jSolIdx][jDim][m_pertSol];
          }
        }
      }
    }
    
    for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
    {
      // local flux point indices in the left and right cell
      const CFuint flxPntIdxThis = (*m_faceFlxPntConnPerOrient)[m_orient][m_pertSide][iFlx];
      const CFuint flxPntIdxOther = (*m_faceFlxPntConnPerOrient)[m_orient][iOtherSide][iFlx];

      m_nbrSolDep = ((*m_flxSolDep)[flxPntIdxThis]).size();
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
                    (*m_faceMappedCoordDir)[m_orient][m_pertSide]*m_mappedFaceNormalDir*m_faceJacobVecs[iFlx][iDim];
         
            m_gradientStateJacobian[iOtherSide][jSolIdxOther][m_pertSide][pertSolIdx][iDim] += 0.5 * m_corrFctDiv[jSolIdxOther][flxPntIdxOther] * (*m_solPolyValsAtFlxPnts)[flxPntIdxThis][pertSolIdx] *
                    (*m_faceMappedCoordDir)[m_orient][iOtherSide]*m_mappedFaceNormalDir*m_faceJacobVecs[iFlx][iDim];
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
      
      // Reset the value of m_nbrFaceFlxPnts in case it is not the same for all faces (Prism) //@todo
      m_nbrFaceFlxPnts = (*m_faceFlxPntConn)[faceIdx].size();

      if ((*m_isFaceOnBoundary[m_pertSide])[faceIdx])
      {  
        for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
        {
          const CFuint currFlxIdx = (*m_faceFlxPntConn)[faceIdx][iFlxPnt];
          
          m_nbrSolDep = ((*m_flxSolDep)[currFlxIdx]).size();
          for (m_pertSol = 0; m_pertSol < m_nbrSolDep; ++m_pertSol)
          {
            const CFuint pertSolIdx = (*m_flxSolDep)[currFlxIdx][m_pertSol];        
        
            for (CFuint jSol = 0; jSol < m_nbrSolDep; ++jSol)
            {
              const CFuint jSolIdx = (*m_flxSolDep)[currFlxIdx][jSol]; 
    
              for (CFuint iDim = 0; iDim < m_dim; ++iDim)
              {
                m_gradientStateJacobian[m_pertSide][jSolIdx][m_pertSide][pertSolIdx][iDim] -= 0.5 * m_corrFctDiv[jSolIdx][currFlxIdx] * (*m_solPolyValsAtFlxPnts)[currFlxIdx][pertSolIdx] *
                          (*m_faceLocalDir)[faceIdx]*flxPntNormals[faceID*m_nbFaceFlxPntsMax*m_dim+iFlxPnt*m_dim+iDim];
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
          
          m_nbrSolDep = ((*m_flxSolDep)[currFlxIdx]).size();
          for (m_pertSol = 0; m_pertSol < m_nbrSolDep; ++m_pertSol)
          {
            const CFuint pertSolIdx = (*m_flxSolDep)[currFlxIdx][m_pertSol];        
        
            for (CFuint jSol = 0; jSol < m_nbrSolDep; ++jSol)
            {
              const CFuint jSolIdx = (*m_flxSolDep)[currFlxIdx][jSol]; 
    
              for (CFuint iDim = 0; iDim < m_dim; ++iDim)
              {
                m_gradientStateJacobian[m_pertSide][jSolIdx][m_pertSide][pertSolIdx][iDim] -= 0.5 * m_corrFctDiv[jSolIdx][currFlxIdx] * (*m_solPolyValsAtFlxPnts)[currFlxIdx][pertSolIdx] *
                          ((*m_faceMappedCoordDir)[orient][cellSide])*flxPntNormals[faceID*m_nbFaceFlxPntsMax*m_dim+iFlxPnt*m_dim+iDim];
              }
            }
          }
        } 
      }
    }

    // Reset the value of m_nbrFaceFlxPnts in case it is not the same for all faces (Prism)
      m_nbrFaceFlxPnts = (*m_faceFlxPntConnPerOrient)[m_orient][0].size();
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

    prepareSolPntFluxComputation((*(m_states[side]))[iSolPnt]->getLocalID());

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
          m_unpertCellDiffRes[side][m_nbrEqs*iSolPnt+iEq] += polyCoef*(m_contFlx[jSolIdx][iDir][iEq]);
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
  addGradientFaceCorrections();
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffJacobFluxReconstruction::computeGradients()
{
  addGradientVolumeTerm();
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
  CombinedJacobFluxReconstruction::setup();

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
  CombinedJacobFluxReconstruction::unsetup();
}


//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD
