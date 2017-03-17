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

#include "FluxReconstructionMethod/DiffRHSFluxReconstruction.hh"
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
  
DiffRHSFluxReconstruction::DiffRHSFluxReconstruction(const std::string& name) :
  FluxReconstructionSolverCom(name),
  socket_gradients("gradients"),
  socket_normals("normals"),
  socket_rhs("rhs"),
  socket_updateCoeff("updateCoeff"),
  socket_faceJacobVecSizeFaceFlxPnts("faceJacobVecSizeFaceFlxPnts"),
  m_diffusiveVarSet(CFNULL),
  m_cellBuilder(CFNULL),
  m_faceBuilder(CFNULL),
  m_solPntsLocalCoords(CFNULL),
  m_faceIntegrationCoefs(CFNULL),
  m_faceMappedCoordDir(CFNULL),
  m_iElemType(),
  m_cell(),
  m_cellStates(),
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
  m_waveSpeedUpd(),
  m_faceLocalDir(),
  m_nbrFlxPnts(),
  m_nbrSolPnts(),
  m_cellGrads(),
  m_cellGradFlxPnt(),
  m_cflConvDiffRatio(),
  m_cellVolume(),
  m_nbrFaceFlxPnts()
  {
    addConfigOptionsTo(this);
  }
  
  
//////////////////////////////////////////////////////////////////////////////

void DiffRHSFluxReconstruction::defineConfigOptions(Config::OptionList& options)
{
}

//////////////////////////////////////////////////////////////////////////////

void DiffRHSFluxReconstruction::configure ( Config::ConfigArgs& args )
{
  FluxReconstructionSolverCom::configure(args);
}  
  
//////////////////////////////////////////////////////////////////////////////

std::vector< Common::SafePtr< BaseDataSocketSink > >
DiffRHSFluxReconstruction::needsSockets()
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

void DiffRHSFluxReconstruction::execute()
{
  CFAUTOTRACE;
  
  CFLog(VERBOSE, "DiffRHSFluxReconstruction::execute()\n");
  
  // get the elementTypeData
  SafePtr< vector<ElementTypeData> > elemType = MeshDataStack::getActive()->getElementTypeData();

  // get InnerCells TopologicalRegionSet
  SafePtr<TopologicalRegionSet> cells = MeshDataStack::getActive()->getTrs("InnerCells");

  // get the geodata of the geometric entity builder and set the TRS
  StdTrsGeoBuilder::GeoData& geoData = m_cellBuilder->getDataGE();
  geoData.trs = cells;
  
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
  
  //// Loop over faces to calculate fluxes and interface fluxes in the flux points
  
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
      geoData2.idx = faceID;
      m_face = m_faceBuilder->buildGE();
      
      // get the neighbouring cells
      m_cells[LEFT ] = m_face->getNeighborGeo(LEFT );
      m_cells[RIGHT] = m_face->getNeighborGeo(RIGHT);

      // get the states in the neighbouring cells
      m_states[LEFT ] = m_cells[LEFT ]->getStates();
      m_states[RIGHT] = m_cells[RIGHT]->getStates();

      setFaceData(m_face->getID());//faceID
      
      // compute volume
      m_cellVolume[LEFT] = m_cells[LEFT]->computeVolume();
      m_cellVolume[RIGHT] = m_cells[RIGHT]->computeVolume();
      
      // if one of the neighbouring cells is parallel updatable, compute the correction flux
      if ((*m_states[LEFT ])[0]->isParUpdatable() || (*m_states[RIGHT])[0]->isParUpdatable())
      {
	// compute the discontinuous fluxes
	computeDiscontinuousFluxes();
	
	// compute FI-FD
	computeInterfaceFlxCorrection();

	// compute the wave speed updates
        computeWaveSpeedUpdates(m_waveSpeedUpd);

        // update the wave speed
        updateWaveSpeed();

	// compute the correction for the left neighbour
	computeCorrection(LEFT, m_divContFlx);

	// update RHS
	updateRHS();
	
	// compute the correction for the right neighbour
	computeCorrection(RIGHT, m_divContFlx);
	
	// update RHS
	updateRHS();
      }

      // release the GeometricEntity
      m_faceBuilder->releaseGE();

    }
  }
  
  // add the correction due to partition faces
  addPartitionFacesCorrection();
  
  //// Loop over the elements to calculate the divergence of the continuous flux
  
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
      geoData.idx = elemIdx;
      m_cell = m_cellBuilder->buildGE();

      // get the states in this cell
      m_cellStates = m_cell->getStates();
      
      // if the states in the cell are parallel updatable, compute the resUpdates (-divFC)
      if ((*m_cellStates)[0]->isParUpdatable())
      {
	// set the cell data
	setCellData();

	// compute the residual updates (-divFC)
	computeResUpdates(m_divContFlx);
      
	// update RHS
        updateRHS();

      } 
      
      // divide by the Jacobian to transform the residuals back to the physical domain
      //divideByJacobDet();
      
      // print out the residual updates for debugging
      if(m_cell->getID() == 49)
      {
	CFLog(VERBOSE, "ID  = " << (*m_cellStates)[0]->getLocalID() << "\n");
        CFLog(VERBOSE, "Update = \n");
        // get the datahandle of the rhs
        DataHandle< CFreal > rhs = socket_rhs.getDataHandle();
        for (CFuint iState = 0; iState < m_nbrSolPnts; ++iState)
        {
          CFuint resID = m_nbrEqs*( (*m_cellStates)[iState]->getLocalID() );
          for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
          {
            CFLog(VERBOSE, "" << rhs[resID+iVar] << " ");
          }
          CFLog(VERBOSE,"\n");
          DataHandle<CFreal> updateCoeff = socket_updateCoeff.getDataHandle();
          CFLog(VERBOSE, "UpdateCoeff: " << updateCoeff[(*m_cellStates)[iState]->getLocalID()] << "\n");
        }
      }
      //release the GeometricEntity
      m_cellBuilder->releaseGE();
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void DiffRHSFluxReconstruction::computeInterfaceFlxCorrection()
{
      
  // Loop over the flux points to calculate FI-FD
  for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
  {
    RealVector avgSol;
    avgSol.resize(m_nbrEqs);
    vector< RealVector* > avgGrad;
    avgGrad.resize(m_nbrEqs);
     
    for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
    {
      avgGrad[iVar] = new RealVector(m_dim);
      *(avgGrad[iVar]) = (*(m_cellGradFlxPnt[LEFT][iFlxPnt][iVar]) + *(m_cellGradFlxPnt[RIGHT][iFlxPnt][iVar]))/2.0;
       
      avgSol[iVar] = ((*(m_cellStatesFlxPnt[LEFT][iFlxPnt]))[iVar] + (*(m_cellStatesFlxPnt[RIGHT][iFlxPnt]))[iVar])/2.0; 
    }
     
    m_flxPntRiemannFlux[iFlxPnt] = m_diffusiveVarSet->getFlux(avgSol,avgGrad,m_unitNormalFlxPnts[iFlxPnt],0);
     
    // compute FI-FD in the mapped coord frame
    m_cellFlx[LEFT][iFlxPnt] = (m_flxPntRiemannFlux[iFlxPnt] - m_cellFlx[LEFT][iFlxPnt])*m_faceJacobVecSizeFlxPnts[iFlxPnt][LEFT];
    m_cellFlx[RIGHT][iFlxPnt] = (m_flxPntRiemannFlux[iFlxPnt] - m_cellFlx[RIGHT][iFlxPnt])*m_faceJacobVecSizeFlxPnts[iFlxPnt][RIGHT];

    if(m_cells[LEFT]->getID() == 150 || m_cells[RIGHT]->getID() == 150)
    {
      CFLog(VERBOSE, "faceID =  " << m_face->getID() << "\n");
      CFLog(VERBOSE, "cell 0 == LEFT <=> " << (m_cells[LEFT]->getID() == 281) << "\n");
      CFLog(VERBOSE, "Unit vector = " << m_unitNormalFlxPnts[iFlxPnt] << "\n");
      CFLog(VERBOSE, "flxIdx LEFT = " << (*m_faceFlxPntConnPerOrient)[m_orient][LEFT][iFlxPnt] << "\n");
      CFLog(VERBOSE, "flxIdx RIGHT = " << (*m_faceFlxPntConnPerOrient)[m_orient][RIGHT][iFlxPnt] << "\n");
      RealVector RiemannL = (m_flxPntRiemannFlux[iFlxPnt]*m_faceJacobVecSizeFlxPnts[iFlxPnt][LEFT]);
      RealVector RiemannR = (m_flxPntRiemannFlux[iFlxPnt]*m_faceJacobVecSizeFlxPnts[iFlxPnt][RIGHT]);
      CFLog(VERBOSE,"RiemannL = " << RiemannL << "\n");
      CFLog(VERBOSE,"RiemannR = " << RiemannR << "\n");
      CFLog(VERBOSE,"delta fluxL = " << m_cellFlx[LEFT][iFlxPnt] << "\n");
      CFLog(VERBOSE,"delta fluxR = " << m_cellFlx[RIGHT][iFlxPnt] << "\n");
      CFLog(VERBOSE,"Jacob L = " << m_faceJacobVecSizeFlxPnts[iFlxPnt][LEFT] << "\n");
      CFLog(VERBOSE,"Jacob R = " << m_faceJacobVecSizeFlxPnts[iFlxPnt][RIGHT] << "\n");
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void DiffRHSFluxReconstruction::setFaceData(CFuint faceID)
{
  // get the local FR data
  vector< FluxReconstructionElementData* >& frLocalData = getMethodData().getFRLocalData();

  // compute flux point coordinates
  SafePtr< vector<RealVector> > flxLocalCoords = frLocalData[0]->getFaceFlxPntsFaceLocalCoords();
    
  // compute face Jacobian vectors
  vector< RealVector > faceJacobVecs = m_face->computeFaceJacobDetVectorAtMappedCoords(*flxLocalCoords);

  // Loop over flux points to set the normal vectors
  for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
  {
    // get face Jacobian vector sizes in the flux points
    DataHandle< vector< CFreal > >
    faceJacobVecSizeFaceFlxPnts = socket_faceJacobVecSizeFaceFlxPnts.getDataHandle();
  
    // get face Jacobian vector size
    m_faceJacobVecAbsSizeFlxPnts[iFlxPnt] = faceJacobVecSizeFaceFlxPnts[faceID][iFlxPnt];

    // set face Jacobian vector size with sign depending on mapped coordinate direction
    m_faceJacobVecSizeFlxPnts[iFlxPnt][LEFT] = m_faceJacobVecAbsSizeFlxPnts[iFlxPnt]*(*m_faceMappedCoordDir)[m_orient][LEFT];
    m_faceJacobVecSizeFlxPnts[iFlxPnt][RIGHT] = m_faceJacobVecAbsSizeFlxPnts[iFlxPnt]*(*m_faceMappedCoordDir)[m_orient][RIGHT];

    // set unit normal vector
    m_unitNormalFlxPnts[iFlxPnt] = faceJacobVecs[iFlxPnt]/m_faceJacobVecAbsSizeFlxPnts[iFlxPnt];
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

void DiffRHSFluxReconstruction::computeDiscontinuousFluxes()
{
  // get the local FR data
  vector< FluxReconstructionElementData* >& frLocalData = getMethodData().getFRLocalData();
  
  // get solution polynomial values at nodes
  vector< vector< CFreal > > solPolyValsAtNodes
       = frLocalData[0]->getSolPolyValsAtNode(*m_flxPntsLocalCoords);

  // loop over flx pnts to extrapolate the states to the flux points and get the 
  // discontinuous normal flux in the flux points and the Riemann flux
  for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
  {     
    // local flux point indices in the left and right cell
    CFuint flxPntIdxL = (*m_faceFlxPntConnPerOrient)[m_orient][LEFT][iFlxPnt];
    CFuint flxPntIdxR = (*m_faceFlxPntConnPerOrient)[m_orient][RIGHT][iFlxPnt]; 
    
    // Loop over solution points to extrapolate the state to the flux points and calculate the discontinuous flux
    State* tempVectorL = new State(solPolyValsAtNodes[flxPntIdxL][0]*(*((*(m_states[LEFT]))[0]->getData())));
    State* tempVectorR = new State(solPolyValsAtNodes[flxPntIdxR][0]*(*((*(m_states[RIGHT]))[0]->getData())));

    m_cellStatesFlxPnt[LEFT][iFlxPnt] = tempVectorL;
    m_cellStatesFlxPnt[RIGHT][iFlxPnt] = tempVectorR;

    for (CFuint iSol = 1; iSol < m_nbrSolPnts; ++iSol)
    {
      //CFLog(VERBOSE, "Left State GlobalID: " << (*(m_states[LEFT]))[iSol]->getGlobalID() << "\n");
      if(m_face->getID() == 624)
      {
	CFLog(DEBUG_MIN, "Left State " << *((*(m_states[LEFT]))[iSol]->getData()) << "\n");
	CFLog(DEBUG_MIN, "Right State = " << *((*(m_states[RIGHT]))[iSol]->getData()) << "\n");
      }
      *(m_cellStatesFlxPnt[LEFT][iFlxPnt]) = (State) (*(m_cellStatesFlxPnt[LEFT][iFlxPnt]->getData()) + 
							  solPolyValsAtNodes[flxPntIdxL][iSol]*(*((*(m_states[LEFT]))[iSol]->getData())));
      *(m_cellStatesFlxPnt[RIGHT][iFlxPnt]) = (State) (*(m_cellStatesFlxPnt[RIGHT][iFlxPnt]->getData()) +
							   solPolyValsAtNodes[flxPntIdxR][iSol]*(*((*(m_states[RIGHT]))[iSol]->getData())));
    }
    if(m_face->getID() == 624)
    {
      CFLog(DEBUG_MIN, "cellID = " << m_cells[LEFT]->getID() << " or " << m_cells[RIGHT]->getID() << "\n");
      CFLog(DEBUG_MIN, "left state in flx pnt = " << *(m_cellStatesFlxPnt[LEFT][iFlxPnt]->getData()) << "\n");
      CFLog(DEBUG_MIN, "right state in flx pnt = " << *(m_cellStatesFlxPnt[RIGHT][iFlxPnt]->getData()) << "\n");
      CFLog(DEBUG_MIN, "Normal: " << m_unitNormalFlxPnts[iFlxPnt] << "\n");
    }
    
    // Loop over flux points to extrapolate the gradients to the flux points
    for (CFuint iSide = 0; iSide < 2; ++iSide)
    {
      for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
      {
        for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
        {
          (m_cellGradFlxPnt[iSide][iFlxPnt][iVar]) = new RealVector(m_dim);
	  *(m_cellGradFlxPnt[iSide][iFlxPnt][iVar]) = 0.0;
        }
        
        for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
        {
          for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
          {
            CFLog(DEBUG_MIN,"cellStates: " << *((*m_cellStates)[iSol]->getData()) << "\n");
            *(m_cellGradFlxPnt[iSide][iFlxPnt][iVar]) += solPolyValsAtNodes[iFlxPnt][iSol]*((*(m_cellGrads[iSide][iSol]))[iVar]);
          }
        }
      }
    }

    // compute the normal flux at the current flux point
    m_cellFlx[LEFT][iFlxPnt] = m_diffusiveVarSet->getFlux(*(m_cellStatesFlxPnt[LEFT][iFlxPnt]->getData()),m_cellGradFlxPnt[LEFT][iFlxPnt],m_unitNormalFlxPnts[iFlxPnt],0);
    m_cellFlx[RIGHT][iFlxPnt] = m_diffusiveVarSet->getFlux(*(m_cellStatesFlxPnt[RIGHT][iFlxPnt]->getData()),m_cellGradFlxPnt[RIGHT][iFlxPnt],m_unitNormalFlxPnts[iFlxPnt],0);

    if(m_cells[LEFT]->getID() == 49 || m_cells[RIGHT]->getID() == 49)
    {
      CFLog(DEBUG_MIN, "Flux Left = " << m_cellFlx[LEFT][iFlxPnt] << "\n");
      CFLog(DEBUG_MIN, "Flux Right = " << m_cellFlx[RIGHT][iFlxPnt] << "\n");
      RealVector fluxL = (m_cellFlx[LEFT][iFlxPnt]*m_faceJacobVecSizeFlxPnts[iFlxPnt][LEFT]);
      RealVector fluxR = (m_cellFlx[RIGHT][iFlxPnt]*m_faceJacobVecSizeFlxPnts[iFlxPnt][RIGHT]);
      CFLog(VERBOSE,"fluxLLocal = " << fluxL << "\n");
      CFLog(VERBOSE,"fluxRLocal = " << fluxR << "\n");
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void DiffRHSFluxReconstruction::computeResUpdates(vector< RealVector >& residuals)
{
  // get the local FR data
  vector< FluxReconstructionElementData* >& frLocalData = getMethodData().getFRLocalData();
  
  // get solution derivative polynomial values at solution points
  vector< vector< vector< CFreal > > > solDerivPolyValsAtNodes
         = frLocalData[m_iElemType]-> getSolPolyDerivsAtNode(*m_solPntsLocalCoords);
   
  // Loop over solution points to calculate the discontinuous flux.
  for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolPnts; ++iSolPnt)
  {
    // calculate the discontinuous flux projected on x, y, z-directions
    for (CFuint iDim = 0; iDim < m_dim; ++iDim)
    {
      // dereference the state
      State& stateSolPnt = *(*m_cellStates)[iSolPnt];
      
      vector< RealVector > temp = *(m_cellGrads[0][iSolPnt]);
      vector< RealVector* > grad;
      grad.resize(m_nbrEqs);
      
      for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
      {
	cf_assert(temp.size() == m_nbrEqs);
	grad[iVar] = & (temp[iVar]);
      }
      
      m_contFlx[iSolPnt][iDim] = m_diffusiveVarSet->getFlux(*(stateSolPnt.getData()),grad,m_cellFluxProjVects[iDim][iSolPnt],0);
    }
  }
         
  // Loop over solution pnts to calculate the divergence of the discontinuous flux
  for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolPnts; ++iSolPnt)
  {
    // reset the divergence of FC
    residuals[iSolPnt] = 0.0;
    // Loop over solution pnt to count factor of all sol pnt polys
    for (CFuint jSolPnt = 0; jSolPnt < m_nbrSolPnts; ++jSolPnt)
    {
      // Loop over deriv directions and sum them to compute divergence
      for (CFuint iDir = 0; iDir < m_dim; ++iDir)
      {
        // Loop over conservative fluxes 
        for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
        {
          // Store divFD in the vector that will be divFC
          residuals[iSolPnt][iEq] -= solDerivPolyValsAtNodes[iSolPnt][iDir][jSolPnt]*(m_contFlx[jSolPnt][iDir][iEq]);
       
	  if (fabs(residuals[iSolPnt][iEq]) < MathTools::MathConsts::CFrealEps())
          {
            residuals[iSolPnt][iEq] = 0;
	  }
	}
      }
    }
    if(m_cell->getID() == 49)
    {
      CFLog(VERBOSE, "state: " << *((*m_cellStates)[iSolPnt]->getData()) << "\n");
      CFLog(VERBOSE, "ID: " << m_cell->getID() << "\n");
      CFLog(VERBOSE, "flx in " << iSolPnt << " : (" << m_contFlx[iSolPnt][0] << " , " << m_contFlx[iSolPnt][1] << "\n");
      CFLog(VERBOSE, "-div FD = " << residuals[iSolPnt] << "\n");
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void DiffRHSFluxReconstruction::setCellData()
{
  // create a list of the dimensions in which the deriv will be calculated
  for (CFuint iDim = 0; iDim < m_dim; ++iDim)
  {
    vector<CFuint> dimList;
    dimList.resize(m_nbrSolPnts);
    for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolPnts; ++iSolPnt)
    {
      dimList[iSolPnt] = iDim;
    }
    m_cellFluxProjVects[iDim] = m_cell->computeMappedCoordPlaneNormalAtMappedCoords(dimList,
                                                                            *m_solPntsLocalCoords);
    if(m_cell->getID() == 49)
    {
      for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
      {
        CFLog(DEBUG_MIN, "normal along " << iDim << " for sol pnt " << iSol << " : " << m_cellFluxProjVects[iDim][iSol] << "\n");
      }
    }
  }
  
  // get the gradients datahandle
  DataHandle< vector< RealVector > > gradients = socket_gradients.getDataHandle();

  for (CFuint iState = 0; iState < m_nbrSolPnts; ++iState)
  {
    const CFuint stateID = (*(m_cellStates))[iState]->getLocalID();
    m_cellGrads[0][iState] = &gradients[stateID];
  }
}

//////////////////////////////////////////////////////////////////////////////

void DiffRHSFluxReconstruction::updateRHS()
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
      rhs[resID+iVar] += resFactor*m_divContFlx[iState][iVar];
    }
  }
}


//////////////////////////////////////////////////////////////////////////////

void DiffRHSFluxReconstruction::updateWaveSpeed()
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
      CFLog(DEBUG_MIN, "updateCoeff = " << updateCoeff[solID] << "\n");
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void DiffRHSFluxReconstruction::computeCorrection(CFuint side, vector< RealVector >& corrections)
{ 
  cf_assert(corrections.size() == m_nbrSolPnts);
  
  // loop over sol pnts to compute the corrections
  for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolPnts; ++iSolPnt)
  {
    // reset the corrections which will be stored in divContFlx in order to be able to reuse updateRHS() 
    corrections[iSolPnt] = 0.0;
    
    cf_assert(corrections[iSolPnt].size() == m_nbrEqs);

    // compute the term due to each flx pnt
    for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
    {
      // the current correction factor (stored in cellFlx)
      RealVector currentCorrFactor = m_cellFlx[side][iFlxPnt];
      cf_assert(currentCorrFactor.size() == m_nbrEqs);
    
      // Fill in the corrections
      for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
      {
        corrections[iSolPnt][iVar] -= currentCorrFactor[iVar] * m_corrFctDiv[iSolPnt][(*m_faceFlxPntConnPerOrient)[m_orient][side][iFlxPnt]]; 
      }
      
      if(m_cells[side]->getID() == 49)
      {
        CFLog(VERBOSE, "FI-FD = " << currentCorrFactor << "\n");
        CFLog(VERBOSE, "iSol: " << iSolPnt << ", flxID = " << (*m_faceFlxPntConnPerOrient)[m_orient][side][iFlxPnt] << "\n");
        CFLog(VERBOSE, "div h = " << m_corrFctDiv[iSolPnt][(*m_faceFlxPntConnPerOrient)[m_orient][side][iFlxPnt]] << "\n");
      }
    }
    
    if(m_cells[side]->getID() == 49)
    {
      CFLog(VERBOSE, "correction = " << corrections[iSolPnt] << "\n");
    }
  }
  
  // in order to use updateRHS, m_cellStates should have the correct states
  m_cellStates = m_cells[side]->getStates();
}

//////////////////////////////////////////////////////////////////////////////

void DiffRHSFluxReconstruction::divideByJacobDet()
{
  // this is achieved by multiplying the update coefs with the Jacobian determinant
  // (and dividing by the cell volume)

  // get the updateCoeff
  DataHandle< CFreal > updateCoeff = socket_updateCoeff.getDataHandle();

  // get the cell volume
  const CFreal invCellVolume = 1.0/m_cell->computeVolume();

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
    updateCoeff[solID] *= jacobDet[iSol]*invCellVolume;
  }
  
}

//////////////////////////////////////////////////////////////////////////////

void DiffRHSFluxReconstruction::setup()
{
  CFAUTOTRACE;
  FluxReconstructionSolverCom::setup();
  
  // boolean telling whether there is a diffusive term
  const bool hasDiffTerm = getMethodData().hasDiffTerm();

  if (hasDiffTerm)
  {
    // get the diffusive varset
    m_diffusiveVarSet = getMethodData().getDiffusiveVar(); 
  }
  
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
  // for now, there should be only one type of element
  cf_assert(frLocalData.size() == 1);
  
  // compute flux point coordinates
  SafePtr< vector<RealVector> > flxLocalCoords = frLocalData[0]->getFaceFlxPntsFaceLocalCoords();
  m_nbrFaceFlxPnts = flxLocalCoords->size();
  
  // number of flux points
  m_nbrFlxPnts = frLocalData[0]->getNbrOfFlxPnts();
  
  // number of sol points
  m_nbrSolPnts = frLocalData[0]->getNbrOfSolPnts();
  
  cf_assert(m_nbrFlxPnts == (frLocalData[0]->getFlxPntsLocalCoords())->size());
  cf_assert(m_nbrSolPnts == (frLocalData[0]->getSolPntsLocalCoords())->size());

  // dimensionality and number of equations
  m_dim   = PhysicalModelStack::getActive()->getDim();
  m_nbrEqs = PhysicalModelStack::getActive()->getNbEq();
  
  // get solution point local coordinates
  m_solPntsLocalCoords = frLocalData[0]->getSolPntsLocalCoords();
    
  // get flux point local coordinates
  m_flxPntsLocalCoords = frLocalData[0]->getFlxPntsLocalCoords();
   
  // get the face - flx pnt connectivity per orient
  m_faceFlxPntConnPerOrient = frLocalData[0]->getFaceFlxPntConnPerOrient();
    
  // get the face connectivity per orientation
  m_faceConnPerOrient = frLocalData[0]->getFaceConnPerOrient();
  
  // get the face integration coefficient
  m_faceIntegrationCoefs = frLocalData[0]->getFaceIntegrationCoefs();
  
  // get flux point mapped coordinate directions per orient
  m_faceMappedCoordDir = frLocalData[0]->getFaceMappedCoordDirPerOrient();
  
  // get flux point mapped coordinate directions
  m_faceLocalDir = frLocalData[0]->getFaceMappedCoordDir();
  
  // get flux point local coordinates
  m_flxPntsLocalCoords1D = frLocalData[0]->getFlxPntsLocalCoord1D();
    
  // get the face - flx pnt connectivity
  m_faceFlxPntConn = frLocalData[0]->getFaceFlxPntConn();
  
  // get convective/diffusive CFL ratio
  m_cflConvDiffRatio = frLocalData[0]->getCFLConvDiffRatio();
  
  // resize the physical data temporary vector
  SafePtr<BaseTerm> convTerm = PhysicalModelStack::getActive()->getImplementor()->getConvectiveTerm(); 
  convTerm->resizePhysicalData(m_pData);
  // Resize vectors
  m_cells.resize(2);
  m_states.resize(2);
  m_cellStatesFlxPnt.resize(2);
  m_cellFlx.resize(2);
  m_cellGrads.resize(2);
  m_cellGradFlxPnt.resize(2);
  m_cellVolume.resize(2);
  m_faceJacobVecAbsSizeFlxPnts.resize(m_nbrFaceFlxPnts);
  m_cellStatesFlxPnt[LEFT].resize(m_nbrFaceFlxPnts);
  m_cellStatesFlxPnt[RIGHT].resize(m_nbrFaceFlxPnts);
  m_cellGrads[LEFT].resize(m_nbrSolPnts);
  m_cellGrads[RIGHT].resize(m_nbrSolPnts);
  m_cellGradFlxPnt[LEFT].resize(m_nbrFaceFlxPnts);
  m_cellGradFlxPnt[RIGHT].resize(m_nbrFaceFlxPnts);
  m_cellFlx[LEFT].resize(m_nbrFaceFlxPnts);
  m_cellFlx[RIGHT].resize(m_nbrFaceFlxPnts);
  m_faceJacobVecSizeFlxPnts.resize(m_nbrFaceFlxPnts);
  m_unitNormalFlxPnts.resize(m_nbrFaceFlxPnts);
  m_flxPntRiemannFlux.resize(m_nbrFaceFlxPnts);
  m_contFlx.resize(m_nbrSolPnts);
  m_divContFlx.resize(m_nbrSolPnts);
  m_corrFctDiv.resize(m_nbrSolPnts);
  m_cellFluxProjVects.resize(m_dim);
  m_flxPntCoords.resize(m_nbrFaceFlxPnts);
  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    m_flxPntCoords[iFlx].resize(m_dim);
    m_faceJacobVecSizeFlxPnts[iFlx].resize(2);
    m_unitNormalFlxPnts[iFlx].resize(m_dim);
    m_cellFlx[LEFT][iFlx].resize(m_nbrEqs);
    m_cellFlx[RIGHT][iFlx].resize(m_nbrEqs);
    m_flxPntRiemannFlux[iFlx].resize(m_nbrEqs);
    m_cellGradFlxPnt[LEFT][iFlx].resize(m_nbrEqs);
    m_cellGradFlxPnt[RIGHT][iFlx].resize(m_nbrEqs);
  }

  for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolPnts; ++iSolPnt)
  {
    m_contFlx[iSolPnt].resize(m_dim);
    m_divContFlx[iSolPnt].resize(m_nbrEqs);
    m_corrFctDiv[iSolPnt].resize(m_nbrFlxPnts);
    for (CFuint iDim = 0; iDim < m_dim; ++iDim)
    {
      m_contFlx[iSolPnt][iDim].resize(m_nbrEqs);
    }
  }
  
  for (CFuint iDim = 0; iDim < m_dim; ++iDim)
  {
    m_cellFluxProjVects[iDim].resize(m_nbrSolPnts);
    for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolPnts; ++iSolPnt)
    {
      m_cellFluxProjVects[iDim][iSolPnt].resize(m_dim);
    }
  }
  
  // compute the divergence of the correction function
  m_corrFctComputer->computeDivCorrectionFunction(frLocalData[0],m_corrFctDiv);
}

//////////////////////////////////////////////////////////////////////////////

void DiffRHSFluxReconstruction::unsetup()
{
  CFAUTOTRACE;
  FluxReconstructionSolverCom::unsetup();
}


//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD

