// Copyright (C) 2016 KU Leuven, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/MethodCommandProvider.hh"

#include "Framework/CFSide.hh"
#include "Framework/MeshData.hh"
#include "Framework/BaseTerm.hh"

#include "MathTools/MathFunctions.hh"

#include "FluxReconstructionMethod/LLAVFluxReconstruction.hh"
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

MethodCommandProvider< LLAVFluxReconstruction,
		       FluxReconstructionSolverData,
		       FluxReconstructionModule >
LLAVFluxReconstructionFluxReconstructionProvider("LLAV");

//////////////////////////////////////////////////////////////////////////////
  
LLAVFluxReconstruction::LLAVFluxReconstruction(const std::string& name) :
  DiffRHSFluxReconstruction(name),
  m_updateVarSet(CFNULL),
  m_order(),
  m_transformationMatrix(),
  m_statesPMinOne(),
  m_epsilon(),
  m_solEpsilons(),
  m_epsilonLR(),
  m_epsilon0(),
  m_s0(),
  m_s0Prev(),
  m_s(),
  m_kappa(),
  m_showrate(),
  m_peclet(),
  m_nodeEpsilons(),
  m_nbNodeNeighbors(),
  m_cellEpsilons(),
  m_cellNodes(),
  m_faceNodes(),
  m_nbrCornerNodes(),
  m_flagComputeNbNghb(),
  m_nodePolyValsAtFlxPnts(),
  m_nodePolyValsAtSolPnts(),
  m_cellNodesConn(CFNULL),
  m_elemIdx(),
  m_cellBuilder(CFNULL),
  m_isFaceOnBoundaryCell(CFNULL),
  m_nghbrCellSideCell(CFNULL),
  m_currCellSideCell(CFNULL),
  m_faceOrientsCell(CFNULL),
  m_faceBCIdxCell(CFNULL),
  m_faces(),
  m_bcStateComputers(CFNULL),
  m_flxPntGhostGrads(),
  m_freezeLimiterRes(),
  m_freezeLimiterIter(),
  m_useMax(),
  m_totalEps(),
  m_totalEpsGlobal(),
  m_Smax(),
  m_SmaxGlobal(),
  m_jacob(false),
  m_nbPosPrev(),
  m_nbPosPrevGlobal(),
  m_subcellRes(),
  socket_artVisc("artVisc"),
  socket_monPhysVar("monPhysVar"),
  socket_smoothness("smoothness"),
  m_maxLambda(),
  m_unitNormalFlxPnts2(),
  m_faceJacobVecSizeFlxPnts2(),
  m_tempSolPntVec(),
  m_tempSolPntVec2()
  {
    addConfigOptionsTo(this);
    
    m_kappa = 5.0;
    setParameter( "Kappa", &m_kappa);
    
    m_showrate = 1;
    setParameter( "ShowRate", &m_showrate);
    
    m_peclet = 2.0;
    setParameter( "Peclet", &m_peclet);
    
    m_s0 = 0.0;
    setParameter( "S0", &m_s0);
    
    m_freezeLimiterRes = -20.0;
    setParameter( "FreezeLimiterRes", &m_freezeLimiterRes);
  
    m_freezeLimiterIter = MathTools::MathConsts::CFuintMax();
    setParameter( "FreezeLimiterIter", &m_freezeLimiterIter);
    
    m_freezeSmoothnessIter = MathTools::MathConsts::CFuintMax();
    setParameter( "FreezeSmoothnessIter", &m_freezeSmoothnessIter);
    
    m_addPosPrev = false;
    setParameter( "AddPositivityPreservation", &m_addPosPrev);
    
    m_minValue = 1.0e-12;
    setParameter( "MinValue", &m_minValue);
    cf_assert(m_minValue > 0.0);
    
    m_monitoredVar = 0;
    setParameter( "MonitoredVar", &m_monitoredVar);
    
    m_viscFactor = 2.0;
    setParameter( "ViscFactor", &m_viscFactor);
    
    m_addUpdCoeff = true;
    setParameter( "AddUpdateCoeff", &m_addUpdCoeff);
    
    m_monitoredPhysVar = MathTools::MathConsts::CFuintMax();
    setParameter( "MonitoredPhysVar", &m_monitoredPhysVar);
  }
  
  
//////////////////////////////////////////////////////////////////////////////

void LLAVFluxReconstruction::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal,Config::DynamicOption<> >("Kappa","Kappa factor of artificial viscosity.");
  
  options.addConfigOption< CFuint,Config::DynamicOption<> >("ShowRate","Showrate of the LLAV information.");
  
  options.addConfigOption< CFreal,Config::DynamicOption<> >("Peclet","Peclet number to be used for artificial viscosity.");
  
  options.addConfigOption< CFreal,Config::DynamicOption<> >("S0","Reference smoothness factor, will be multiplied by -log(P).");
  
  options.addConfigOption< CFreal,Config::DynamicOption<> >("FreezeLimiterRes","Residual after which to freeze the residual.");
  
  options.addConfigOption< CFuint,Config::DynamicOption<> >("FreezeLimiterIter","Iteration after which to freeze the residual.");
  
  options.addConfigOption< CFuint,Config::DynamicOption<> >("FreezeSmoothnessIter","Iteration after which to freeze the reference smoothness.");
  
  options.addConfigOption< bool >("AddPositivityPreservation","Bool telling whether extra viscosity needs to be added for positivity preservation.");
  
  options.addConfigOption< CFreal,Config::DynamicOption<> >("MinValue","Minimum value at which point positivity preservation is added.");
  
  options.addConfigOption< CFreal >("ViscFactor","Maximum factor applied to viscosity for positivity preservation.");
  
  options.addConfigOption< CFuint,Config::DynamicOption<> >("MonitoredVar","Index of the monitored var for positivity preservation.");
  
  options.addConfigOption< bool >("AddUpdateCoeff","Boolean telling whether the update coefficient based on the artificial flux is added.");
  
  options.addConfigOption< CFuint,Config::DynamicOption<> >("MonitoredPhysVar","Index of the monitored physical var for positivity preservation, if not specified MonitoredVar is used instead.");
}

//////////////////////////////////////////////////////////////////////////////

void LLAVFluxReconstruction::configure ( Config::ConfigArgs& args )
{
  FluxReconstructionSolverCom::configure(args);
}  

//////////////////////////////////////////////////////////////////////////////

std::vector< Common::SafePtr< BaseDataSocketSource > >
  LLAVFluxReconstruction::providesSockets()
{
  std::vector< Common::SafePtr< BaseDataSocketSource > > result;
  result.push_back(&socket_artVisc);
  result.push_back(&socket_monPhysVar);
  result.push_back(&socket_smoothness);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

void LLAVFluxReconstruction::execute()
{
  CFAUTOTRACE;
  
  CFLog(VERBOSE, "LLAVFluxReconstruction::execute()\n");
  
  // get the elementTypeData
  SafePtr< vector<ElementTypeData> > elemType = MeshDataStack::getActive()->getElementTypeData();

  // get InnerCells TopologicalRegionSet
  SafePtr<TopologicalRegionSet> cells = MeshDataStack::getActive()->getTrs("InnerCells");

  // get the geodata of the geometric entity builder and set the TRS
  CellToFaceGEBuilder::GeoData& geoDataCell = m_cellBuilder->getDataGE();
  geoDataCell.trs = cells;
  
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
  
  const CFreal residual = SubSystemStatusStack::getActive()->getResidual();
  
  const CFuint iter = SubSystemStatusStack::getActive()->getNbIter();
  
  m_useMax = residual < m_freezeLimiterRes || iter > m_freezeLimiterIter;
  m_totalEps = 0.0;
  
//  if (iter < m_freezeSmoothnessIter) 
//  {
//    m_s0 = 0.2*(m_Smax - m_kappa) + 0.8*m_s0;
//    m_s0Prev = m_s0;
//  }
//  else
//  {
//    m_s0 = m_s0Prev;
//  }
  //m_peclet = m_minValue/pow(10,m_s0);
  
  m_Smax = -100.0;
  
  m_nodeEpsilons = 0.0;
  
  m_nbPosPrev = 0;
  m_maxLambda = 0.;
  
  //// Loop over the elements to compute the artificial viscosities
  
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
      
      //CFLog(INFO, "state: " << *((*m_cellStates)[0]) << "\n");
      
      // get the nodes in this cell
      m_cellNodes  = m_cell->getNodes();
      
//       // if the states in the cell are parallel updatable, compute the resUpdates (-divFC)
//       if ((*m_cellStates)[0]->isParUpdatable())
//       {
      // compute the states projected on order P-1
      computeProjStates(m_statesPMinOne);
	
      //CFLog(INFO, "projstate: " << m_statesPMinOne[0] << "\n");
	
      // compute the artificial viscosity
      computeEpsilon();
	
      // store epsilon
      storeEpsilon();
//       } 
      
      //release the GeometricEntity
      m_cellBuilder->releaseGE();
    }
  }
  
  //CFLog(INFO, "max lambda: " << m_maxLambda << "\n");
  
  const std::string nsp = this->getMethodData().getNamespace();
  
#ifdef CF_HAVE_MPI
    MPI_Comm comm = PE::GetPE().GetCommunicator(nsp);
    PE::GetPE().setBarrier(nsp);
    const CFuint count = 1;
    MPI_Allreduce(&m_totalEps, &m_totalEpsGlobal, count, MPI_DOUBLE, MPI_SUM, comm);
    MPI_Allreduce(&m_Smax, &m_SmaxGlobal, count, MPI_DOUBLE, MPI_MAX, comm);
#endif
    
  if (PE::GetPE().GetRank(nsp) == 0 && iter%m_showrate == 0) 
  {
    // print total artificial viscosity
    CFLog(INFO, "total eps: " << m_totalEpsGlobal << ", Smax: " << m_SmaxGlobal << "\n");
  }

  PE::GetPE().setBarrier(nsp);
  
  m_Smax = m_SmaxGlobal;
  
  m_flagComputeNbNghb = false;
  
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
      geoDataFace.idx = faceID;
      m_face = m_faceBuilder->buildGE();
      
      // get the neighbouring cells
      m_cells[LEFT ] = m_face->getNeighborGeo(LEFT );
      m_cells[RIGHT] = m_face->getNeighborGeo(RIGHT);

      // get the states in the neighbouring cells
      m_states[LEFT ] = m_cells[LEFT ]->getStates();
      m_states[RIGHT] = m_cells[RIGHT]->getStates();
      
      // compute volume
      m_cellVolume[LEFT] = m_cells[LEFT]->computeVolume();
      m_cellVolume[RIGHT] = m_cells[RIGHT]->computeVolume();
      
      // if one of the neighbouring cells is parallel updatable, compute the correction flux
  //    if ((*m_states[LEFT ])[0]->isParUpdatable() || (*m_states[RIGHT])[0]->isParUpdatable())
  //    {
	// set the face data
	setFaceData(m_face->getID());//faceID
	
	// compute the left and right states and gradients in the flx pnts
	computeFlxPntStatesAndGrads();
	
	// compute the common interface flux
	computeInterfaceFlxCorrection();

	if (m_addUpdCoeff)
	{
	  // compute the wave speed updates
          computeWaveSpeedUpdates(m_waveSpeedUpd);

          // update the wave speed
          updateWaveSpeed();
	}

	// compute the correction for the left neighbour
	computeCorrection(LEFT, m_divContFlx);

	// update RHS
	updateRHS();
	
	// compute the correction for the right neighbour
	computeCorrection(RIGHT, m_divContFlx);
	
	// update RHS
	updateRHS();
   //   }

      // release the GeometricEntity
      m_faceBuilder->releaseGE();
    }
  }

  //// Loop over the elements to calculate the divergence of the continuous flux
  
  // loop over element types, for the moment there should only be one
  for (m_iElemType = 0; m_iElemType < nbrElemTypes; ++m_iElemType)
  {
    // get start and end indexes for this type of element
    const CFuint startIdx = (*elemType)[m_iElemType].getStartIdx();
    const CFuint endIdx   = (*elemType)[m_iElemType].getEndIdx();

    // loop over cells
    for (m_elemIdx = startIdx; m_elemIdx < endIdx; ++m_elemIdx)
    {
      // build the GeometricEntity
      geoDataCell.idx = m_elemIdx;
      m_cell = m_cellBuilder->buildGE();

      // get the states in this cell
      m_cellStates = m_cell->getStates();
      
      // if the states in the cell are parallel updatable, compute the resUpdates (-divFC)
  //    if ((*m_cellStates)[0]->isParUpdatable())
   //   {
	// get the neighbouring faces
        m_faces = m_cell->getNeighborGeos();
      
	// set the cell data
	setCellData();

	// compute the divergence of the discontinuous flux (-divFD+divhFD)
	computeDivDiscontFlx(m_divContFlx);
      
	// update RHS
        updateRHS();
  //    } 
      
      // divide by the Jacobian to transform the residuals back to the physical domain
      //divideByJacobDet();
      
//       for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
//       {
//         (*((*m_cellStates)[iSol]))[1] = m_solEpsilons[iSol];
//       }
//       DataHandle<CFreal> updateCoeff = socket_updateCoeff.getDataHandle();
//       for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
//       {
//         (*((*m_cellStates)[iSol]))[3] = updateCoeff[(*m_cellStates)[iSol]->getLocalID()];
//       }
      
      // print out the residual updates for debugging
      if(m_cell->getID() == 1988) //
      {
	CFLog(VERBOSE, "ID = " << (*m_cellStates)[0]->getLocalID() << "\n");
        CFLog(VERBOSE, "TotalUpdate = \n");
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

void LLAVFluxReconstruction::computeInterfaceFlxCorrection()
{ 
  // Loop over the flux points to calculate FI
  for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
  { 
    const CFreal epsilon = 0.5*(m_epsilonLR[LEFT][iFlxPnt]+m_epsilonLR[RIGHT][iFlxPnt]);
    
    // compute the average sol and grad to use the BR2 scheme
    for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
    {
      *(m_avgGrad[iVar]) = (*(m_cellGradFlxPnt[LEFT][iFlxPnt][iVar]) + *(m_cellGradFlxPnt[RIGHT][iFlxPnt][iVar]))/2.0;
    }
    
    m_flxPntRiemannFlux[iFlxPnt] = 0.0;
    
    for (CFuint iDim = 0; iDim < m_dim; ++iDim)
    {
      for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
      {
        m_flxPntRiemannFlux[iFlxPnt][iVar] += epsilon*((*(m_avgGrad[iVar]))[iDim])*m_unitNormalFlxPnts[iFlxPnt][iDim];
      }
    }
     
    // compute FI in the mapped coord frame
    m_cellFlx[LEFT][iFlxPnt] = (m_flxPntRiemannFlux[iFlxPnt])*m_faceJacobVecSizeFlxPnts[iFlxPnt][LEFT];
    m_cellFlx[RIGHT][iFlxPnt] = (m_flxPntRiemannFlux[iFlxPnt])*m_faceJacobVecSizeFlxPnts[iFlxPnt][RIGHT];
  }
}

//////////////////////////////////////////////////////////////////////////////

void LLAVFluxReconstruction::setFaceData(CFuint faceID)
{
  DiffRHSFluxReconstruction::setFaceData(faceID);
  
  m_faceNodes = m_face->getNodes();
  
  // loop over flx pnts to extrapolate the states to the flux points
  for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
  {   
    //m_epsilonLR[LEFT][iFlxPnt] = m_cellEpsilons[m_cells[LEFT]->getID()];
    //m_epsilonLR[RIGHT][iFlxPnt] = m_cellEpsilons[m_cells[RIGHT]->getID()];
    
    // local flux point indices in the left and right cell
    const CFuint flxPntIdxL = (*m_faceFlxPntConnPerOrient)[m_orient][LEFT][iFlxPnt];
    const CFuint flxPntIdxR = (*m_faceFlxPntConnPerOrient)[m_orient][RIGHT][iFlxPnt]; 
    
    // reset the states in the flx pnts
    m_epsilonLR[LEFT][iFlxPnt] = 0.0;
    m_epsilonLR[RIGHT][iFlxPnt] = 0.0;
    
    for (CFuint iSide = 0; iSide < 2; ++iSide)
    {
      m_cellNodes = m_cells[iSide]->getNodes();
      CFuint flxIdx;
      iSide == LEFT ? flxIdx = flxPntIdxL : flxIdx = flxPntIdxR;

      // loop over the sol pnts to compute the states and grads in the flx pnts
      for (CFuint iNode = 0; iNode < m_faceNodes->size(); ++iNode)
      {
	for (CFuint iNodeCell = 0; iNodeCell < m_nbrCornerNodes; ++iNodeCell)
        {
	  if ((*m_faceNodes)[iNode]->getLocalID() == (*m_cellNodes)[iNodeCell]->getLocalID())
	  {
            //const CFuint nodeIdx = (*m_faceNodes)[iNode]->getLocalID();
	    // get node local index
            const CFuint nodeIdx = (*m_cellNodesConn)(m_cells[iSide]->getID(),iNodeCell);
	    
            m_epsilonLR[iSide][iFlxPnt] += m_nodePolyValsAtFlxPnts[flxIdx][iNodeCell]*m_nodeEpsilons[nodeIdx]/m_nbNodeNeighbors[nodeIdx];
	  }
	}
      }
      //CFLog(INFO, "interEps: " << m_epsilonLR[iSide][iFlxPnt] << "\n");
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void LLAVFluxReconstruction::computeWaveSpeedUpdates(vector< CFreal >& waveSpeedUpd)
{
  // compute the wave speed updates for the neighbouring cells
  cf_assert(waveSpeedUpd.size() == 2);
  CFreal visc = 1.0;
  
  for (CFuint iSide = 0; iSide < 2; ++iSide)
  {
    waveSpeedUpd[iSide] = 0.0;
    for (CFuint iFlx = 0; iFlx < m_cellFlx[iSide].size(); ++iFlx)
    {
      const CFreal jacobXJacobXIntCoef = m_faceJacobVecAbsSizeFlxPnts[iFlx]*
                                         m_faceJacobVecAbsSizeFlxPnts[iFlx]*
                                         (*m_faceIntegrationCoefs)[iFlx]*
                                         m_cflConvDiffRatio;
      //const CFreal rho = (*(m_cellStatesFlxPnt[iSide][iFlx]))[0];
      const CFreal epsilon = 0.5*(m_epsilonLR[LEFT][iFlx]+m_epsilonLR[RIGHT][iFlx]);
      const CFreal viscCoef = computeViscCoef(m_cellStatesFlxPnt[iSide][iFlx]);
      visc = epsilon*viscCoef;
      
      // transform update states to physical data to calculate eigenvalues
      waveSpeedUpd[iSide] += visc*jacobXJacobXIntCoef/m_cellVolume[iSide];
    }
    //if (waveSpeedUpd[iSide] > 10.0) CFLog(INFO, "wvspLLAV: " << waveSpeedUpd[iSide] << "\n");
  }
}

//////////////////////////////////////////////////////////////////////////////

void LLAVFluxReconstruction::computeDivDiscontFlx(vector< RealVector >& residuals)
{
  // reset the extrapolated fluxes
  for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrTotalFlxPnts; ++iFlxPnt)
  {
    m_extrapolatedFluxes[iFlxPnt] = 0.0;
  }

  // Loop over solution points to calculate the discontinuous flux.
  for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolPnts; ++iSolPnt)
  { 
    vector< RealVector >& temp = *(m_cellGrads[0][iSolPnt]);

    // calculate the discontinuous flux projected on x, y, z-directions
    for (CFuint iDim = 0; iDim < m_dim+m_ndimplus; ++iDim)
    {
      m_contFlx[iSolPnt][iDim] = 0.0;

      for (CFuint iDim2 = 0; iDim2 < m_dim; ++iDim2)
      {
        for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
        {
          m_contFlx[iSolPnt][iDim][iVar] += m_solEpsilons[iSolPnt]*(temp[iVar][iDim2])*m_cellFluxProjVects[iDim][iSolPnt][iDim2];
        }
      }
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
    residuals[iSolPnt] = 0.0;

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
          residuals[iSolPnt][iEq] += polyCoef*(m_contFlx[jSolIdx][iDir+m_ndimplus][iEq]);
	}
      }
    }
  }

  // add the contribution of the faces
  const CFuint nbrFaces = m_cell->nbNeighborGeos();

  for (CFuint iFace = 0; iFace < nbrFaces; ++iFace)
  {
    if (!((*m_isFaceOnBoundaryCell)[iFace]))
    {
      for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
      {
        const CFuint currFlxIdx = (*m_faceFlxPntConn)[iFace][iFlxPnt];

        for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolDep; ++iSolPnt)
        {
          const CFuint solIdx = (*m_flxSolDep)[currFlxIdx][iSolPnt];
          const CFreal divh = m_corrFctDiv[solIdx][currFlxIdx];
   
          // Fill in the corrections
          for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
          {
            residuals[solIdx][iVar] += -m_extrapolatedFluxes[currFlxIdx][iVar] * divh; 
          }
        }
      }
    }
    else
    {
      m_faceNodes = (*m_faces)[iFace]->getNodes();
      m_face = (*m_faces)[iFace];
      m_cellNodes = m_cell->getNodes();
      
      // get the datahandle of the update coefficients
      DataHandle<CFreal> updateCoeff = socket_updateCoeff.getDataHandle();
  
      // compute flux point coordinates
      for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
      {
        m_flxPntCoords[iFlx] = m_face->computeCoordFromMappedCoord((*m_flxLocalCoords)[iFlx]);	
      }
          
      // compute face Jacobian vectors
      m_faceJacobVecs = m_face->computeFaceJacobDetVectorAtMappedCoords(*m_flxLocalCoords);
  
      // get face Jacobian vector sizes in the flux points
      DataHandle< vector< CFreal > > faceJacobVecSizeFaceFlxPnts = socket_faceJacobVecSizeFaceFlxPnts.getDataHandle();
  
      // Loop over flux points to compute the unit normals
      for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
      {
        // get face Jacobian vector size
        CFreal faceJacobVecAbsSizeFlxPnts = faceJacobVecSizeFaceFlxPnts[m_face->getID()][iFlxPnt];
	
	// set face Jacobian vector size with sign depending on mapped coordinate direction
        m_faceJacobVecSizeFlxPnts2[iFlxPnt] = faceJacobVecAbsSizeFlxPnts*((*m_faceLocalDir)[iFace]);
 
	// set unit normal vector
        m_unitNormalFlxPnts2[iFlxPnt] = (m_faceJacobVecs[iFlxPnt]/faceJacobVecAbsSizeFlxPnts);
      }
	
      for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
      {
        const CFuint currFlxIdx = (*m_faceFlxPntConn)[iFace][iFlxPnt];
    
        // reset the grads in the flx pnts
        for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
        {
          *(m_cellGradFlxPnt[0][iFlxPnt][iVar]) = 0.0;
        }
        
        *(m_cellStatesFlxPnt[0][iFlxPnt]) = 0.0;

        // loop over the sol pnts to compute the states and grads in the flx pnts
        for (CFuint iSol = 0; iSol < m_nbrSolDep; ++iSol)
        {
          const CFuint solIdx = (*m_flxSolDep)[currFlxIdx][iSol];
          const CFreal coeff = (*m_solPolyValsAtFlxPnts)[currFlxIdx][solIdx];

	  *(m_cellStatesFlxPnt[0][iFlxPnt]) += coeff*(*((*(m_cellStates))[solIdx]));
	  
          for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
          {
            *(m_cellGradFlxPnt[0][iFlxPnt][iVar]) += coeff*((*(m_cellGrads[0][solIdx]))[iVar]);
          }
        }
      }
	
      // compute ghost gradients
      if ((getMethodData().getUpdateVarStr() == "Cons" || getMethodData().getUpdateVarStr() == "RhoivtTv") && getMethodData().hasDiffTerm())
      {
	for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
        {
	  for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
          {
	    *(m_flxPntGhostGrads[iFlxPnt][iVar]) = *(m_cellGradFlxPnt[0][iFlxPnt][iVar]);
	  }
	}
      }
      else
      {
	(*m_bcStateComputers)[(*m_faceBCIdxCell)[iFace]]->computeGhostGradients(m_cellGradFlxPnt[0],m_flxPntGhostGrads,m_unitNormalFlxPnts2,m_flxPntCoords);
      }
	
      for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
      {
	const CFuint currFlxIdx = (*m_faceFlxPntConn)[iFace][iFlxPnt];
	
	CFreal epsilon = 0.0;
	
	// loop over the sol pnts to compute the states and grads in the flx pnts
        for (CFuint iNode = 0; iNode < m_faceNodes->size(); ++iNode)
        {
	  for (CFuint iNodeCell = 0; iNodeCell < m_nbrCornerNodes; ++iNodeCell)
          {
	    if ((*m_faceNodes)[iNode]->getLocalID() == (*m_cellNodes)[iNodeCell]->getLocalID())
	    {
              //const CFuint nodeIdx = (*m_faceNodes)[iNode]->getLocalID();
	      // get node local index
              const CFuint nodeIdx = (*m_cellNodesConn)(m_cell->getID(),iNodeCell);
	  
              epsilon += m_nodePolyValsAtFlxPnts[currFlxIdx][iNodeCell]*m_nodeEpsilons[nodeIdx]/m_nbNodeNeighbors[nodeIdx];
	    }
	  }
        }
	
	if (!m_jacob && m_addUpdCoeff)
	{
	  // adding updateCoeff
	  CFreal visc = 1.0;
  
          m_waveSpeedUpd[0] = 0.0;

          const CFreal jacobXJacobXIntCoef = m_faceJacobVecSizeFlxPnts2[iFlxPnt]*
                                             m_faceJacobVecSizeFlxPnts2[iFlxPnt]*
                                             (*m_faceIntegrationCoefs)[iFlxPnt]*
                                             m_cflConvDiffRatio;
          //const CFreal rho = (*(m_cellStatesFlxPnt[0][iFlxPnt]))[0];
          const CFreal viscCoef = computeViscCoef(m_cellStatesFlxPnt[0][iFlxPnt]);
          visc = epsilon*viscCoef;
      
          // transform update states to physical data to calculate eigenvalues
          m_waveSpeedUpd[0] += visc*jacobXJacobXIntCoef/m_cell->computeVolume();

          // loop over the sol pnts of both sides to update the wave speeds
          for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
          {
            const CFuint solID = (*m_cellStates)[iSol]->getLocalID();
            updateCoeff[solID] += m_waveSpeedUpd[0];
          }
          //if (m_waveSpeedUpd[0] > 10.0) CFLog(INFO, "wvspLLAVBnd: " << m_waveSpeedUpd[0] << "\n");
	}
	
        // compute the average sol and grad to use the BR2 scheme
        for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
        {
	  if (m_cell->getID() == 1092) CFLog(VERBOSE, "var: " << iVar << ", grad: " << *(m_cellGradFlxPnt[0][iFlxPnt][iVar]) << ", ghost: " << *(m_flxPntGhostGrads[iFlxPnt][iVar]) << "\n");
          *(m_avgGrad[iVar]) = (*(m_cellGradFlxPnt[0][iFlxPnt][iVar]) + *(m_flxPntGhostGrads[iFlxPnt][iVar]))/2.0;
        }
              
        m_flxPntRiemannFlux[iFlxPnt] = 0.0;
	      
        for (CFuint iDim = 0; iDim < m_dim; ++iDim)
        {
          for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
          {
            m_flxPntRiemannFlux[iFlxPnt][iVar] += epsilon*((*(m_avgGrad[iVar]))[iDim])*m_unitNormalFlxPnts2[iFlxPnt][iDim];
	    if (m_cell->getID() == 1092) CFLog(VERBOSE, "avgrad: " << (*(m_avgGrad[iVar]))[iDim] << "\n");
          }
        }

        // compute FI in the mapped coord frame
        m_cellFlx[0][iFlxPnt] = (m_flxPntRiemannFlux[iFlxPnt])*m_faceJacobVecSizeFlxPnts2[iFlxPnt]; 
	if (m_cell->getID() == 1092) CFLog(VERBOSE, "riemannunit: " << m_flxPntRiemannFlux[iFlxPnt] << "jacob: " << m_faceJacobVecSizeFlxPnts2[iFlxPnt] << "\n");
	
	for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolDep; ++iSolPnt)
        {  
          const CFuint solIdx = (*m_flxSolDep)[currFlxIdx][iSolPnt];
          const CFreal divh = m_corrFctDiv[solIdx][currFlxIdx];
   
          // Fill in the corrections
          for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
          {
            residuals[solIdx][iVar] += (m_cellFlx[0][iFlxPnt][iVar] - m_extrapolatedFluxes[currFlxIdx][iVar]) * divh; 
          }
        }
      }
    }
}
}

//////////////////////////////////////////////////////////////////////////////

CFreal LLAVFluxReconstruction::computeViscCoef(RealVector* state)
{
  return 1.;
}

//////////////////////////////////////////////////////////////////////////////

void LLAVFluxReconstruction::setCellData()
{
  DiffRHSFluxReconstruction::setCellData();
  
  m_cellNodes = m_cell->getNodes();
  
  // get datahandle
  DataHandle< CFreal > artVisc = socket_artVisc.getDataHandle();
  
  // loop over flx pnts to extrapolate the states to the flux points
  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {   
    //m_solEpsilons[iSol] = m_cellEpsilons[m_cell->getID()];
    
    // reset the states in the flx pnts
    m_solEpsilons[iSol] = 0.0;

    // loop over the sol pnts to compute the states and grads in the flx pnts
    for (CFuint iNode = 0; iNode < m_nbrCornerNodes; ++iNode)
    {
      // get node local index
      //const CFuint nodeIdx = (*m_cellNodesConn)(m_elemIdx,iNode);
      
      const CFuint nodeIdx = (*m_cellNodes)[iNode]->getLocalID();
      
      m_solEpsilons[iSol] += m_nodePolyValsAtSolPnts[iSol][iNode]*m_nodeEpsilons[nodeIdx]/m_nbNodeNeighbors[nodeIdx];
      //CFLog(VERBOSE, "solEps: " << m_solEpsilons[iSol] << ", nodeEps: " << m_nodeEpsilons[nodeIdx] << ", nghb: " << m_nbNodeNeighbors[nodeIdx] << "\n");
    }
    
    //CFLog(INFO, "eps: " << m_solEpsilons[iSol] << "\n");
    
    artVisc[(((*m_cellStates)[iSol]))->getLocalID()] = m_solEpsilons[iSol];
    
//       CFuint ID = m_cell->getID();
//   bool cond = ID == 51 || ID == 233 || ID == 344 || ID == 345 || ID == 389 || ID == 3431 || ID == 3432 || ID == 3544 || ID == 3545;
//   if (cond)
//   {
//     CFLog(INFO, "nodal eps = " << m_solEpsilons[iSol] << "\n");
//   }
  }
}

//////////////////////////////////////////////////////////////////////////////

void LLAVFluxReconstruction::computeProjStates(std::vector< RealVector >& projStates)
{
  cf_assert(m_nbrSolPnts == projStates.size());
  
  if (m_order != 1)
  {
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
      {
        m_tempSolPntVec[iSol] = (*((*m_cellStates)[iSol]))[iEq];
      }

      m_tempSolPntVec2 = m_transformationMatrix*m_tempSolPntVec;

      for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
      {
        projStates[iSol][iEq] = m_tempSolPntVec2[iSol];
      }
    }
  }
  else
  {
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      CFreal stateSum = 0.0;
      
      for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
      {
        stateSum += (*((*m_cellStates)[iSol]))[iEq];
      }

      stateSum /= m_nbrSolPnts;

      for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
      {
        projStates[iSol][iEq] = stateSum;
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void LLAVFluxReconstruction::computeEpsilon()
{
  computeEpsilon0();
  
  computeSmoothness();
  
  if (m_s < m_s0 - m_kappa)
  {
    m_epsilon = 0.0;
  }
  else if (m_s > m_s0 + m_kappa)
  {
    m_epsilon = m_epsilon0;
  }
  else
  {
    m_epsilon = m_epsilon0*0.5*(1.0 + sin(0.5*MathTools::MathConsts::CFrealPi()*(m_s-m_s0)/m_kappa));
  }
  
  //CFLog(INFO, "eps0: " << m_epsilon0 << ", S: " << m_s << ", eps: " << m_epsilon << "\n");
  
//   for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
//       {
//         (*((*m_cellStates)[iSol]))[2] = m_epsilon0;
//       }
  if (m_epsilon < 0.0 || m_epsilon != m_epsilon) 
  {
      CFLog(INFO, "eps: " << m_epsilon << ", eps0: " << m_epsilon0 << ", s: " << m_s << "\n");
      m_epsilon = 0.0;
  }
  //cf_assert(m_epsilon > -1.0e-8);
  //if (m_epsilon >= m_epsilon0*0.0001) CFLog(INFO, "eps0 = " << m_epsilon0 << ", DS: " << m_s-m_s0 << "\n");
//   CFuint ID = m_cell->getID();
//   bool cond = ID == 51 || ID == 233 || ID == 344 || ID == 345 || ID == 389 || ID == 3431 || ID == 3432 || ID == 3544 || ID == 3545;
//   if (cond)
//   {
//     CFLog(INFO, "eps = " << m_epsilon << ", DS: " << m_s-m_s0 << ", eps0: " << m_epsilon0 << "\n");
//   }
}

//////////////////////////////////////////////////////////////////////////////

void LLAVFluxReconstruction::computeEpsilon0()
{ 
  // get the datahandle of the update coefficients
  DataHandle<CFreal> updateCoeff = socket_updateCoeff.getDataHandle();
  
  const CFreal wavespeed = updateCoeff[(*m_cellStates)[0]->getLocalID()];
  
  //if (wavespeed < 1e-8) CFLog(INFO, "wvsp: " << wavespeed << "\n");
  //if (wavespeed > 1e5) CFLog(INFO, "state HW: " << *((*m_cellStates)[0]) << "\n");
  
  m_maxLambda = max(m_maxLambda,wavespeed);
  
  //const CFreal deltaKsi = 2.0/(m_order+2.0);
  
  const CFreal peclet = computePeclet();
  
  m_epsilon0 = max(wavespeed*(2.0/peclet - m_subcellRes/peclet),0.0);
  
  //CFLog(INFO, "eps0: " << m_epsilon0 << "\n");
  
  //if (m_epsilon0 < -1.0e-8) CFLog(INFO, "wvsp: " << wavespeed << "\n");
  if (m_epsilon0 != m_epsilon0) CFLog(INFO, "wvsp: " << wavespeed << ", state: " << *((*m_cellStates)[0]) << ", coord: " << (*m_cellStates)[0]->getCoordinates() << "\n");
  
  if (m_addPosPrev) addPositivityPreservation();

  //CFLog(INFO, "lambda: " << wavespeed << "\n");
  
  //CFreal vol = m_cell->computeVolume();
  //m_epsilon0 = m_peclet*sqrt(vol);

}

//////////////////////////////////////////////////////////////////////////////

CFreal LLAVFluxReconstruction::computePeclet()
{
  return m_peclet;
}

//////////////////////////////////////////////////////////////////////////////

void LLAVFluxReconstruction::addPositivityPreservation()
{ 
  //CFreal minState = MathTools::MathConsts::CFrealMax();
  
  DataHandle< CFreal > posPrev = socket_posPrev.getDataHandle();
  
  const CFreal posPrevValue = posPrev[m_cell->getID()];
  
//   for (CFuint iFlxPnt = 0; iFlxPnt < m_flxPntsLocalCoords->size(); ++iFlxPnt)
//   {
//     RealVector extrapolatedState;
//     extrapolatedState.resize(m_nbrEqs);
//     extrapolatedState = 0.0;
// 
//     // loop over the sol pnts to compute the states and grads in the flx pnts
//     for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
//     {
//       extrapolatedState += (*m_solPolyValsAtFlxPnts)[iFlxPnt][iSol]*(*((*(m_cellStates))[iSol]));
//     }
//     
//     CFreal rho = extrapolatedState[0];
//     CFreal rhoU = extrapolatedState[1];
//     CFreal rhoV = extrapolatedState[2];
//     CFreal rhoE = extrapolatedState[3];
//     CFreal press = 0.4*(rhoE - 0.5*(rhoU*rhoU+rhoV*rhoV)/rho);
//       
//     minState = min(minState, press);
//   }
//   
//   cf_assert(minState > 0.0);
  
  if (posPrevValue < m_minValue)  //minState < m_minValue
  {
    //m_epsilon0 *= m_minValue/minState;
    //m_epsilon0 *= (m_viscFactor - 1.0)/(m_minValue*m_minValue)*minState*minState - 2.0*(m_viscFactor - 1.0)/m_minValue*minState + m_viscFactor;
    const CFreal factor = min(10*m_viscFactor, (m_viscFactor - 1.0)/(m_minValue*m_minValue)*posPrevValue*posPrevValue - 2.0*(m_viscFactor - 1.0)/m_minValue*posPrevValue + m_viscFactor);
    m_epsilon0 *= factor;
    m_nbPosPrev++;
  }
}

//////////////////////////////////////////////////////////////////////////////

void LLAVFluxReconstruction::computeSmoothness()
{ 
  CFreal sNum = 0.0;
  
  CFreal sDenom = 0.0;
  
  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {
    CFreal stateP = (*((*m_cellStates)[iSol]))[m_monitoredVar];
    CFreal diffStatesPPMinOne = stateP - m_statesPMinOne[iSol][m_monitoredVar];
    sNum += diffStatesPPMinOne*diffStatesPPMinOne;
    sDenom += stateP*stateP;
  }
  if (sNum <= MathTools::MathConsts::CFrealEps() || sDenom <= MathTools::MathConsts::CFrealEps())
  {
    m_s = -100.0;
  }
  else
  {
    m_s = log10(sNum/sDenom);
  }
  
  // get datahandle
  DataHandle< CFreal > smoothness = socket_smoothness.getDataHandle();
  
  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {
    smoothness[(((*m_cellStates)[iSol]))->getLocalID()] = m_s;
  }
  
  if (m_s > m_Smax)
  {
      m_Smax = m_s;
  }
//   CFuint ID = m_cell->getID();
//   bool cond = ID == 51 || ID == 233 || ID == 344 || ID == 345 || ID == 389 || ID == 3431 || ID == 3432 || ID == 3544 || ID == 3545;
//   if (cond)
//   {
//     CFLog(INFO, "S = " << m_s << ", num = " << sNum << ", denom = " << sDenom << "\n");
//   }
}

//////////////////////////////////////////////////////////////////////////////

void LLAVFluxReconstruction::storeEpsilon()
{ 
  for (CFuint iNode = 0; iNode < m_nbrCornerNodes; ++iNode)
  {
    // get node ID
    const CFuint nodeID = (*m_cellNodes)[iNode]->getLocalID();

    if (!m_useMax) 
    {
      m_nodeEpsilons[nodeID] += m_epsilon;
      m_cellEpsilons[m_cell->getID()] = m_epsilon;
      m_totalEps += m_epsilon;
    }
    else
    {
      const CFreal maxEps = max(m_epsilon, m_cellEpsilons[m_cell->getID()]);
      m_nodeEpsilons[nodeID] += maxEps;
      m_cellEpsilons[m_cell->getID()] = maxEps;
      m_totalEps += maxEps;
    }
    
    if (m_flagComputeNbNghb)
    {
      m_nbNodeNeighbors[nodeID] += 1.0;
    }
  }
  
  //CFLog(INFO, "solEps: " << m_cellEpsilons[m_cell->getID()] << "\n");
  
  if (m_epsilon > 0.5*m_epsilon0) CFLog(VERBOSE, "cellID eps: " << m_cell->getID() << "\n");
  CFLog(VERBOSE, "eps0 = " << m_epsilon0 << ", eps = " << m_epsilon << ", S = " << m_s << ", S0 = " << m_s0 << "\n");
}

//////////////////////////////////////////////////////////////////////////////

void LLAVFluxReconstruction::setup()
{
  CFAUTOTRACE;
  
  // setup parent class
  DiffRHSFluxReconstruction::setup();

  // get the update varset
  m_updateVarSet = getMethodData().getUpdateVar();
  
  // get CellToFaceGeBuilder
  m_cellBuilder = getMethodData().getCellBuilder();
  m_isFaceOnBoundaryCell = m_cellBuilder->getGeoBuilder()->getIsFaceOnBoundary();
  m_nghbrCellSideCell    = m_cellBuilder->getGeoBuilder()->getNeighbrCellSide ();
  m_currCellSideCell     = m_cellBuilder->getGeoBuilder()->getCurrentCellSide ();
  m_faceOrientsCell      = m_cellBuilder->getGeoBuilder()->getFaceOrient      ();
  m_faceBCIdxCell        = m_cellBuilder->getGeoBuilder()->getFaceBCIdx       ();
  
  // get BC state computers
  m_bcStateComputers = getMethodData().getBCStateComputers();
  
  // get cell-node connectivity
  m_cellNodesConn = MeshDataStack::getActive()->getConnectivity("cellNodes_InnerCells");
  
  // get the local spectral FD data
  vector< FluxReconstructionElementData* >& frLocalData = getMethodData().getFRLocalData();
  cf_assert(frLocalData.size() > 0);
  // for now, there should be only one type of element
  cf_assert(frLocalData.size() == 1);
  
  const CFPolyOrder::Type order = frLocalData[0]->getPolyOrder();
  
  m_subcellRes = frLocalData[0]->getSubcellResolution();
  
  // get the coefs for extrapolation of the node artificial viscosities to the flx pnts
  m_nodePolyValsAtFlxPnts = frLocalData[0]->getNodePolyValsAtPnt(*(frLocalData[0]->getFlxPntsLocalCoords()));
  
  // get the coefs for extrapolation of the node artificial viscosities to the sol pnts
  m_nodePolyValsAtSolPnts = frLocalData[0]->getNodePolyValsAtPnt(*(frLocalData[0]->getSolPntsLocalCoords()));
  
  m_order = static_cast<CFuint>(order);

  // number of cell corner nodes
  /// @note in the future, hanging nodes should be taken into account here
  m_nbrCornerNodes = frLocalData[0]->getNbrCornerNodes();
  
  // get the number of nodes in the mesh
  const CFuint nbrNodes = MeshDataStack::getActive()->getNbNodes();
  
  // get the elementTypeData
  SafePtr< vector<ElementTypeData> > elemType = MeshDataStack::getActive()->getElementTypeData();
  
  // get the number of cells in the mesh
  const CFuint nbrCells = (*elemType)[0].getEndIdx();
  
  // get datahandle
  DataHandle< CFreal > artVisc = socket_artVisc.getDataHandle();
  DataHandle< CFreal > monPhysVar = socket_monPhysVar.getDataHandle();
  DataHandle< CFreal > smoothness = socket_smoothness.getDataHandle();
  
  const CFuint nbStates = nbrCells*m_nbrSolPnts;

  // resize socket
  artVisc.resize(nbStates);
  monPhysVar.resize(nbStates);
  smoothness.resize(nbStates);
  
  m_nodeEpsilons.resize(nbrNodes);
  m_nbNodeNeighbors.resize(nbrNodes);
  m_cellEpsilons.resize(nbrCells);
  m_solEpsilons.resize(m_nbrSolPnts);
  m_epsilonLR.resize(2);
  m_epsilonLR[LEFT].resize(m_nbrFaceFlxPnts);
  m_epsilonLR[RIGHT].resize(m_nbrFaceFlxPnts);
  m_unitNormalFlxPnts2.resize(m_nbrFaceFlxPnts);
  m_faceJacobVecSizeFlxPnts2.resize(m_nbrFaceFlxPnts);
  m_tempSolPntVec.resize(m_nbrSolPnts);
  m_tempSolPntVec2.resize(m_nbrSolPnts);
  
  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {
    RealVector temp(m_nbrEqs);
    temp = 0.0;
    m_statesPMinOne.push_back(temp);
  }
  
  SafePtr<RealMatrix> vdm = frLocalData[0]->getVandermondeMatrix();
  
  SafePtr<RealMatrix> vdmInv = frLocalData[0]->getVandermondeMatrixInv();
  
  RealMatrix temp(m_nbrSolPnts,m_nbrSolPnts);
  temp = 0.0;
  if (m_dim == 2)
  {
    if (m_ndimplus==3){  //if Triag
      for (CFuint idx = 0; idx < (m_order)*(m_order+1)/2; ++idx)
      {
        temp(idx,idx) = 1.0;
      }
    }
    else{
      for (CFuint idx = 0; idx < (m_order)*(m_order); ++idx)
      {
        temp(idx,idx) = 1.0;
      }
    }  
  }
  else if (m_dim == 3)
  {
    if (m_ndimplus==4){  //if Tetra
      for (CFuint idx = 0; idx < (m_order)*(m_order+1)*(m_order+2)/6; ++idx)
      {
        temp(idx,idx) = 1.0;
      }
    }
    else{
       for (CFuint idx = 0; idx < (m_order)*(m_order)*(m_order); ++idx)
      {
        temp(idx,idx) = 1.0;
      }
    }  
  }
  
  m_transformationMatrix.resize(m_nbrSolPnts,m_nbrSolPnts);
  
  m_transformationMatrix = (*vdm)*temp*(*vdmInv);
  
  //m_s0 = -m_s0*log10(static_cast<CFreal>(m_order)) + m_s0;
  
  m_Smax = m_s0 + m_kappa;
  
  m_SmaxGlobal = m_Smax;
  
  m_nbNodeNeighbors = 0.0;
  
  m_flagComputeNbNghb = true;
  
  m_flxPntGhostGrads.resize(m_nbrFaceFlxPnts);
  
  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    m_unitNormalFlxPnts2[iFlx].resize(m_dim);

    for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
    {
      m_flxPntGhostGrads[iFlx].push_back(new RealVector(m_dim));
    }
  }
  
  cf_assert(m_monitoredVar < m_nbrEqs);
}

//////////////////////////////////////////////////////////////////////////////

void LLAVFluxReconstruction::unsetup()
{
  CFAUTOTRACE;
  
  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    for (CFuint iGrad = 0; iGrad < m_nbrEqs; ++iGrad)
    {
      deletePtr(m_flxPntGhostGrads[iFlx][iGrad]); 
    }
    m_flxPntGhostGrads[iFlx].clear();
  }
  m_flxPntGhostGrads.clear();
  
  // unsetup parent class
  DiffRHSFluxReconstruction::unsetup();
}


//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD

