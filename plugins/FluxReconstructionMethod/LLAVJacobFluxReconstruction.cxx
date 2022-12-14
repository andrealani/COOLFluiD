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

#include "FluxReconstructionMethod/LLAVJacobFluxReconstruction.hh"
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

MethodCommandProvider< LLAVJacobFluxReconstruction,
		       FluxReconstructionSolverData,
		       FluxReconstructionModule >
LLAVJacobFluxReconstructionFluxReconstructionProvider("LLAVJacob");

//////////////////////////////////////////////////////////////////////////////
  
LLAVJacobFluxReconstruction::LLAVJacobFluxReconstruction(const std::string& name) :
  DiffRHSJacobFluxReconstruction(name),
  m_updateVarSet(CFNULL),
  m_order(),
  m_transformationMatrix(),
  m_statesPMinOne(),
  m_epsilon(),
  m_solEpsilons(),
  m_epsilonLR(),
  m_epsilon0(),
  m_s0(),
  m_s(),
  m_kappa(),
  m_showrate(),
  m_peclet(),
  m_nodeEpsilons(),
  m_nbNodeNeighbors(),
  m_cellEpsilons(),
  m_cellNodes(),
  m_nbrCornerNodes(),
  m_faceNodes(),
  m_flagComputeNbNghb(),
  m_nodePolyValsAtFlxPnts(),
  m_nodePolyValsAtSolPnts(),
  m_cellNodesConn(CFNULL),
  m_elemIdx(),
  m_facesCell(),
  m_jacob(),
  m_freezeLimiterRes(),
  m_freezeLimiterIter(),
  m_useMax(),
  m_totalEps(),
  m_totalEpsGlobal(),
  m_Smax(),
  m_SmaxGlobal(),
  m_subcellRes(),
  m_unitNormalFlxPnts2(),
  socket_artVisc("artVisc"),
  socket_monPhysVar("monPhysVar"),
  socket_smoothness("smoothness"),
  m_epsBackUp(),
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
    
    m_dampingCoeff = 1.0;
    setParameter( "DampingCoeff", &m_dampingCoeff);
    
    m_freezeLimiterRes = -20.0;
    setParameter( "FreezeLimiterRes", &m_freezeLimiterRes);
  
    m_freezeLimiterIter = MathTools::MathConsts::CFuintMax();
    setParameter( "FreezeLimiterIter", &m_freezeLimiterIter);
    
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

void LLAVJacobFluxReconstruction::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal,Config::DynamicOption<> >("Kappa","Kappa factor of artificial viscosity.");
  
  options.addConfigOption< CFuint,Config::DynamicOption<> >("ShowRate","Showrate of LLAV information.");
  
  options.addConfigOption< CFreal,Config::DynamicOption<> >("Peclet","Peclet number to be used for artificial viscosity.");
  
  options.addConfigOption< CFreal,Config::DynamicOption<> >("S0","Reference smoothness factor, will be multiplied by -log(P).");
  
  options.addConfigOption< CFreal,Config::DynamicOption<> >("DampingCoeff","Damping coefficient for reculculation of eps (0<coeff<1).");
  
  options.addConfigOption< CFreal,Config::DynamicOption<> >("FreezeLimiterRes","Residual after which to freeze the residual.");
  
  options.addConfigOption< CFuint,Config::DynamicOption<> >("FreezeLimiterIter","Iteration after which to freeze the residual.");
  
  options.addConfigOption< CFreal >("ViscFactor","Maximum factor applied to viscosity for positivity preservation.");
  
  options.addConfigOption< CFuint,Config::DynamicOption<> >("MonitoredVar","Index of the monitored var for positivity preservation.");
  
  options.addConfigOption< bool >("AddUpdateCoeff","Boolean telling whether the update coefficient based on the artificial flux is added.");
  
  options.addConfigOption< CFuint,Config::DynamicOption<> >("MonitoredPhysVar","Index of the monitored physical var for positivity preservation, if not specified MonitoredVar is used instead.");
}

//////////////////////////////////////////////////////////////////////////////

void LLAVJacobFluxReconstruction::configure ( Config::ConfigArgs& args )
{
  FluxReconstructionSolverCom::configure(args);
}  

//////////////////////////////////////////////////////////////////////////////

std::vector< Common::SafePtr< BaseDataSocketSource > >
  LLAVJacobFluxReconstruction::providesSockets()
{
  std::vector< Common::SafePtr< BaseDataSocketSource > > result;
  result.push_back(&socket_artVisc);
  result.push_back(&socket_monPhysVar);
  result.push_back(&socket_smoothness);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

void LLAVJacobFluxReconstruction::execute()
{
  CFAUTOTRACE;
  
  CFLog(VERBOSE, "LLAVJacobFluxReconstruction::execute()\n");
  
  // get the elementTypeData
  SafePtr< vector<ElementTypeData> > elemType = MeshDataStack::getActive()->getElementTypeData();

  // get InnerCells TopologicalRegionSet
  SafePtr<TopologicalRegionSet> cells = MeshDataStack::getActive()->getTrs("InnerCells");

  // get the geodata of the geometric entity builder and set the TRS
  CellToFaceGEBuilder::GeoData& geoDataCell = m_cellBuilder->getDataGE();
  geoDataCell.trs = cells;
  
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
  
  m_nodeEpsilons = 0.0;
  
  const CFreal residual = SubSystemStatusStack::getActive()->getResidual();
  
  const CFuint iter = SubSystemStatusStack::getActive()->getNbIter();
  
  m_useMax = residual < m_freezeLimiterRes || iter > m_freezeLimiterIter;
  m_Smax = -100.0;
  m_totalEps = 0.0;
  
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
      
      // get the nodes in this cell
      m_cellNodes  = m_cell->getNodes();
      
//       // if the states in the cell are parallel updatable, compute the resUpdates (-divFC)
//       if ((*m_cellStates)[0]->isParUpdatable())
//       {
	// compute the states projected on order P-1
	computeProjStates(m_statesPMinOne);
	
	// compute the artificial viscosity
	computeEpsilon();
	
	// store epsilon
	storeEpsilon();
//       } 
      
      //release the GeometricEntity
      m_cellBuilder->releaseGE();
    }
  }
  
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
//      if ((*m_states[LEFT ])[0]->isParUpdatable() || (*m_states[RIGHT])[0]->isParUpdatable())
//      {
	// build the neighbouring cells
        const CFuint cellIDL = m_face->getNeighborGeo(LEFT)->getID();
        geoDataCBL.idx = cellIDL;
        m_cells[LEFT] = m_cellBuilders[LEFT ]->buildGE();
        const CFuint cellIDR = m_face->getNeighborGeo(RIGHT)->getID();
        geoDataCBR.idx = cellIDR;
        m_cells[RIGHT] = m_cellBuilders[RIGHT]->buildGE();
	m_jacob = false;

	// set the face data
	setFaceData(m_face->getID());//faceID

	// compute the left and right states and gradients in the flx pnts
	computeFlxPntStatesAndGrads();

	// compute FI
	computeInterfaceFlxCorrection();

	if (m_addUpdCoeff)
	{
	  // compute the wave speed updates
          computeWaveSpeedUpdates(m_waveSpeedUpd);

          // update the wave speed
          updateWaveSpeed();
	}

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
	
//      // if one of the neighbouring cells is parallel updatable, compute the correction flux
//      if ((*m_states[LEFT ])[0]->isParUpdatable() || (*m_states[RIGHT])[0]->isParUpdatable())
//      {
        
        // get all the faces neighbournig the cells
        m_faces[LEFT ] = m_cells[LEFT ]->getNeighborGeos();
        m_faces[RIGHT] = m_cells[RIGHT]->getNeighborGeos();
        
        // set the local indexes of the other faces than the current faces
        setOtherFacesLocalIdxs();

        setCellData(LEFT);
        setCellData(RIGHT);

	// compute needed cell contributions
        if (!m_cellFlags[cellIDL] && (*m_states[LEFT ])[0]->isParUpdatable())
        {
          m_pertSide = LEFT;

          computeUnpertCellDiffResiduals(LEFT);
          
          m_unpertAllCellDiffRes[cellIDL] = m_unpertCellDiffRes[LEFT];

	  // update RHS
	  updateRHSUnpertCell(LEFT);
        }
        if (!m_cellFlags[cellIDR] && (*m_states[RIGHT])[0]->isParUpdatable())
        {
          m_pertSide = RIGHT;

          computeUnpertCellDiffResiduals(RIGHT);

          m_unpertAllCellDiffRes[cellIDR] = m_unpertCellDiffRes[RIGHT];

	  // update RHS
	  updateRHSUnpertCell(RIGHT);
        }

        m_jacob = true;

	// make a back up of the grads and put the perturbed and unperturbed corrections in the correct format
        for (CFuint iState = 0; iState < m_nbrSolPnts; ++iState)
        {
          for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
          {
            m_cellGradsBackUp[LEFT][iState][iVar] = (*m_cellGrads[LEFT][iState])[iVar];
            m_cellGradsBackUp[RIGHT][iState][iVar] = (*m_cellGrads[RIGHT][iState])[iVar];
            
            m_resUpdates[LEFT][m_nbrEqs*iState+iVar] = m_divContFlxL[iState][iVar];
            m_resUpdates[RIGHT][m_nbrEqs*iState+iVar] = m_divContFlxR[iState][iVar];
          }
        }

	for (CFuint iSide = 0; iSide < 2; ++iSide)
        {
          // compute solution points Jacobian determinants
          m_solJacobDet[iSide] = m_cells[iSide]->computeGeometricShapeFunctionJacobianDeterminant(*m_solPntsLocalCoords);
	}

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

        // release the cells
        m_cellBuilders[LEFT ]->releaseGE();
        m_cellBuilders[RIGHT]->releaseGE();
        
        m_cellFlags[cellIDL] = true;
        m_cellFlags[cellIDR] = true;

//      }

      // release the GeometricEntity
      m_faceBuilder->releaseGE();

    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void LLAVJacobFluxReconstruction::computePerturbedGradientsAnalytical(const CFuint side)
{ 

  // Add the discontinuous gradient
  *m_cellStates = *(m_states[side]);

  const CFreal eps = m_numJacob->getEps();
  
  // Loop over solution pnts to calculate the grad updates
  for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolSolDep; ++iSolPnt)
  {
    const CFuint iSolIdx = (*m_solSolDep)[m_pertSol][iSolPnt];
    
    m_affectedSolPnts[side][iSolIdx] = true;
    
    // inverse Jacobian determinant
    const CFreal invJacobDet = 1.0/m_solJacobDet[side][iSolIdx];
    
    // Loop over gradient directions
    for (CFuint iDir = 0; iDir < m_dim; ++iDir)
    {
      m_projectedCorrL = eps * m_neighbCellFluxProjVects[m_pertSide][iDir+m_ndimplus][m_pertSol];
	  
      // compute the grad updates
      (*m_cellGrads[side][iSolIdx])[m_pertVar] += (*m_solPolyDerivAtSolPnts)[iSolIdx][iDir][m_pertSol]*m_projectedCorrL*invJacobDet;
    }
  }
  
  // get face Jacobian vector sizes in the flux points
  DataHandle< vector< CFreal > > faceJacobVecSizeFaceFlxPnts = socket_faceJacobVecSizeFaceFlxPnts.getDataHandle();
  
  // Perturbed flx pnt idx and cell wide idx
  CFuint pertFlxPnt;
  CFuint pertFlxPntIdx;
  
  // Add the contribution of the correction to the gradients for each face
  // compute other face contributions to the gradients
  const CFuint nbrOtherFaces = m_otherFaceLocalIdxs[side].size();
  
  for (CFuint iFace = 0; iFace < nbrOtherFaces; ++iFace)
  {
    // get local face index
    const CFuint faceIdx = m_otherFaceLocalIdxs[side][iFace];

    //Don't add anything for bnd face, corresponding to enforcing zero AV on bnd
//    if ((*m_isFaceOnBoundary[side])[faceIdx])
//    {
//
//      // compute face Jacobian vectors
//      m_faceJacobVecs = (*m_faces[side])[faceIdx]->computeFaceJacobDetVectorAtMappedCoords(*m_flxLocalCoords);
//        
//      // Loop over flux points to set the normal vectors
//      for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
//      {
//        const CFuint currFlxIdx = (*m_faceFlxPntConn)[faceIdx][iFlxPnt];
//      
//        for (CFuint jFlxPnt = 0; jFlxPnt < m_nbrFlxDep; ++jFlxPnt)
//        {
//          if ((*m_solFlxDep)[m_pertSol][jFlxPnt] == currFlxIdx)
//          {
//            pertFlxPnt = iFlxPnt;
//            pertFlxPntIdx = currFlxIdx;
//            
//            break;
//          }
//        }
//        
//        // get face Jacobian vector size
//        m_faceJacobVecAbsSizeFlxPnts[iFlxPnt] = faceJacobVecSizeFaceFlxPnts[(*m_faces[side])[faceIdx]->getID()][iFlxPnt];
//
//        // set unit normal vector
//        m_unitNormalFlxPnts[iFlxPnt] = m_faceJacobVecs[iFlxPnt]/m_faceJacobVecAbsSizeFlxPnts[iFlxPnt];
//        
//        m_flxPntCoords[iFlxPnt] = (*m_faces[side])[faceIdx]->computeCoordFromMappedCoord((*m_flxLocalCoords)[iFlxPnt]);
//        
//        *(m_cellStatesFlxPnt[0][iFlxPnt]) = 0.0;
//        
//        for (CFuint iSol = 0; iSol < m_nbrSolDep; ++iSol)
//        {
//          const CFuint solIdx = (*m_flxSolDep)[currFlxIdx][iSol];
//            
//          *(m_cellStatesFlxPnt[0][iFlxPnt]) += (*m_solPolyValsAtFlxPnts)[currFlxIdx][solIdx]*(*((*m_states[side])[solIdx]));
//        }
//      }
//      // compute ghost states with pert
//      (*m_bcStateComputers)[(*m_faceBCIdx[side])[faceIdx]]->setFace(((*m_faces[side])[faceIdx]));
//
//      (*m_bcStateComputers)[(*m_faceBCIdx[side])[faceIdx]]->computeGhostStates(m_cellStatesFlxPnt[0],m_flxPntGhostSol,m_unitNormalFlxPnts,m_flxPntCoords);
//      
//      (*(m_cellStatesFlxPnt[0][pertFlxPnt]))[m_pertVar] -= eps * (*m_solPolyValsAtFlxPnts)[pertFlxPntIdx][m_pertSol];
//      // compute ghost states without pert
//      (*m_bcStateComputers)[(*m_faceBCIdx[side])[faceIdx]]->computeGhostStates(m_cellStatesFlxPnt[0],m_cellStatesFlxPnt[1],m_unitNormalFlxPnts,m_flxPntCoords);
//      
//      for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
//      {
//        CFreal temp = (*(m_flxPntGhostSol[pertFlxPnt]))[iEq]-(*(m_cellStatesFlxPnt[1][pertFlxPnt]))[iEq];
//        
//        if (iEq == m_pertVar) temp -= eps*(*m_solPolyValsAtFlxPnts)[pertFlxPntIdx][m_pertSol];
//        
//        ///@todo check if faceLocalDir is ok & faceFlxPntConn
//        m_projectedCorrL = 0.5*temp*(m_faceJacobVecAbsSizeFlxPnts[pertFlxPnt]*(*m_faceLocalDir)[faceIdx])*m_unitNormalFlxPnts[pertFlxPnt];
//        
//        // Loop over solution pnts to calculate the grad updates
//        for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolDep; ++iSolPnt)
//        {
//          const CFuint iSolIdx = (*m_flxSolDep)[pertFlxPntIdx][iSolPnt];
//        
//          // inverse Jacobian determinant
//          const CFreal invJacobDet = 1.0/m_solJacobDet[side][iSolIdx];
//
//          /// @todo Check if this is also OK for triangles!!
//          (*m_cellGrads[side][iSolIdx])[iEq] += m_projectedCorrL*m_corrFctDiv[iSolIdx][pertFlxPntIdx]*invJacobDet;
////          if (m_cells[m_pertSide]->getID() == 1) 
////	  {
////          RealVector temp = m_projectedCorrL*m_corrFctDiv[iSolIdx][pertFlxPntIdx]*invJacobDet;
////            CFLog(INFO,"Ana Bnd: " << iSolIdx << ", "  << temp << "\n");
////          }
//        }
//      }
//    }
    if (!(*m_isFaceOnBoundary[side])[faceIdx])
    {
      // Get orientation of face
      const CFuint orient = (*m_faceOrients[side])[faceIdx];
        
      // cell side with respect to this face
      const CFuint cellSide = (*m_currCellSide[side])[faceIdx];
      
      // compute face Jacobian vectors
      m_faceJacobVecs = (*m_faces[side])[faceIdx]->computeFaceJacobDetVectorAtMappedCoords(*m_flxLocalCoords);
      
      // Loop over flux points to set the normal vectors
      for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
      {
        const CFuint currFlxIdx = (*m_faceFlxPntConnPerOrient)[orient][cellSide][iFlxPnt];
      
        for (CFuint jFlxPnt = 0; jFlxPnt < m_nbrFlxDep; ++jFlxPnt)
        {
          if ((*m_solFlxDep)[m_pertSol][jFlxPnt] == currFlxIdx)
          {
            pertFlxPnt = iFlxPnt;
            pertFlxPntIdx = currFlxIdx;
            
            // get face Jacobian vector size
            m_faceJacobVecAbsSizeFlxPnts[iFlxPnt] = faceJacobVecSizeFaceFlxPnts[(*m_faces[side])[faceIdx]->getID()][iFlxPnt];

            // set unit normal vector
            m_unitNormalFlxPnts[iFlxPnt] = m_faceJacobVecs[iFlxPnt]/m_faceJacobVecAbsSizeFlxPnts[iFlxPnt];
        
            break;
          }
        }
      }
      
      ///@todo check if faceLocalDir is ok & faceFlxPntConn
      m_projectedCorrL = -0.5*eps*(*m_solPolyValsAtFlxPnts)[pertFlxPntIdx][m_pertSol]*(m_faceJacobVecAbsSizeFlxPnts[pertFlxPnt]*(*m_faceMappedCoordDir)[orient][cellSide])*m_unitNormalFlxPnts[pertFlxPnt];

      // Loop over solution pnts to calculate the grad updates
      for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolDep; ++iSolPnt)
      {
        const CFuint iSolIdx = (*m_flxSolDep)[pertFlxPntIdx][iSolPnt];
        
        // inverse Jacobian determinant
        const CFreal invJacobDet = 1.0/m_solJacobDet[side][iSolIdx];

        /// @todo Check if this is also OK for triangles!!
        (*m_cellGrads[side][iSolIdx])[m_pertVar] += m_projectedCorrL*m_corrFctDiv[iSolIdx][pertFlxPntIdx]*invJacobDet;
//          if (m_cells[m_pertSide]->getID() == 5) 
//	  {
//            RealVector temp = m_projectedCorrL*m_corrFctDiv[iSolIdx][pertFlxPntIdx]*invJacobDet;
//              CFLog(INFO,"Ana otherFace: " << iSolIdx << ", "  << temp << "\n");
//          }
      }
    }
  }
  
  // Add the contribution of the correction of the gradients for this face
  
  // compute face Jacobian vectors
  m_faceJacobVecs = m_face->computeFaceJacobDetVectorAtMappedCoords(*m_flxLocalCoords);
  
  CFuint pertFlxPntIdxOtherSide;
  const CFuint otherSide = (side == LEFT) ? RIGHT : LEFT;
      
  // Loop over flux points to set the normal vectors
  for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
  {
    const CFuint currFlxIdx = (*m_faceFlxPntConnPerOrient)[m_orient][side][iFlxPnt];
     
    for (CFuint jFlxPnt = 0; jFlxPnt < m_nbrFlxDep; ++jFlxPnt)
    {
      if ((*m_solFlxDep)[m_pertSol][jFlxPnt] == currFlxIdx)
      {
        pertFlxPnt = iFlxPnt;
        pertFlxPntIdx = currFlxIdx;
        pertFlxPntIdxOtherSide = (*m_faceFlxPntConnPerOrient)[m_orient][otherSide][iFlxPnt];
        
        break;
      }
    }
      
    // get face Jacobian vector size
    m_faceJacobVecAbsSizeFlxPnts[iFlxPnt] = faceJacobVecSizeFaceFlxPnts[m_face->getID()][iFlxPnt];

    // set unit normal vector
    m_unitNormalFlxPnts[iFlxPnt] = m_faceJacobVecs[iFlxPnt]/m_faceJacobVecAbsSizeFlxPnts[iFlxPnt];
  }
      
  ///@todo check if faceLocalDir is ok & faceFlxPntConn
  m_projectedCorrL = -0.5*eps*(*m_solPolyValsAtFlxPnts)[pertFlxPntIdx][m_pertSol]*(m_faceJacobVecAbsSizeFlxPnts[pertFlxPnt]*(*m_faceMappedCoordDir)[m_orient][side])*m_unitNormalFlxPnts[pertFlxPnt];
  m_projectedCorrR = 0.5*eps*(*m_solPolyValsAtFlxPnts)[pertFlxPntIdx][m_pertSol]*(m_faceJacobVecAbsSizeFlxPnts[pertFlxPnt]*(*m_faceMappedCoordDir)[m_orient][otherSide])*m_unitNormalFlxPnts[pertFlxPnt];

  // Loop over solution pnts to calculate the grad updates
  for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolDep; ++iSolPnt)
  {
    const CFuint iSolIdx = (*m_flxSolDep)[pertFlxPntIdx][iSolPnt];
    const CFuint iSolIdxOtherSide = (*m_flxSolDep)[pertFlxPntIdxOtherSide][iSolPnt];
      
    m_affectedSolPnts[otherSide][iSolIdxOtherSide] = true;
    
    // inverse Jacobian determinant
    const CFreal invJacobDet = 1.0/m_solJacobDet[side][iSolIdx];
    const CFreal invJacobDetOtherSide = 1.0/m_solJacobDet[otherSide][iSolIdxOtherSide];

    /// @todo Check if this is also OK for triangles!!
    (*m_cellGrads[side][iSolIdx])[m_pertVar] += m_projectedCorrL*m_corrFctDiv[iSolIdx][pertFlxPntIdx]*invJacobDet;
    (*m_cellGrads[otherSide][iSolIdxOtherSide])[m_pertVar] += m_projectedCorrR*m_corrFctDiv[iSolIdxOtherSide][pertFlxPntIdxOtherSide]*invJacobDetOtherSide;
  }
  
  // compute the states projected on order P-1
  computeProjStates(m_statesPMinOne, side);
  
  computeEpsilon0(side);
  
  computeSmoothness(side);
  
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
  
  const CFreal cellEps = m_cellEpsilons[m_cells[side]->getID()];
  
  // loop over flx pnts to extrapolate the states to the flux points
  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {   
    //if (m_solEpsilons[iSol] > 0.001) CFLog(INFO, "eps before: " << m_solEpsilons[iSol] << "\n");
    // reset the states in the flx pnts
    //m_solEpsilons[iSol] = 0.0;

    // loop over the sol pnts to compute the states and grads in the flx pnts
    for (CFuint iNode = 0; iNode < m_nbrCornerNodes; ++iNode)
    {
      // get node local index
      const CFuint nodeIdx = (*m_cellNodes)[iNode]->getLocalID();
      
      const CFreal nodeEps = m_nodeEpsilons[nodeIdx] + m_epsilon - cellEps;
      
      //m_solEpsilons[side][iSol] += m_nodePolyValsAtSolPnts[iSol][iNode]*nodeEps/m_nbNodeNeighbors[nodeIdx];
    }
    //if (m_solEpsilons[side][iSol] > 0.001) CFLog(INFO, "eps after: " << m_solEpsilons[iSol] << "\n");
  }
  
  // loop over flx pnts to extrapolate the states to the flux points
  for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
  {   
    // local flux point indices in the left and right cell
    const CFuint flxPntIdx = (*m_faceFlxPntConnPerOrient)[m_orient][side][iFlxPnt];
    //if (m_epsilonLR[side][iFlxPnt] > 0.001) CFLog(INFO, "eps before: " << m_epsilonLR[side][iFlxPnt] << "\n");
    // reset the states in the flx pnts
    m_epsBackUp = m_epsilonLR[side][iFlxPnt];
    //m_epsilonLR[side][iFlxPnt] = 0.0;
    
    m_cellNodes = m_cells[side]->getNodes();

    // loop over the sol pnts to compute the states and grads in the flx pnts
    for (CFuint iNode = 0; iNode < m_faceNodes->size(); ++iNode)
    {
      for (CFuint iNodeCell = 0; iNodeCell < m_nbrCornerNodes; ++iNodeCell)
      {
	if ((*m_faceNodes)[iNode]->getLocalID() == (*m_cellNodes)[iNodeCell]->getLocalID())
	{
          //const CFuint nodeIdx = (*m_faceNodes)[iNode]->getLocalID();
	  // get node local index
          const CFuint nodeIdx = (*m_cellNodesConn)(m_cells[side]->getID(),iNodeCell);
          
          const CFreal nodeEps = m_nodeEpsilons[nodeIdx] + m_epsilon - cellEps;
	    
          //m_epsilonLR[side][iFlxPnt] += m_nodePolyValsAtFlxPnts[flxPntIdx][iNodeCell]*nodeEps/m_nbNodeNeighbors[nodeIdx];
	}
      }
    }
    //if (m_epsilonLR[side][iFlxPnt] > 0.001) CFLog(INFO, "eps after: " << m_epsilonLR[side][iFlxPnt] << "\n");
  }
}

//////////////////////////////////////////////////////////////////////////////

void LLAVJacobFluxReconstruction::computeInterfaceFlxCorrection()
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

void LLAVJacobFluxReconstruction::setFaceData(CFuint faceID)
{
  DiffRHSJacobFluxReconstruction::setFaceData(faceID);
  
  m_faceNodes = m_face->getNodes();
  
  // loop over flx pnts to extrapolate the states to the flux points
  for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
  {   
//     m_epsilonLR[LEFT][iFlxPnt] = m_cellEpsilons[m_cells[LEFT]->getID()];
//     m_epsilonLR[RIGHT][iFlxPnt] = m_cellEpsilons[m_cells[RIGHT]->getID()];
    
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
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void LLAVJacobFluxReconstruction::computeWaveSpeedUpdates(vector< CFreal >& waveSpeedUpd)
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
      const CFreal rho = 1.0;
      const CFreal epsilon = 0.5*(m_epsilonLR[LEFT][iFlx]+m_epsilonLR[RIGHT][iFlx]);
      visc = epsilon;
      
      // transform update states to physical data to calculate eigenvalues
      waveSpeedUpd[iSide] += visc*jacobXJacobXIntCoef/m_cellVolume[iSide];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void LLAVJacobFluxReconstruction::computeDivDiscontFlx(vector< RealVector >& residuals)
{
  // reset the extrapolated fluxes
  for (CFuint iFlxPnt = 0; iFlxPnt < m_flxPntsLocalCoords->size(); ++iFlxPnt)
  {
    m_extrapolatedFluxes[iFlxPnt] = 0.0;
  }

  // Loop over solution points to calculate the discontinuous flux.
  for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolPnts; ++iSolPnt)
  { 
    for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
    {
      *(m_tempGrad[iVar]) = (*(m_cellGrads[0][iSolPnt]))[iVar];
    }

    // calculate the discontinuous flux projected on x, y, z-directions
    for (CFuint iDim = 0; iDim < m_dim+m_ndimplus; ++iDim)
    { 
      m_contFlx[iSolPnt][iDim] = 0.0;
    
      for (CFuint iDim2 = 0; iDim2 < m_dim; ++iDim2)
      {
        for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
        {
          m_contFlx[iSolPnt][iDim][iVar] += m_solEpsilons[m_pertSide][iSolPnt]*((*m_tempGrad[iVar])[iDim2])*m_cellFluxProjVects[iDim][iSolPnt][iDim2];
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

  const CFuint nbrFaces = m_cells[m_pertSide]->nbNeighborGeos();
  
  for (CFuint iFace = 0; iFace < nbrFaces; ++iFace)
  {
    if (!((*m_isFaceOnBoundary[m_pertSide])[iFace]))
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

      //m_face=((*m_faces[m_pertSide])[iFace]);

      m_faceNodes = (*m_faces[m_pertSide])[iFace]->getNodes();
      //m_face = (*m_faces)[iFace];
      m_cellNodes = m_cells[m_pertSide]->getNodes();

      const CFreal cellEps = m_cellEpsilons[m_cells[m_pertSide]->getID()];

      // get the datahandle of the update coefficients
      DataHandle<CFreal> updateCoeff = socket_updateCoeff.getDataHandle();
  
      vector< CFreal > faceJacobVecSizeFlxPnts;
      faceJacobVecSizeFlxPnts.resize(m_nbrFaceFlxPnts);

      // compute flux point coordinates
      for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
      {
        m_flxPntCoords[iFlx] = (*m_faces[m_pertSide])[iFace]->computeCoordFromMappedCoord((*m_flxLocalCoords)[iFlx]);	
      }

      // compute face Jacobian vectors
      vector< RealVector > faceJacobVecs = (*m_faces[m_pertSide])[iFace]->computeFaceJacobDetVectorAtMappedCoords(*m_flxLocalCoords);

      // get face Jacobian vector sizes in the flux points
      DataHandle< vector< CFreal > > faceJacobVecSizeFaceFlxPnts = socket_faceJacobVecSizeFaceFlxPnts.getDataHandle();

      // Loop over flux points to compute the unit normals
      for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
      {
        // get face Jacobian vector size
        CFreal faceJacobVecAbsSizeFlxPnts = faceJacobVecSizeFaceFlxPnts[(*m_faces[m_pertSide])[iFace]->getID()][iFlxPnt];

	// set face Jacobian vector size with sign depending on mapped coordinate direction
        faceJacobVecSizeFlxPnts[iFlxPnt] = faceJacobVecAbsSizeFlxPnts*((*m_faceLocalDir)[iFace]);

	// set unit normal vector
	m_unitNormalFlxPnts2[iFlxPnt] = (faceJacobVecs[iFlxPnt]/faceJacobVecAbsSizeFlxPnts);
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

	  *(m_cellStatesFlxPnt[0][iFlxPnt]) += coeff*(*((*(m_states[m_pertSide]))[solIdx]));
	  
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
        (*m_bcStateComputers)[(*m_faceBCIdx[m_pertSide])[iFace]]->setFace(((*m_faces[m_pertSide])[iFace]));

	      (*m_bcStateComputers)[(*m_faceBCIdx[m_pertSide])[iFace]]->computeGhostGradients(m_cellGradFlxPnt[0],m_flxPntGhostGrads,m_unitNormalFlxPnts2,m_flxPntCoords);
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
              const CFuint nodeIdx = (*m_cellNodesConn)(m_cells[m_pertSide]->getID(),iNodeCell);
              
              CFreal nodeEps = m_nodeEpsilons[nodeIdx];
              
              //if (m_jacob) nodeEps = nodeEps + m_epsilon - cellEps;
	  
              epsilon += m_nodePolyValsAtFlxPnts[currFlxIdx][iNodeCell]*nodeEps/m_nbNodeNeighbors[nodeIdx];
	    }
	  }
        }

        if (!m_jacob && m_addUpdCoeff)
	{
	  // adding updateCoeff
	  CFreal visc = 1.0;
  
          m_waveSpeedUpd[0] = 0.0;

          const CFreal jacobXJacobXIntCoef = faceJacobVecSizeFlxPnts[iFlxPnt]*
                                             faceJacobVecSizeFlxPnts[iFlxPnt]*
                                             (*m_faceIntegrationCoefs)[iFlxPnt]*
                                             m_cflConvDiffRatio;
          visc = epsilon;
      
          // transform update states to physical data to calculate eigenvalues
          m_waveSpeedUpd[0] += visc*jacobXJacobXIntCoef/m_cells[m_pertSide]->computeVolume();

          // loop over the sol pnts of both sides to update the wave speeds
          for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
          {
            const CFuint solID = (*m_states[m_pertSide])[iSol]->getLocalID();
            updateCoeff[solID] += m_waveSpeedUpd[0];
          }
	}

        // compute the average sol and grad to use the BR2 scheme
        for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
        {
	  if (m_cells[m_pertSide]->getID() == 1092) CFLog(VERBOSE, "var: " << iVar << ", grad: " << *(m_cellGradFlxPnt[0][iFlxPnt][iVar]) << ", ghost: " << *(m_flxPntGhostGrads[iFlxPnt][iVar]) << "\n");
          *(m_avgGrad[iVar]) = (*(m_cellGradFlxPnt[0][iFlxPnt][iVar]) + *(m_flxPntGhostGrads[iFlxPnt][iVar]))/2.0;
        }
              
        m_flxPntRiemannFlux[iFlxPnt] = 0.0;
	      
        for (CFuint iDim = 0; iDim < m_dim; ++iDim)
        {
          for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
          {
            m_flxPntRiemannFlux[iFlxPnt][iVar] += epsilon*((*(m_avgGrad[iVar]))[iDim])*m_unitNormalFlxPnts2[iFlxPnt][iDim];
	    if (m_cells[m_pertSide]->getID() == 1092) CFLog(VERBOSE, "avgrad: " << (*(m_avgGrad[iVar]))[iDim] << "\n");
          }
        }

        // compute FI in the mapped coord frame
        m_cellFlx[0][iFlxPnt] = (m_flxPntRiemannFlux[iFlxPnt])*faceJacobVecSizeFlxPnts[iFlxPnt]; 
	if (m_cells[m_pertSide]->getID() == 1092) CFLog(VERBOSE, "riemannunit: " << m_flxPntRiemannFlux[iFlxPnt] << "jacob: " << faceJacobVecSizeFlxPnts[iFlxPnt] << "\n");
	
	for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolDep; ++iSolPnt)
        {  
	  const CFuint solIdx = (*m_flxSolDep)[currFlxIdx][iSolPnt];
          const CFreal divh = m_corrFctDiv[solIdx][currFlxIdx];
  
          // Fill in the corrections
          for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
          {
            residuals[solIdx][iVar] += (m_cellFlx[0][iFlxPnt][iVar] - m_extrapolatedFluxes[currFlxIdx][iVar]) * divh; 
	    if (m_cells[m_pertSide]->getID() == 1092) CFLog(VERBOSE, "riemann: " << m_cellFlx[0][iFlxPnt][iVar] << ", extr: " << m_extrapolatedFluxes[currFlxIdx][iVar] << "\n");
          }
        }
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void LLAVJacobFluxReconstruction::computeDivDiscontFlxNeighb(RealVector& residuals, const CFuint side)
{
  // reset the extrapolated fluxes
  for (CFuint iFlxPnt = 0; iFlxPnt < m_flxPntsLocalCoords->size(); ++iFlxPnt)
  {
    m_extrapolatedFluxes[iFlxPnt] = 0.0;
  }

  // Loop over solution points to calculate the discontinuous flux.
  for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolPnts; ++iSolPnt)
  { 
    if (m_affectedSolPnts[side][iSolPnt])
    {    
      for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
      {
        *(m_tempGrad[iVar]) = (*(m_cellGrads[side][iSolPnt]))[iVar];
      }
    
      // calculate the discontinuous flux projected on x, y, z-directions
      for (CFuint iDim = 0; iDim < m_dim+m_ndimplus; ++iDim)
      { 
        m_contFlxNeighb[side][iSolPnt][iDim] = 0.0;
    
        for (CFuint iDim2 = 0; iDim2 < m_dim; ++iDim2)
        {
          for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
          {
            m_contFlxNeighb[side][iSolPnt][iDim][iVar] += m_solEpsilons[side][iSolPnt]*((*m_tempGrad[iVar])[iDim2])*m_neighbCellFluxProjVects[side][iDim][iSolPnt][iDim2];
          }
        }
      }
    }
    
    for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFlxDep; ++iFlxPnt)
    {
      const CFuint flxIdx = (*m_solFlxDep)[iSolPnt][iFlxPnt];
      const CFuint dim = (*m_flxPntFlxDim)[flxIdx];
     
      m_extrapolatedFluxes[flxIdx] += (*m_solPolyValsAtFlxPnts)[flxIdx][iSolPnt]*(m_contFlxNeighb[side][iSolPnt][dim]);
    }
  }

  // Loop over solution pnts to calculate the divergence of the discontinuous flux
  for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolPnts; ++iSolPnt)
  {
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      // reset the divergence of FC
      residuals[m_nbrEqs*iSolPnt+iEq] = 0.0;
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
          residuals[m_nbrEqs*iSolPnt+iEq] += polyCoef*(m_contFlxNeighb[side][jSolIdx][iDir+m_ndimplus][iEq]);
	}
      }
    }
  }

  const CFuint nbrFaces = m_cells[side]->nbNeighborGeos();
  for (CFuint iFace = 0; iFace < nbrFaces; ++iFace)
  {
    if (!((*m_isFaceOnBoundary[side])[iFace]))
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
            residuals[m_nbrEqs*solIdx+iVar] += -m_extrapolatedFluxes[currFlxIdx][iVar] * divh; 
          }
        }
      }
    }
    else
    {

      m_faceNodes = (*m_faces[side])[iFace]->getNodes();
      //m_face = (*m_faces[side])[iFace];
      //(*m_bcStateComputers)[(*m_faceBCIdx[side])[iFace]]->setFace(m_face);
      m_cellNodes = m_cells[side]->getNodes();
    

      const CFreal cellEps = m_cellEpsilons[m_cells[side]->getID()];
      
      // get the datahandle of the update coefficients
      DataHandle<CFreal> updateCoeff = socket_updateCoeff.getDataHandle();
  
      vector< CFreal > faceJacobVecSizeFlxPnts;
      faceJacobVecSizeFlxPnts.resize(m_nbrFaceFlxPnts);

      // compute flux point coordinates
      for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
      {
        m_flxPntCoords[iFlx] = (*m_faces[side])[iFace]->computeCoordFromMappedCoord((*m_flxLocalCoords)[iFlx]);	
      }

      // compute face Jacobian vectors
      vector< RealVector > faceJacobVecs = (*m_faces[side])[iFace]->computeFaceJacobDetVectorAtMappedCoords(*m_flxLocalCoords);

      // get face Jacobian vector sizes in the flux points
      DataHandle< vector< CFreal > > faceJacobVecSizeFaceFlxPnts = socket_faceJacobVecSizeFaceFlxPnts.getDataHandle();

      // Loop over flux points to compute the unit normals
      for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
      {
        // get face Jacobian vector size
        CFreal faceJacobVecAbsSizeFlxPnts = faceJacobVecSizeFaceFlxPnts[(*m_faces[side])[iFace]->getID()][iFlxPnt];

	// set face Jacobian vector size with sign depending on mapped coordinate direction
        faceJacobVecSizeFlxPnts[iFlxPnt] = faceJacobVecAbsSizeFlxPnts*((*m_faceLocalDir)[iFace]);

	// set unit normal vector
	m_unitNormalFlxPnts2[iFlxPnt] = (faceJacobVecs[iFlxPnt]/faceJacobVecAbsSizeFlxPnts);
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

	  *(m_cellStatesFlxPnt[0][iFlxPnt]) += coeff*(*((*(m_states[side]))[solIdx]));
	  
          for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
          {
            *(m_cellGradFlxPnt[0][iFlxPnt][iVar]) += coeff*((*(m_cellGrads[side][solIdx]))[iVar]);
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

        (*m_bcStateComputers)[(*m_faceBCIdx[side])[iFace]]->setFace(((*m_faces[side])[iFace]));


	      (*m_bcStateComputers)[(*m_faceBCIdx[side])[iFace]]->computeGhostGradients(m_cellGradFlxPnt[0],m_flxPntGhostGrads,m_unitNormalFlxPnts2,m_flxPntCoords);
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
              const CFuint nodeIdx = (*m_cellNodesConn)(m_cells[side]->getID(),iNodeCell);

              CFreal nodeEps = m_nodeEpsilons[nodeIdx];
              
              //if (m_jacob) nodeEps = nodeEps + m_epsilon - cellEps;
	  
              epsilon += m_nodePolyValsAtFlxPnts[currFlxIdx][iNodeCell]*nodeEps/m_nbNodeNeighbors[nodeIdx];
	    }
	  }
        }

        if (!m_jacob && m_addUpdCoeff)
	{
	  // adding updateCoeff
	  CFreal visc = 1.0;
  
          m_waveSpeedUpd[0] = 0.0;

          const CFreal jacobXJacobXIntCoef = faceJacobVecSizeFlxPnts[iFlxPnt]*
                                             faceJacobVecSizeFlxPnts[iFlxPnt]*
                                             (*m_faceIntegrationCoefs)[iFlxPnt]*
                                             m_cflConvDiffRatio;
          visc = epsilon;
      
          // transform update states to physical data to calculate eigenvalues
          m_waveSpeedUpd[0] += visc*jacobXJacobXIntCoef/m_cells[side]->computeVolume();

          // loop over the sol pnts of both sides to update the wave speeds
          for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
          {
            const CFuint solID = (*m_states[side])[iSol]->getLocalID();
            updateCoeff[solID] += m_waveSpeedUpd[0];
          }
	}
    
        // compute the average sol and grad to use the BR2 scheme
        for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
        {
	  if (m_cells[side]->getID() == 1092) CFLog(VERBOSE, "var: " << iVar << ", grad: " << *(m_cellGradFlxPnt[0][iFlxPnt][iVar]) << ", ghost: " << *(m_flxPntGhostGrads[iFlxPnt][iVar]) << "\n");
          *(m_avgGrad[iVar]) = (*(m_cellGradFlxPnt[0][iFlxPnt][iVar]) + *(m_flxPntGhostGrads[iFlxPnt][iVar]))/2.0;
        }
              
        m_flxPntRiemannFlux[iFlxPnt] = 0.0;
	      
        for (CFuint iDim = 0; iDim < m_dim; ++iDim)
        {
          for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
          {
            m_flxPntRiemannFlux[iFlxPnt][iVar] += epsilon*((*(m_avgGrad[iVar]))[iDim])*m_unitNormalFlxPnts2[iFlxPnt][iDim];
	    if (m_cells[side]->getID() == 1092) CFLog(VERBOSE, "avgrad: " << (*(m_avgGrad[iVar]))[iDim] << "\n");
          }
        }
     
        // compute FI in the mapped coord frame
        m_cellFlx[0][iFlxPnt] = (m_flxPntRiemannFlux[iFlxPnt])*faceJacobVecSizeFlxPnts[iFlxPnt]; 
	if (m_cells[side]->getID() == 1092) CFLog(VERBOSE, "riemannunit: " << m_flxPntRiemannFlux[iFlxPnt] << "jacob: " << faceJacobVecSizeFlxPnts[iFlxPnt] << "\n");
	
	for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolDep; ++iSolPnt)
        {  
	  const CFuint solIdx = (*m_flxSolDep)[currFlxIdx][iSolPnt];
          const CFreal divh = m_corrFctDiv[solIdx][currFlxIdx];
  
          // Fill in the corrections
          for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
          {
            residuals[m_nbrEqs*solIdx+iVar] += (m_cellFlx[0][iFlxPnt][iVar] - m_extrapolatedFluxes[currFlxIdx][iVar]) * divh; 
	    if (m_cells[side]->getID() == 1092) CFLog(VERBOSE, "riemann: " << m_cellFlx[0][iFlxPnt][iVar] << ", extr: " << m_extrapolatedFluxes[currFlxIdx][iVar] << "\n");
          }
        }
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void LLAVJacobFluxReconstruction::setCellData(const CFuint side)
{ 
  m_cellNodes = m_cells[side]->getNodes();
  
  DataHandle< CFreal > artVisc = socket_artVisc.getDataHandle();
  
  // loop over flx pnts to extrapolate the states to the flux points
  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {   
//     m_solEpsilons[iSol] = m_cellEpsilons[m_cell->getID()];
    
    // reset the states in the flx pnts
    m_solEpsilons[side][iSol] = 0.0;

    // loop over the sol pnts to compute the states and grads in the flx pnts
    for (CFuint iNode = 0; iNode < m_nbrCornerNodes; ++iNode)
    {
      // get node local index
      //const CFuint nodeIdx = (*m_cellNodesConn)(m_elemIdx,iNode);
      
      const CFuint nodeIdx = (*m_cellNodes)[iNode]->getLocalID();
      
      m_solEpsilons[side][iSol] += m_nodePolyValsAtSolPnts[iSol][iNode]*m_nodeEpsilons[nodeIdx]/m_nbNodeNeighbors[nodeIdx];
    }
    
    artVisc[(((*m_states[side])[iSol]))->getLocalID()] = m_solEpsilons[side][iSol];
  }
}

//////////////////////////////////////////////////////////////////////////////

void LLAVJacobFluxReconstruction::computeProjStates(std::vector< RealVector >& projStates)
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

void LLAVJacobFluxReconstruction::computeProjStates(std::vector< RealVector >& projStates, const CFuint side)
{
  cf_assert(m_nbrSolPnts == projStates.size());
    
  if (m_order != 1)
  {
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      for (CFuint iSol = 0; iSol < projStates.size(); ++iSol)
      {
        m_tempSolPntVec[iSol] = (*((*m_states[side])[iSol]))[iEq];
      }

      m_tempSolPntVec = m_transformationMatrix*m_tempSolPntVec;

      for (CFuint iSol = 0; iSol < projStates.size(); ++iSol)
      {
        projStates[iSol][iEq] = m_tempSolPntVec[iSol];
      }
    }
  }
  else
  {
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      CFreal stateSum = 0.0;
      
      for (CFuint iSol = 0; iSol < projStates.size(); ++iSol)
      {
        stateSum += (*((*m_states[side])[iSol]))[iEq];
      }

      stateSum /= projStates.size();

      for (CFuint iSol = 0; iSol < projStates.size(); ++iSol)
      {
        projStates[iSol][iEq] = stateSum;
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void LLAVJacobFluxReconstruction::computeEpsilon()
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
  
  if (m_epsilon < 0.0 || m_epsilon != m_epsilon) 
  {
      CFLog(INFO, "eps: " << m_epsilon << ", eps0: " << m_epsilon0 << ", s: " << m_s << "\n");
      m_epsilon = 0.0;
  }
}

//////////////////////////////////////////////////////////////////////////////

void LLAVJacobFluxReconstruction::computeEpsilon0()
{ 
  // get the datahandle of the update coefficients
  DataHandle<CFreal> updateCoeff = socket_updateCoeff.getDataHandle();
  
  const CFreal wavespeed = updateCoeff[(*m_cellStates)[0]->getLocalID()];
  
  //const CFreal deltaKsi = 2.0/(m_order+2.0);
  
  const CFreal peclet = computePeclet();
  
  m_epsilon0 = max(wavespeed*(2.0/peclet - m_subcellRes/peclet),0.0);
}

//////////////////////////////////////////////////////////////////////////////

void LLAVJacobFluxReconstruction::computeEpsilon0(const CFuint side)
{ 
  // get the datahandle of the update coefficients
  DataHandle<CFreal> updateCoeff = socket_updateCoeff.getDataHandle();
  
  const CFreal wavespeed = updateCoeff[(*m_states[side])[0]->getLocalID()];
  
  //const CFreal deltaKsi = 2.0/(m_order+2.0);
  
  const CFreal peclet = computePeclet();
  
  m_epsilon0 = max(wavespeed*(2.0/peclet - m_subcellRes/peclet),0.0);
}

//////////////////////////////////////////////////////////////////////////////

CFreal LLAVJacobFluxReconstruction::computePeclet()
{
  return m_peclet;
}

//////////////////////////////////////////////////////////////////////////////

void LLAVJacobFluxReconstruction::computeSmoothness()
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
}

//////////////////////////////////////////////////////////////////////////////

void LLAVJacobFluxReconstruction::computeSmoothness(const CFuint side)
{ 
  CFreal sNum = 0.0;
  
  CFreal sDenom = 0.0;
  
  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {
    CFreal stateP = (*((*m_states[side])[iSol]))[m_monitoredVar];
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
  
  if (m_s > m_Smax)
  {
      m_Smax = m_s;
  }
}

//////////////////////////////////////////////////////////////////////////////

void LLAVJacobFluxReconstruction::storeEpsilon()
{
  for (CFuint iNode = 0; iNode < m_nbrCornerNodes; ++iNode)
  {
    // get node ID
    const CFuint nodeID = (*m_cellNodes)[iNode]->getLocalID();

    if (!m_useMax) 
    {
      const CFreal newEps = (1.0-m_dampingCoeff)*m_cellEpsilons[m_cell->getID()] + m_dampingCoeff*m_epsilon;
      m_nodeEpsilons[nodeID] += newEps;
      m_cellEpsilons[m_cell->getID()] = newEps;
      m_totalEps += newEps;
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
}

//////////////////////////////////////////////////////////////////////////////

void LLAVJacobFluxReconstruction::computeFlux(const RealVector& values, const std::vector< RealVector* >& gradients, const RealVector& normal, const CFreal& radius, RealVector& flux)
{
  //const CFreal epsilon = 0.5*(m_epsilonLR[LEFT][m_currFlx]+m_epsilonLR[RIGHT][m_currFlx]);
    
  const CFreal epsilon = m_solEpsilons[m_pertSide][m_pertSol];

  flux = 0.0;

  for (CFuint iDim = 0; iDim < m_dim; ++iDim)
  {
    for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
    {
      flux[iVar] += epsilon*((*(gradients[iVar]))[iDim])*normal[iDim];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void LLAVJacobFluxReconstruction::setup()
{
  CFAUTOTRACE;
  CFLog(VERBOSE, "LLAVJacob setup\n");
  // setup parent class
  DiffRHSJacobFluxReconstruction::setup();

  // get the update varset
  m_updateVarSet = getMethodData().getUpdateVar();
  
  // get cell-node connectivity
  m_cellNodesConn = MeshDataStack::getActive()->getConnectivity("cellNodes_InnerCells");
  
  // get the local spectral FD data
  vector< FluxReconstructionElementData* >& frLocalData = getMethodData().getFRLocalData();
  cf_assert(frLocalData.size() > 0);
  // for now, there should be only one type of element
  cf_assert(frLocalData.size() == 1);
  
  const CFPolyOrder::Type order = frLocalData[0]->getPolyOrder();
  
  m_subcellRes = frLocalData[0]->getSubcellResolution();
  
  m_order = static_cast<CFuint>(order);
  
  // get the coefs for extrapolation of the node artificial viscosities to the flx pnts
  m_nodePolyValsAtFlxPnts = frLocalData[0]->getNodePolyValsAtPnt(*(frLocalData[0]->getFlxPntsLocalCoords()));
  
  // get the coefs for extrapolation of the node artificial viscosities to the sol pnts
  m_nodePolyValsAtSolPnts = frLocalData[0]->getNodePolyValsAtPnt(*(frLocalData[0]->getSolPntsLocalCoords()));
  
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
  m_solEpsilons.resize(2);
  m_solEpsilons[LEFT].resize(m_nbrSolPnts);
  m_solEpsilons[RIGHT].resize(m_nbrSolPnts);
  m_epsilonLR.resize(2);
  m_epsilonLR[LEFT].resize(m_nbrFaceFlxPnts);
  m_epsilonLR[RIGHT].resize(m_nbrFaceFlxPnts);
  m_unitNormalFlxPnts2.resize(m_nbrFaceFlxPnts);
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
  
  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    m_unitNormalFlxPnts2[iFlx].resize(m_dim);
  }
  
  m_transformationMatrix.resize(m_nbrSolPnts,m_nbrSolPnts);
  
  m_transformationMatrix = (*vdm)*temp*(*vdmInv);
  
  //m_s0 = -m_s0*log10(static_cast<CFreal>(m_order));
  
  m_Smax = m_s0 + m_kappa;
  
  m_SmaxGlobal = m_Smax;
  
  m_nbNodeNeighbors = 0.0;
  
  m_flagComputeNbNghb = true;
  
  cf_assert(m_monitoredVar < m_nbrEqs);
}

//////////////////////////////////////////////////////////////////////////////

void LLAVJacobFluxReconstruction::unsetup()
{
  CFAUTOTRACE;
  
  // unsetup parent class
  DiffRHSJacobFluxReconstruction::unsetup();
}


//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD
