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

#include "FluxReconstructionMethod/LLAVDiffJacobFluxReconstruction.hh"
#include "FluxReconstructionMethod/DiffRHSJacobFluxReconstruction.hh"
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

MethodCommandProvider< LLAVDiffJacobFluxReconstruction,
		       FluxReconstructionSolverData,
		       FluxReconstructionModule >
LLAVDiffJacobFluxReconstructionFluxReconstructionProvider("LLAVDiffJacob");

//////////////////////////////////////////////////////////////////////////////
  
LLAVDiffJacobFluxReconstruction::LLAVDiffJacobFluxReconstruction(const std::string& name) :
  LLAVJacobFluxReconstruction(name),
  m_pData()
  {
    addConfigOptionsTo(this);
  }
  
  
//////////////////////////////////////////////////////////////////////////////

void LLAVDiffJacobFluxReconstruction::defineConfigOptions(Config::OptionList& options)
{
}

//////////////////////////////////////////////////////////////////////////////

void LLAVDiffJacobFluxReconstruction::configure ( Config::ConfigArgs& args )
{
  FluxReconstructionSolverCom::configure(args);
}  

//////////////////////////////////////////////////////////////////////////////

void LLAVDiffJacobFluxReconstruction::execute()
{
  CFAUTOTRACE;
  
  CFLog(VERBOSE, "LLAVDiffJacobFluxReconstruction::execute()\n");
  
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
  
  // get the geodata of the cell builders and set the TRS
  CellToFaceGEBuilder::GeoData& geoDataCBL = m_cellBuilders[LEFT]->getDataGE();
  geoDataCBL.trs = cells;
  CellToFaceGEBuilder::GeoData& geoDataCBR = m_cellBuilders[RIGHT]->getDataGE();
  geoDataCBR.trs = cells;
  
  m_nodeEpsilons = 0.0;
  
  const CFreal residual = SubSystemStatusStack::getActive()->getResidual();
  
  const CFuint iter = SubSystemStatusStack::getActive()->getNbIter();
  
  m_useMax = residual < m_freezeLimiterRes || iter > m_freezeLimiterIter;
  
  //m_peclet = m_minValue/pow(10,m_s0);
  m_Smax = -100.0;
  
  m_totalEps = 0.0;
  m_nbPosPrev = 0;
  
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
    
  if (PE::GetPE().GetRank(nsp) == 0) 
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
      
      cf_assert(m_cellVolume[LEFT] > 0.0);
      cf_assert(m_cellVolume[RIGHT] > 0.0);
      
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
        
        // compute needed cell contributions
        if (!m_cellFlags[cellIDL] && (*m_states[LEFT ])[0]->isParUpdatable())
        {
          DiffRHSJacobFluxReconstruction::computeUnpertCellDiffResiduals(LEFT);
          m_unpertAllCellDiffRes[cellIDL] = m_unpertCellDiffRes[LEFT];

	  // update RHS
	  updateRHSUnpertCell(LEFT);
        }
        if (!m_cellFlags[cellIDR] && (*m_states[RIGHT])[0]->isParUpdatable())
        {
          DiffRHSJacobFluxReconstruction::computeUnpertCellDiffResiduals(RIGHT);
          m_unpertAllCellDiffRes[cellIDR] = m_unpertCellDiffRes[RIGHT];

	  // update RHS
	  updateRHSUnpertCell(RIGHT);
        }
	
	m_jacob = true;
	
//      // if one of the neighbouring cells is parallel updatable, compute the correction flux
//      if ((*m_states[LEFT ])[0]->isParUpdatable() || (*m_states[RIGHT])[0]->isParUpdatable())
//      {

	// get all the faces neighbouring the cells
        m_faces[LEFT ] = m_cells[LEFT ]->getNeighborGeos();
        m_faces[RIGHT] = m_cells[RIGHT]->getNeighborGeos();

        // set the local indexes of the other faces than the current faces
        setOtherFacesLocalIdxs();

	// make a back up of the grads
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
          computeOneJacobDiffFaceTerm(LEFT);
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

void LLAVDiffJacobFluxReconstruction::computePerturbedGradients(const CFuint side)
{ 
  DiffRHSJacobFluxReconstruction::computePerturbedGradients(side);
  
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
  
  CFuint otherSide;
  side == LEFT ? otherSide = RIGHT : otherSide = LEFT;
  
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

void LLAVDiffJacobFluxReconstruction::computeInterfaceFlxCorrection()
{ 
  // Loop over the flux points to calculate FI
  for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
  { 
    const CFreal epsilon = 0.5*(m_epsilonLR[LEFT][iFlxPnt]+m_epsilonLR[RIGHT][iFlxPnt]);
    
    // compute the average sol and grad to use the BR2 scheme
    for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
    {
      *(m_avgGrad[iVar]) = (*(m_cellGradFlxPnt[LEFT][iFlxPnt][iVar]) + *(m_cellGradFlxPnt[RIGHT][iFlxPnt][iVar]))/2.0;
      
      m_avgSol[iVar] = ((*(m_cellStatesFlxPnt[LEFT][iFlxPnt]))[iVar] + (*(m_cellStatesFlxPnt[RIGHT][iFlxPnt]))[iVar])/2.0;
    }
    
    prepareFluxComputation();
     
    // compute the Riemann flux
    DiffRHSJacobFluxReconstruction::computeFlux(m_avgSol,m_avgGrad,m_unitNormalFlxPnts[iFlxPnt],0,m_flxPntRiemannFlux[iFlxPnt]);
    
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

void LLAVDiffJacobFluxReconstruction::computeWaveSpeedUpdates(vector< CFreal >& waveSpeedUpd)
{
  // compute the wave speed updates for the neighbouring cells
  cf_assert(waveSpeedUpd.size() == 2);
  CFreal visc = 1.0;
  
  if (m_addUpdCoeff)
  {
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
}

//////////////////////////////////////////////////////////////////////////////

void LLAVDiffJacobFluxReconstruction::computeDivDiscontFlx(vector< RealVector >& residuals)
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

    m_avgSol = *((*m_cellStates)[iSolPnt]->getData());
    
    prepareFluxComputation();

    // calculate the discontinuous flux projected on x, y, z-directions
    for (CFuint iDim = 0; iDim < m_dim+m_ndimplus; ++iDim)
    {
       DiffRHSJacobFluxReconstruction::computeFlux(m_avgSol,m_tempGrad,m_cellFluxProjVects[iDim][iSolPnt],0,m_contFlx[iSolPnt][iDim]);
    }
    
    // calculate the discontinuous flux projected on x, y, z-directions
    for (CFuint iDim = 0; iDim < m_dim+m_ndimplus; ++iDim)
    { 
      for (CFuint iDim2 = 0; iDim2 < m_dim; ++iDim2)
      {
        for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
        {
          m_contFlx[iSolPnt][iDim][iVar] += m_solEpsilons[iSolPnt]*((*(m_tempGrad[iVar]))[iDim2])*m_cellFluxProjVects[iDim][iSolPnt][iDim2];
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
      m_faceNodes = (*m_facesCell)[iFace]->getNodes();
      //m_face = (*m_faces)[iFace];
      m_cellNodes = m_cell->getNodes();
      
      const CFreal cellEps = m_cellEpsilons[m_cell->getID()];
      
      // get the datahandle of the update coefficients
      DataHandle<CFreal> updateCoeff = socket_updateCoeff.getDataHandle();
  
      vector< CFreal > faceJacobVecSizeFlxPnts;
      faceJacobVecSizeFlxPnts.resize(m_nbrFaceFlxPnts);

      // compute flux point coordinates
      for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
      {
        m_flxPntCoords[iFlx] = (*m_facesCell)[iFace]->computeCoordFromMappedCoord((*m_flxLocalCoords)[iFlx]);	
      }

      // compute face Jacobian vectors
      vector< RealVector > faceJacobVecs = (*m_facesCell)[iFace]->computeFaceJacobDetVectorAtMappedCoords(*m_flxLocalCoords);

      // get face Jacobian vector sizes in the flux points
      DataHandle< vector< CFreal > > faceJacobVecSizeFaceFlxPnts = socket_faceJacobVecSizeFaceFlxPnts.getDataHandle();

      // Loop over flux points to compute the unit normals
      for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
      {
        // get face Jacobian vector size
        CFreal faceJacobVecAbsSizeFlxPnts = faceJacobVecSizeFaceFlxPnts[(*m_facesCell)[iFace]->getID()][iFlxPnt];

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
          m_waveSpeedUpd[0] += visc*jacobXJacobXIntCoef/m_cell->computeVolume();

          // loop over the sol pnts of both sides to update the wave speeds
          for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
          {
            const CFuint solID = (*m_cellStates)[iSol]->getLocalID();
            updateCoeff[solID] += m_waveSpeedUpd[0];
          }
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
        m_cellFlx[0][iFlxPnt] = (m_flxPntRiemannFlux[iFlxPnt])*faceJacobVecSizeFlxPnts[iFlxPnt]; 
	if (m_cell->getID() == 1092) CFLog(VERBOSE, "riemannunit: " << m_flxPntRiemannFlux[iFlxPnt] << "jacob: " << faceJacobVecSizeFlxPnts[iFlxPnt] << "\n");
	
	for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolDep; ++iSolPnt)
        {  
	  const CFuint solIdx = (*m_flxSolDep)[currFlxIdx][iSolPnt];
          const CFreal divh = m_corrFctDiv[solIdx][currFlxIdx];
  
          // Fill in the corrections
          for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
          {
            residuals[solIdx][iVar] += (m_cellFlx[0][iFlxPnt][iVar] - m_extrapolatedFluxes[currFlxIdx][iVar]) * divh; 
	    if (m_cell->getID() == 1092) CFLog(VERBOSE, "riemann: " << m_cellFlx[0][iFlxPnt][iVar] << ", extr: " << m_extrapolatedFluxes[currFlxIdx][iVar] << "\n");
          }
        }
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void LLAVDiffJacobFluxReconstruction::computeFlux(const RealVector& values, const std::vector< RealVector* >& gradients, const RealVector& normal, const CFreal& radius, RealVector& flux)
{
  const CFreal epsilon = 0.5*(m_epsilonLR[LEFT][m_currFlx]+m_epsilonLR[RIGHT][m_currFlx]);

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

void LLAVDiffJacobFluxReconstruction::computeUnpertCellDiffResiduals()
{
  // put the perturbed and unperturbed corrections in the correct format
  for (CFuint iState = 0; iState < m_nbrSolPnts; ++iState)
  {
    for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
    {
      m_resUpdates[LEFT][m_nbrEqs*iState+iVar] = m_divContFlxL[iState][iVar];
      m_resUpdates[RIGHT][m_nbrEqs*iState+iVar] = m_divContFlxR[iState][iVar];
    }
  }

  for (CFuint iSide = 0; iSide < 2; ++iSide)
  {
    m_cellNodes = m_cells[iSide]->getNodes();

    // loop over sol pnts
    for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
    {   
      // reset the states in the flx pnts
      m_solEpsilons[iSol] = 0.0;

      // loop over the sol pnts to compute the states and grads in the flx pnts
      for (CFuint iNode = 0; iNode < m_nbrCornerNodes; ++iNode)
      {
        // get node local index
        //const CFuint nodeIdx = (*m_cellNodesConn)(m_elemIdx,iNode);
      
        const CFuint nodeIdx = (*m_cellNodes)[iNode]->getLocalID();

        m_solEpsilons[iSol] += m_nodePolyValsAtSolPnts[iSol][iNode]*m_nodeEpsilons[nodeIdx]/m_nbNodeNeighbors[nodeIdx];
      }
    }

    // create a list of the dimensions in which the deriv will be calculated
    for (CFuint iDim = 0; iDim < m_dim+m_ndimplus; ++iDim)
    {
      m_cellFluxProjVects[iDim] = m_cells[iSide]->computeMappedCoordPlaneNormalAtMappedCoords(m_dimList[iDim],*m_solPntsLocalCoords);
    }

    // set the states
    *m_cellStates = *(m_states[iSide]);

    // make a backup of the grads if necessary
    vector< vector< RealVector >* > gradsBackup;
    gradsBackup.resize(m_nbrSolPnts);
    if (iSide == RIGHT)
    {
      for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
      {
        gradsBackup[iSol] = m_cellGrads[0][iSol];
      }
      m_cellGrads[0] = m_cellGrads[1];
    }

    m_cell = m_cells[iSide];
    m_facesCell = m_cells[iSide]->getNeighborGeos();
    m_isFaceOnBoundaryCell = m_isFaceOnBoundary[iSide];
    m_faceBCIdxCell = m_faceBCIdx[iSide];
    // compute the volume term
    computeDivDiscontFlx(m_pertDivContFlx[0]);

    m_isFaceOnBoundaryCell = m_cellBuilder->getGeoBuilder()->getIsFaceOnBoundary();
    m_faceBCIdxCell        = m_cellBuilder->getGeoBuilder()->getFaceBCIdx       ();

    // put the unpert discontinuous diff residual in the correct format
    for (CFuint iState = 0; iState < m_nbrSolPnts; ++iState)
    {
      for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
      {
        m_unpertCellDiffRes[iSide][m_nbrEqs*iState+iVar] = m_pertDivContFlx[0][iState][iVar];
      }
    }

    // restore grads if necessary
    if (iSide == RIGHT)
    {
      for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
      {
        m_cellGrads[0][iSol] = gradsBackup[iSol];
      }
    }

    // current face term
    m_unpertCellDiffRes[iSide] += m_resUpdates[iSide];

    // other face terms
    const CFuint nbrOtherFaces = m_otherFaceLocalIdxs[iSide].size();
    for (CFuint iFace = 0; iFace < nbrOtherFaces; ++iFace)
    {
      // get local face index
      const CFuint faceIdx = m_otherFaceLocalIdxs[iSide][iFace];

      if (!(*m_isFaceOnBoundary[iSide])[faceIdx])
      {
        // compute the internal face contribution to the diffusive residuals
        // using m_pertResUpdates because the values stored in m_resUpdates should be preserved
        // cell side with respect to this face
        computeFaceRes(iSide, faceIdx, iFace, m_pertDivContFlx[0]);

        // put the perturbed and unperturbed corrections in the correct format
        for (CFuint iState = 0; iState < m_nbrSolPnts; ++iState)
        {
          for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
          {
            m_pertResUpdates[0][m_nbrEqs*iState+iVar] = m_pertDivContFlx[0][iState][iVar];
          }
        }

        // add internal face term
        m_unpertCellDiffRes[iSide] += m_pertResUpdates[0];
      }
    }
    if (m_cells[iSide]->getID() == 1944) CFLog(VERBOSE, "unpert res: " << m_unpertCellDiffRes[iSide] << "\n");
  }
}

//////////////////////////////////////////////////////////////////////////////

void LLAVDiffJacobFluxReconstruction::computePertCellDiffResiduals(const CFuint side)
{
//   for (CFuint iState = 0; iState < m_nbrSolPnts; ++iState)
//   {
//     m_solEpsilons[iState] = m_epsilonLR[side][0];
//   }
  
  m_cellNodes = m_cells[side]->getNodes();
  
  // loop over sol pnts
  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {   
    // reset the states in the flx pnts
    m_solEpsilons[iSol] = 0.0;

    // loop over the sol pnts to compute the states and grads in the flx pnts
    for (CFuint iNode = 0; iNode < m_nbrCornerNodes; ++iNode)
    {
      // get node local index
      //const CFuint nodeIdx = (*m_cellNodesConn)(m_elemIdx,iNode);
      
      const CFuint nodeIdx = (*m_cellNodes)[iNode]->getLocalID();
      
      m_solEpsilons[iSol] += m_nodePolyValsAtSolPnts[iSol][iNode]*m_nodeEpsilons[nodeIdx]/m_nbNodeNeighbors[nodeIdx];
    }
  }
  
  // create a list of the dimensions in which the deriv will be calculated
  for (CFuint iDim = 0; iDim < m_dim+m_ndimplus; ++iDim)
  {
    m_cellFluxProjVects[iDim] = m_cells[side]->computeMappedCoordPlaneNormalAtMappedCoords(m_dimList[iDim],*m_solPntsLocalCoords);
  }
  
  // set the states
  *m_cellStates = *(m_states[side]);

  // make a backup of the grads if necessary
  vector< vector< RealVector >* > gradsBackup;
  gradsBackup.resize(m_nbrSolPnts);
  if (side == RIGHT)
  {
    for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
    {
      gradsBackup[iSol] = m_cellGrads[0][iSol];
    }
    m_cellGrads[0] = m_cellGrads[1];
  }
  
  m_cell = m_cells[side];
  m_facesCell = m_cells[side]->getNeighborGeos();
  m_isFaceOnBoundaryCell = m_isFaceOnBoundary[side];
  m_faceBCIdxCell = m_faceBCIdx[side];

  // compute the volume term
  computeDivDiscontFlx(m_pertDivContFlx[0]);
  
  m_isFaceOnBoundaryCell = m_cellBuilder->getGeoBuilder()->getIsFaceOnBoundary();
  m_faceBCIdxCell        = m_cellBuilder->getGeoBuilder()->getFaceBCIdx       ();

  // put the perturbed and unperturbed corrections in the correct format
  for (CFuint iState = 0; iState < m_nbrSolPnts; ++iState)
  {
    for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
    {
      m_pertCellDiffRes[m_nbrEqs*iState+iVar] = m_pertDivContFlx[0][iState][iVar];
    }
  }

  // add current face diffusive fluxes (m_pertResUpdates is set outside this function)
  m_pertCellDiffRes += m_pertResUpdates[side];
  
  // restore grads if necessary
  if (side == RIGHT)
  {
    for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
    {
      m_cellGrads[0][iSol] = gradsBackup[iSol];
    }
  }

  // add other face diffusive fluxes
  const CFuint nbrOtherFaces = m_otherFaceLocalIdxs[side].size();
  for (CFuint iFace = 0; iFace < nbrOtherFaces; ++iFace)
  {
    // get local face index
    const CFuint faceIdx = m_otherFaceLocalIdxs[side][iFace];

    if (!(*m_isFaceOnBoundary[side])[faceIdx])
    {
      computeFaceRes(side, faceIdx, iFace, m_pertDivContFlx[0]);

      // put the perturbed and unperturbed corrections in the correct format
      for (CFuint iState = 0; iState < m_nbrSolPnts; ++iState)
      {
        for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
        {
          m_pertResUpdates[0][m_nbrEqs*iState+iVar] = m_pertDivContFlx[0][iState][iVar];
        }
      }

      // add the contribution to the diffusive residuals
      m_pertCellDiffRes += m_pertResUpdates[0];
    }
  }
  if (m_cells[side]->getID() == 1944) CFLog(VERBOSE, "pert res: " << m_pertCellDiffRes << "\n");
}

//////////////////////////////////////////////////////////////////////////////

void LLAVDiffJacobFluxReconstruction::setup()
{
  CFAUTOTRACE;
  CFLog(VERBOSE, "LLAVDiffJacob setup\n");
  // setup parent class
  LLAVJacobFluxReconstruction::setup();
}

//////////////////////////////////////////////////////////////////////////////

void LLAVDiffJacobFluxReconstruction::unsetup()
{
  CFAUTOTRACE;
  
  // unsetup parent class
  LLAVJacobFluxReconstruction::unsetup();
}


//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD

