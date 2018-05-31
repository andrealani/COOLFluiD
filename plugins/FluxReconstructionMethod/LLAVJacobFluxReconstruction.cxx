// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
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
  m_nbPosPrev(),
  m_nbPosPrevGlobal(),
  m_subcellRes()
  {
    addConfigOptionsTo(this);
    
    m_kappa = 5.0;
    setParameter( "Kappa", &m_kappa);
    
    m_peclet = 2.0;
    setParameter( "Peclet", &m_peclet);
    
    m_freezeLimiterRes = -20.0;
    setParameter( "FreezeLimiterRes", &m_freezeLimiterRes);
  
    m_freezeLimiterIter = MathTools::MathConsts::CFuintMax();
    setParameter( "FreezeLimiterIter", &m_freezeLimiterIter);
    
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
  }
  
  
//////////////////////////////////////////////////////////////////////////////

void LLAVJacobFluxReconstruction::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal >("Kappa","Kappa factor of artificial viscosity.");
  
  options.addConfigOption< CFreal >("Peclet","Peclet number to be used for artificial viscosity.");
  
  options.addConfigOption< CFreal >("FreezeLimiterRes","Residual after which to freeze the residual.");
  
  options.addConfigOption< CFuint >("FreezeLimiterIter","Iteration after which to freeze the residual.");
  
  options.addConfigOption< bool >("AddPositivityPreservation","Bool telling whether extra viscosity needs to be added for positivity preservation.");
  
  options.addConfigOption< CFreal >("MinValue","Minimum value at which point positivity preservation is added.");
  
  options.addConfigOption< CFreal >("ViscFactor","Maximum factor applied to viscosity for positivity preservation.");
  
  options.addConfigOption< CFuint >("MonitoredVar","Index of the monitored var for positivity preservation.");
  
  options.addConfigOption< bool >("AddUpdateCoeff","Boolean telling whether the update coefficient based on the artificial flux is added.");
}

//////////////////////////////////////////////////////////////////////////////

void LLAVJacobFluxReconstruction::configure ( Config::ConfigArgs& args )
{
  FluxReconstructionSolverCom::configure(args);
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
      
      // if the states in the cell are parallel updatable, compute the resUpdates (-divFC)
      if ((*m_cellStates)[0]->isParUpdatable())
      {
	// compute the states projected on order P-1
	computeProjStates(m_statesPMinOne);
	
	// compute the artificial viscosity
	computeEpsilon();
	
	// store epsilon
	storeEpsilon();
      } 
      
      //release the GeometricEntity
      m_cellBuilder->releaseGE();
    }
  }
  
  const std::string nsp = this->getMethodData().getNamespace();
  
#ifdef CF_HAVE_MPI
    MPI_Comm comm = PE::GetPE().GetCommunicator(nsp);
    const CFuint count = 1;
    MPI_Allreduce(&m_totalEps, &m_totalEpsGlobal, count, MPI_DOUBLE, MPI_SUM, comm);
    MPI_Allreduce(&m_nbPosPrev, &m_nbPosPrevGlobal, count, MPI_UNSIGNED, MPI_SUM, comm);
#endif
    
  if (PE::GetPE().GetRank(nsp) == 0) 
  {
    // print total artificial viscosity and number of positivity preservations
    CFLog(INFO, "total eps: " << m_totalEpsGlobal << ", number of times positivity preserved: " << m_nbPosPrevGlobal << "\n");
  }

  PE::GetPE().setBarrier(nsp);
  
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
      if ((*m_states[LEFT ])[0]->isParUpdatable() || (*m_states[RIGHT])[0]->isParUpdatable())
      {
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
	
	m_jacob = true;

	// get all the faces neighbouring the cells
        m_faces[LEFT ] = m_cells[LEFT ]->getNeighborGeos();
        m_faces[RIGHT] = m_cells[RIGHT]->getNeighborGeos();

        // set the local indexes of the other faces than the current faces
        setOtherFacesLocalIdxs();

        // get the neigbouring states of the other faces
        setFaceNeighbourStates();

        // get the neigbouring gradients of the other faces
        setFaceNeighbourGradients();

	// make a back up of the grads
        for (CFuint iState = 0; iState < m_nbrSolPnts; ++iState)
        {
          for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
          {
            m_cellGradsBackUp[LEFT][iState][iVar] = (*m_cellGrads[LEFT][iState])[iVar];
            m_cellGradsBackUp[RIGHT][iState][iVar] = (*m_cellGrads[RIGHT][iState])[iVar];
          }
        }

	for (CFuint iSide = 0; iSide < 2; ++iSide)
        {
          // compute solution points Jacobian determinants
          m_solJacobDet[iSide] = m_cells[iSide]->computeGeometricShapeFunctionJacobianDeterminant(*m_solPntsLocalCoords);
	}

	if (!m_freezeGrads)
	{
          // compute auxiliary term for the perturbed gradient reconstruction from current face
          computeCellGradsMinusFaceTerm();

          // compute the unperturbed cell diffusive residuals
          computeUnpertCellDiffResiduals();
	}

        // compute the diffusive face term contribution to the jacobian
        if ((*m_states[LEFT])[0]->isParUpdatable() && (*m_states[RIGHT])[0]->isParUpdatable())
        {
	  if (!m_freezeGrads)
	  {
            // compute auxiliary terms for the perturbed gradient reconstruction from other faces
            computeCellGradsMinusOtherFaceTerms(LEFT );

            computeCellGradsMinusOtherFaceTerms(RIGHT);
	  }

          computeBothJacobsDiffFaceTerm();
        }
        else if ((*m_states[LEFT])[0]->isParUpdatable())
        {
	  if (!m_freezeGrads)
	  {
            // compute auxiliary terms for the perturbed gradient reconstruction from other faces
            computeCellGradsMinusOtherFaceTerms(RIGHT);
	    
	    computeCellGradsMinusOtherFaceTerms(LEFT );
	  }

          //computeOneJacobDiffFaceTerm(LEFT );
          
	  computeBothJacobsDiffFaceTerm();
        }
        else if ((*m_states[RIGHT])[0]->isParUpdatable())
        {
	  if (!m_freezeGrads)
	  {
            // compute auxiliary terms for the perturbed gradient reconstruction from other faces
            computeCellGradsMinusOtherFaceTerms(LEFT );
	    
	    computeCellGradsMinusOtherFaceTerms(RIGHT);
	  }

          //computeOneJacobDiffFaceTerm(RIGHT);
          
	  computeBothJacobsDiffFaceTerm();
        }


        // release the cells
        m_cellBuilders[LEFT ]->releaseGE();
        m_cellBuilders[RIGHT]->releaseGE();
      }
      
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
    for (CFuint elemIdx = startIdx; elemIdx < endIdx; ++elemIdx)
    {
      // build the GeometricEntity
      geoDataCell.idx = elemIdx;
      m_cell = m_cellBuilder->buildGE();

      // get the states in this cell
      m_cellStates = m_cell->getStates();
      
      // if the states in the cell are parallel updatable, compute the resUpdates (-divFC)
      if ((*m_cellStates)[0]->isParUpdatable())
      {
	// get the neighbouring faces
        m_facesCell = m_cell->getNeighborGeos();
	m_faces[0] = m_cell->getNeighborGeos();
	m_jacob = false;
      
	// set the cell data
	setCellData();

	// compute the divergence of the discontinuous flux (-divFD+divhFD)
	computeDivDiscontFlx(m_divContFlx);

	// update RHS
        updateRHS();
	
	m_jacob = true;

	// compute the contribution to the jacobian
        computeJacobDiffVolTerm();
      } 
      
      // divide by the Jacobian to transform the residuals back to the physical domain
      //divideByJacobDet();
      
//       for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
//       {
//         (*((*m_cellStates)[iSol]))[0] = m_solEpsilons[iSol];
//       }
      
      // print out the residual updates for debugging
      if(m_cell->getID() == 1944)
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
    if (m_cells[LEFT]->getID() == 1944) CFLog(VERBOSE, "FI: " << m_cellFlx[LEFT][iFlxPnt] << ", e: " << epsilon << ", grad: " << (*(m_avgGrad[0])) << "\n");
    if (m_cells[RIGHT]->getID() == 1944) CFLog(VERBOSE, "FI: " << m_cellFlx[RIGHT][iFlxPnt] << ", e: " << epsilon << ", grad: " << (*(m_avgGrad[0])) << "\n");
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
      visc = epsilon/rho;
      
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
    vector< RealVector > temp = *(m_cellGrads[0][iSolPnt]);
    vector< RealVector* > grad;
    grad.resize(m_nbrEqs);

    for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
    {
      cf_assert(temp.size() == m_nbrEqs);
      grad[iVar] = & (temp[iVar]);
    }
    
    // calculate the discontinuous flux projected on x, y, z-directions
    for (CFuint iDim = 0; iDim < m_dim; ++iDim)
    { 
      m_contFlx[iSolPnt][iDim] = 0.0;
    
      for (CFuint iDim2 = 0; iDim2 < m_dim; ++iDim2)
      {
        for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
        {
          m_contFlx[iSolPnt][iDim][iVar] += m_solEpsilons[iSolPnt]*((*(grad[iVar]))[iDim2])*m_cellFluxProjVects[iDim][iSolPnt][iDim2];
        }
      }
    }
    
    for (CFuint iFlxPnt = 0; iFlxPnt < m_flxPntsLocalCoords->size(); ++iFlxPnt)
    {
      CFuint dim = (*m_flxPntFlxDim)[iFlxPnt];
      m_extrapolatedFluxes[iFlxPnt] += (*m_solPolyValsAtFlxPnts)[iFlxPnt][iSolPnt]*(m_contFlx[iSolPnt][dim]);
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
          residuals[iSolPnt][iEq] += (*m_solPolyDerivAtSolPnts)[iSolPnt][iDir][jSolPnt]*(m_contFlx[jSolPnt][iDir][iEq]);

	  if (fabs(residuals[iSolPnt][iEq]) < MathTools::MathConsts::CFrealEps())
          {
            residuals[iSolPnt][iEq] = 0.0;
	  }
	}
      }
    }
  }
    
  const CFuint nbrFaces = m_cell->nbNeighborGeos();
  for (CFuint iFace = 0; iFace < nbrFaces; ++iFace)
  {
    if (!((*m_isFaceOnBoundaryCell)[iFace]))
    {
      for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolPnts; ++iSolPnt)
      {
        for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
        {
          const CFuint currFlxIdx = (*m_faceFlxPntConn)[iFace][iFlxPnt];
          const CFreal divh = m_corrFctDiv[iSolPnt][currFlxIdx];

          if (fabs(divh) > MathTools::MathConsts::CFrealEps())
          {   
            // Fill in the corrections
            for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
            {
              residuals[iSolPnt][iVar] += -m_extrapolatedFluxes[currFlxIdx][iVar] * divh; 
            }
          }
        }
      }
    }
    else
    {
      m_faceNodes = (*m_facesCell)[iFace]->getNodes();
      //m_face = (*m_faces)[iFace];
      m_cellNodes = m_cell->getNodes();
      
      // get the datahandle of the update coefficients
      DataHandle<CFreal> updateCoeff = socket_updateCoeff.getDataHandle();
	
      vector< RealVector > unitNormalFlxPnts;
  
      vector< CFreal > faceJacobVecSizeFlxPnts;
      faceJacobVecSizeFlxPnts.resize(m_nbrFaceFlxPnts);
	
      // get the local FR data
      vector< FluxReconstructionElementData* >& frLocalData = getMethodData().getFRLocalData();
    
      // compute flux point coordinates
      SafePtr< vector<RealVector> > flxLocalCoords = frLocalData[0]->getFaceFlxPntsFaceLocalCoords();
  
      // compute flux point coordinates
      for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
      {
        m_flxPntCoords[iFlx] = (*m_facesCell)[iFace]->computeCoordFromMappedCoord((*flxLocalCoords)[iFlx]);	
      }
          
      // compute face Jacobian vectors
      vector< RealVector > faceJacobVecs = (*m_facesCell)[iFace]->computeFaceJacobDetVectorAtMappedCoords(*flxLocalCoords);
  
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
        unitNormalFlxPnts.push_back(faceJacobVecs[iFlxPnt]/faceJacobVecAbsSizeFlxPnts);
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
        for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
        {
	  *(m_cellStatesFlxPnt[0][iFlxPnt]) += (*m_solPolyValsAtFlxPnts)[currFlxIdx][iSol]*(*((*(m_cellStates))[iSol]));
	  
          for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
          {
            *(m_cellGradFlxPnt[0][iFlxPnt][iVar]) += (*m_solPolyValsAtFlxPnts)[currFlxIdx][iSol]*((*(m_cellGrads[0][iSol]))[iVar]);
          }
        }
      }
      
      // compute ghost gradients
      (*m_bcStateComputers)[(*m_faceBCIdxCell)[iFace]]->computeGhostGradients(m_cellGradFlxPnt[0],m_flxPntGhostGrads,unitNormalFlxPnts,m_flxPntCoords);
      
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

          const CFreal jacobXJacobXIntCoef = faceJacobVecSizeFlxPnts[iFlxPnt]*
                                             faceJacobVecSizeFlxPnts[iFlxPnt]*
                                             (*m_faceIntegrationCoefs)[iFlxPnt]*
                                             m_cflConvDiffRatio;
          const CFreal rho = (*(m_cellStatesFlxPnt[0][iFlxPnt]))[0];
          visc = epsilon/rho;
      
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
            m_flxPntRiemannFlux[iFlxPnt][iVar] += epsilon*((*(m_avgGrad[iVar]))[iDim])*unitNormalFlxPnts[iFlxPnt][iDim];
	    if (m_cell->getID() == 1092) CFLog(VERBOSE, "avgrad: " << (*(m_avgGrad[iVar]))[iDim] << "\n");
          }
        }
     
        // compute FI in the mapped coord frame
        m_cellFlx[0][iFlxPnt] = (m_flxPntRiemannFlux[iFlxPnt])*faceJacobVecSizeFlxPnts[iFlxPnt]; 
	if (m_cell->getID() == 1092) CFLog(VERBOSE, "riemannunit: " << m_flxPntRiemannFlux[iFlxPnt] << "jacob: " << faceJacobVecSizeFlxPnts[iFlxPnt] << "\n");
	
	for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolPnts; ++iSolPnt)
        {  
          const CFreal divh = m_corrFctDiv[iSolPnt][currFlxIdx];

          if (fabs(divh) > MathTools::MathConsts::CFrealEps())
          {   
            // Fill in the corrections
            for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
            {
              residuals[iSolPnt][iVar] += (m_cellFlx[0][iFlxPnt][iVar] - m_extrapolatedFluxes[currFlxIdx][iVar]) * divh; 
	      if (m_cell->getID() == 1092) CFLog(VERBOSE, "riemann: " << m_cellFlx[0][iFlxPnt][iVar] << ", extr: " << m_extrapolatedFluxes[currFlxIdx][iVar] << "\n");
            }
          }
        }
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void LLAVJacobFluxReconstruction::setCellData()
{
  DiffRHSJacobFluxReconstruction::setCellData();
  
  m_cellNodes = m_cell->getNodes();
  
  // loop over flx pnts to extrapolate the states to the flux points
  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {   
//     m_solEpsilons[iSol] = m_cellEpsilons[m_cell->getID()];
    
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
}

//////////////////////////////////////////////////////////////////////////////

void LLAVJacobFluxReconstruction::computeProjStates(std::vector< RealVector >& projStates)
{
  if (m_order != 1)
  {
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      RealVector temp(projStates.size());
    
      for (CFuint iSol = 0; iSol < projStates.size(); ++iSol)
      {
        temp[iSol] = (*((*m_cellStates)[iSol]))[iEq];
      }

      RealVector tempProj(projStates.size());

      tempProj = m_transformationMatrix*temp;

      for (CFuint iSol = 0; iSol < projStates.size(); ++iSol)
      {
        projStates[iSol][iEq] = tempProj[iSol];
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
        stateSum += (*((*m_cellStates)[iSol]))[iEq];
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
}

//////////////////////////////////////////////////////////////////////////////

void LLAVJacobFluxReconstruction::computeEpsilon0()
{ 
  // get the datahandle of the update coefficients
  DataHandle<CFreal> updateCoeff = socket_updateCoeff.getDataHandle();
  
  const CFreal wavespeed = updateCoeff[(*m_cellStates)[0]->getLocalID()];
  
  //const CFreal deltaKsi = 2.0/(m_order+2.0);
  
  const CFreal peclet = computePeclet();
  
  m_epsilon0 = wavespeed*(2.0/peclet - m_subcellRes/peclet);
  
  if (m_addPosPrev) addPositivityPreservation();
}

//////////////////////////////////////////////////////////////////////////////

CFreal LLAVJacobFluxReconstruction::computePeclet()
{
  return m_peclet;
}

//////////////////////////////////////////////////////////////////////////////

void LLAVJacobFluxReconstruction::addPositivityPreservation()
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
  
  if (true)  //minState < m_minValue
  {
    //m_epsilon0 *= m_minValue/minState;
    //m_epsilon0 *= (m_viscFactor - 1.0)/(m_minValue*m_minValue)*minState*minState - 2.0*(m_viscFactor - 1.0)/m_minValue*minState + m_viscFactor;
    //const CFreal factor = min(100.0*m_viscFactor, (m_viscFactor - 1.0)/(m_minValue*m_minValue)*posPrevValue*posPrevValue - 2.0*(m_viscFactor - 1.0)/m_minValue*posPrevValue + m_viscFactor);
    const CFreal a =  0.;//2.*m_minValue; //m_viscFactor*m_minValue;
    const CFreal k = (m_minValue + a)/(tan(3.1415927/(1.0 - m_viscFactor)*(1.1 - (1.0 + m_viscFactor)/2.0)));
    const CFreal factor = (1.0 - m_viscFactor)/3.1415927 * atan((posPrevValue + a)/k) + (1.0 + m_viscFactor)/2.0;
    m_epsilon0 *= factor;
    if (factor > 1.1) m_nbPosPrev++;
    //if (factor > 1.1) CFLog(INFO, "posPrev: " << posPrevValue << ", factor: " << factor << "\n");
    
//     // only needed to plot where extra viscosity is added
//     for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
//     {
//       (*((*(m_cellStates))[iSol]))[0] = 10.0;
//     }
  }
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
}

//////////////////////////////////////////////////////////////////////////////

void LLAVJacobFluxReconstruction::computeFlux(const RealVector& sol, const std::vector< RealVector* >& grad, const RealVector& normals, RealVector& flux)
{

  const CFreal epsilon = 0.5*(m_epsilonLR[LEFT][m_currFlx]+m_epsilonLR[RIGHT][m_currFlx]);
  
  flux = 0.0;
    
  for (CFuint iDim = 0; iDim < m_dim; ++iDim)
  {
    for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
    {
      flux[iVar] += epsilon*((*(grad[iVar]))[iDim])*normals[iDim];
    }
  }

}

//////////////////////////////////////////////////////////////////////////////

void LLAVJacobFluxReconstruction::computeUnpertCellDiffResiduals()
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
    for (CFuint iDim = 0; iDim < m_dim; ++iDim)
    {
      vector<CFuint> dimList;
      dimList.resize(m_nbrSolPnts);
      for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolPnts; ++iSolPnt)
      {
        dimList[iSolPnt] = iDim;
      }
      m_cellFluxProjVects[iDim] = m_cells[iSide]->computeMappedCoordPlaneNormalAtMappedCoords(dimList,*m_solPntsLocalCoords);
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

void LLAVJacobFluxReconstruction::computePertCellDiffResiduals(const CFuint side)
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
  for (CFuint iDim = 0; iDim < m_dim; ++iDim)
  {
    vector<CFuint> dimList;
    dimList.resize(m_nbrSolPnts);
    for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolPnts; ++iSolPnt)
    {
      dimList[iSolPnt] = iDim;
    }
    m_cellFluxProjVects[iDim] = m_cells[side]->computeMappedCoordPlaneNormalAtMappedCoords(dimList,*m_solPntsLocalCoords);
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

void LLAVJacobFluxReconstruction::setup()
{
  CFAUTOTRACE;
  
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
  
  m_nodeEpsilons.resize(nbrNodes);
  m_nbNodeNeighbors.resize(nbrNodes);
  m_cellEpsilons.resize(nbrCells);
  m_solEpsilons.resize(m_nbrSolPnts);
  m_epsilonLR.resize(2);
  m_epsilonLR[LEFT].resize(m_nbrFaceFlxPnts);
  m_epsilonLR[RIGHT].resize(m_nbrFaceFlxPnts);
  
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
  for (CFuint idx = 0; idx < (m_order)*(m_order); ++idx)
  {
    temp(idx,idx) = 1.0;
  }
  
  m_transformationMatrix.resize(m_nbrSolPnts,m_nbrSolPnts);
  
  m_transformationMatrix = (*vdm)*temp*(*vdmInv);
  
  m_s0 = -3.0*log10(static_cast<CFreal>(m_order));
  
  m_nbNodeNeighbors = 0.0;
  
  m_flagComputeNbNghb = true;
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

