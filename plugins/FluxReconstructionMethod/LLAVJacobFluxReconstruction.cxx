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
  m_s0(),
  m_s(),
  m_epsilon0(),
  m_kappa(),
  m_cellNodes(),
  m_nbrCornerNodes(),
  m_nodeEpsilons(),
  m_cellEpsilons(),
  m_epsilonLR()
  {
    addConfigOptionsTo(this);
  }
  
  
//////////////////////////////////////////////////////////////////////////////

void LLAVJacobFluxReconstruction::defineConfigOptions(Config::OptionList& options)
{
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

	// set the face data
	setFaceData(m_face->getID());//faceID

	// compute the left and right states and gradients in the flx pnts
	computeFlxPntStatesAndGrads();

	// compute FI
	computeInterfaceFlxCorrection();

	// compute the wave speed updates
        //computeWaveSpeedUpdates(m_waveSpeedUpd);

        // update the wave speed
        //updateWaveSpeed();

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
      
      // get the neighbouring faces
      m_faces[0] = m_cell->getNeighborGeos();

      // if the states in the cell are parallel updatable, compute the resUpdates (-divFC)
      if ((*m_cellStates)[0]->isParUpdatable())
      {
	// set the cell data
	setCellData();

	// compute the divergence of the discontinuous flux (-divFD+divhFD)
	computeDivDiscontFlx(m_divContFlx);
      
	// update RHS
        updateRHS();

	// compute the contribution to the jacobian
        computeJacobDiffVolTerm();
      } 
      
      // divide by the Jacobian to transform the residuals back to the physical domain
      //divideByJacobDet();
      
      // print out the residual updates for debugging
      if(m_cell->getID() == 191)
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
  const CFreal epsilon = 0.5*(m_epsilonLR[LEFT]+m_epsilonLR[RIGHT]);
  
  // Loop over the flux points to calculate FI
  for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
  { 
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
  
  m_epsilonLR[LEFT] = m_cellEpsilons[m_cells[LEFT]->getID()];

  m_epsilonLR[RIGHT] = m_cellEpsilons[m_cells[RIGHT]->getID()];
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

    // calculate the discontinuous flux projected on x, y, z-directions
    for (CFuint iDim = 0; iDim < m_dim; ++iDim)
    { 
      m_contFlx[iSolPnt][iDim] = 0.0;
    
      for (CFuint iDim2 = 0; iDim2 < m_dim; ++iDim2)
      {
        for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
        {
          m_contFlx[iSolPnt][iDim][iVar] += m_epsilon*((*(grad[iVar]))[iDim2])*m_cellFluxProjVects[iDim][iSolPnt][iDim2];
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

    for (CFuint iFlxPnt = 0; iFlxPnt < m_flxPntsLocalCoords->size(); ++iFlxPnt)
    {
      const CFreal divh = m_corrFctDiv[iSolPnt][iFlxPnt];

      if (fabs(divh) > MathTools::MathConsts::CFrealEps())
      {   
        // Fill in the corrections
        for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
        {
          residuals[iSolPnt][iVar] += -m_extrapolatedFluxes[iFlxPnt][iVar] * divh; 
        }
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void LLAVJacobFluxReconstruction::setCellData()
{
  DiffRHSJacobFluxReconstruction::setCellData();
  
  m_epsilon = m_cellEpsilons[m_cell->getID()];
}

//////////////////////////////////////////////////////////////////////////////

void LLAVJacobFluxReconstruction::computeProjStates(std::vector< RealVector >& projStates)
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

//////////////////////////////////////////////////////////////////////////////

void LLAVJacobFluxReconstruction::computeEpsilon()
{
  computeEpsilon0();
  
  computeSmoothness();
  
  m_kappa = 0.4;
  
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
  
  const CFreal deltaKsi = 1.0/(m_order+2.0);
  
  const CFreal Pe = 10.0;//2.0
  
  m_epsilon0 = wavespeed*(2.0/Pe - deltaKsi/Pe);

}

//////////////////////////////////////////////////////////////////////////////

void LLAVJacobFluxReconstruction::computeSmoothness()
{ 
  CFreal sNum = 0.0;
  
  CFreal sDenom = 0.0;
  
  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {
    CFreal stateP = (*((*m_cellStates)[iSol]))[0];
    CFreal diffStatesPPMinOne = stateP - m_statesPMinOne[iSol][0];
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
  m_cellEpsilons[m_cell->getID()] = m_epsilon;
  
  for (CFuint iNode = 0; iNode < m_nbrCornerNodes; ++iNode)
  {
    // get node ID
    const CFuint nodeID = (*m_cellNodes)[iNode]->getLocalID();

    m_nodeEpsilons[nodeID] += m_epsilon;
  }
}

//////////////////////////////////////////////////////////////////////////////

void LLAVJacobFluxReconstruction::computeFlux(const RealVector& sol, const std::vector< RealVector* >& grad, const RealVector& normals, RealVector& flux)
{

  const CFreal epsilon = 0.5*(m_epsilonLR[LEFT]+m_epsilonLR[RIGHT]);
  
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
    m_epsilon = m_epsilonLR[iSide];
    
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

    // compute the volume term
    computeDivDiscontFlx(m_pertDivContFlx[0]);

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

      if ((*m_isFaceOnBoundary[iSide])[faceIdx])
      {
        // compute the boundary face contribution to the diffusive residuals
	computeBndRes(iSide, faceIdx, m_pertDivContFlx[0]);

	// put the perturbed and unperturbed corrections in the correct format
        // using m_pertResUpdates because the values stored in m_resUpdates should be preserved
        for (CFuint iState = 0; iState < m_nbrSolPnts; ++iState)
        {
          for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
          {
            m_pertResUpdates[0][m_nbrEqs*iState+iVar] = m_pertDivContFlx[0][iState][iVar];
          }
        }

        // add boundary face term
        m_unpertCellDiffRes[iSide] += m_pertResUpdates[0];
      }
      else
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
  }
}

//////////////////////////////////////////////////////////////////////////////

void LLAVJacobFluxReconstruction::computePertCellDiffResiduals(const CFuint side)
{
  m_epsilon = m_epsilonLR[side];
  
  DiffRHSJacobFluxReconstruction::computePertCellDiffResiduals(side);
}

//////////////////////////////////////////////////////////////////////////////

void LLAVJacobFluxReconstruction::setup()
{
  CFAUTOTRACE;
  DiffRHSJacobFluxReconstruction::setup();

  // get the update varset
  m_updateVarSet = getMethodData().getUpdateVar();
  
  // get the local spectral FD data
  vector< FluxReconstructionElementData* >& frLocalData = getMethodData().getFRLocalData();
  cf_assert(frLocalData.size() > 0);
  // for now, there should be only one type of element
  cf_assert(frLocalData.size() == 1);
  
  const CFPolyOrder::Type order = frLocalData[0]->getPolyOrder();
  
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
  
  m_nodeEpsilons.resize(nbrNodes);
  m_cellEpsilons.resize(nbrCells);
  m_epsilonLR.resize(2);
  
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
}

//////////////////////////////////////////////////////////////////////////////

void LLAVJacobFluxReconstruction::unsetup()
{
  CFAUTOTRACE;
  
  DiffRHSJacobFluxReconstruction::unsetup();
}


//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD

