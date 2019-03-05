// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
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

MethodCommandProvider< DiffRHSJacobFluxReconstruction,
		       FluxReconstructionSolverData,
		       FluxReconstructionModule >
diffRHSJacobFluxReconstructionProvider("DiffRHSJacob");

//////////////////////////////////////////////////////////////////////////////
  
DiffRHSJacobFluxReconstruction::DiffRHSJacobFluxReconstruction(const std::string& name) :
  DiffRHSFluxReconstruction(name),
  m_cellBuilders(),
  m_cellBuilder(CFNULL),
  m_lss(CFNULL),
  m_numJacob(CFNULL),
  m_acc(CFNULL),
  m_accSC(CFNULL),
  m_faces(),
  m_faceNghbrStates(),
  m_faceNghbrGrads(),
  m_pertResUpdates(),
  m_derivResUpdates(),
  m_gradUpdates(),
  m_pertGrads(),
  m_cellGradsMinusFaceTerm(),
  m_cellGradsMinusOtherFaceTerm(),
  m_unpertCellDiffRes(),
  m_pertCellDiffRes(),
  m_derivCellDiffRes(),
  m_solJacobDet(),
  m_otherFaceLocalIdxs(),
  m_isFaceOnBoundary(),
  m_nghbrCellSide(),
  m_currCellSide(),
  m_faceOrients(),
  m_faceBCIdx(),
  m_bcStateComputers(CFNULL),
  m_flxPntGhostSol(),
  m_divContFlxL(),
  m_divContFlxR(),
  m_resUpdates(),
  m_cellGradsBackUp(),
  m_pertDivContFlx(),
  m_pertCorrections(),
  m_isFaceOnBoundaryCell(CFNULL),
  m_nghbrCellSideCell(CFNULL),
  m_currCellSideCell(CFNULL),
  m_faceOrientsCell(CFNULL),
  m_faceBCIdxCell(CFNULL),
  m_flxPntGhostGrads(),
  m_currFlx(),
  m_dimList(),
  m_gradTermL(),
  m_gradTermR(),
  m_gradTermTemp(),
  m_gradTerm(),
  m_gradTermBefore(),
  m_projectedCorrL(),
  m_projectedCorrR(),
  m_pertSide(),
  m_pertSol(),
  m_pertVar(),
  m_cellGradFlxPntBackup(),
  m_eps(),
  m_cellFlags(),
  m_unpertAllCellDiffRes(),
  m_neighbCellFluxProjVects()
  {
  }

//////////////////////////////////////////////////////////////////////////////

void DiffRHSJacobFluxReconstruction::configure ( Config::ConfigArgs& args )
{
  DiffRHSFluxReconstruction::configure(args);
}  

//////////////////////////////////////////////////////////////////////////////

void DiffRHSJacobFluxReconstruction::execute()
{
  CFAUTOTRACE;
  
  CFLog(VERBOSE, "DiffRHSJacobFluxReconstruction::execute()\n");
  
  // get the elementTypeData
  SafePtr< vector<ElementTypeData> > elemType = MeshDataStack::getActive()->getElementTypeData();

  // get InnerCells TopologicalRegionSet
  SafePtr<TopologicalRegionSet> cells = MeshDataStack::getActive()->getTrs("InnerCells");

  // get the geodata of the geometric entity builder and set the TRS
  CellToFaceGEBuilder::GeoData& geoDataCell = m_cellBuilder->getDataGE();
  geoDataCell.trs = cells;
  
  // reset cell flags
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
        
        // compute needed cell contributions
        if (!m_cellFlags[cellIDL] && (*m_states[LEFT ])[0]->isParUpdatable())
        {
          computeUnpertCellDiffResiduals(LEFT);
          m_unpertAllCellDiffRes[cellIDL] = m_unpertCellDiffRes[LEFT];

	  // update RHS
	  updateRHSUnpertCell(LEFT);
        }
        if (!m_cellFlags[cellIDR] && (*m_states[RIGHT])[0]->isParUpdatable())
        {
          computeUnpertCellDiffResiduals(RIGHT);
          m_unpertAllCellDiffRes[cellIDR] = m_unpertCellDiffRes[RIGHT];

	  // update RHS
	  updateRHSUnpertCell(RIGHT);
        }

	// get all the faces neighbouring the cells
        m_faces[LEFT ] = m_cells[LEFT ]->getNeighborGeos();
        m_faces[RIGHT] = m_cells[RIGHT]->getNeighborGeos();

        // set the local indexes of the other faces than the current faces
        setOtherFacesLocalIdxs();

        // get the neigbouring states of the other faces
        //setFaceNeighbourStates();

        // get the neigbouring gradients of the other faces
        //setFaceNeighbourGradients();

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

        // compute auxiliary term for the perturbed gradient reconstruction from current face
        //computeCellGradsMinusFaceTerm();

        // compute the unperturbed cell diffusive residuals
        //computeUnpertCellDiffResiduals();

        // compute the diffusive face term contribution to the jacobian
        if ((*m_states[LEFT])[0]->isParUpdatable() && (*m_states[RIGHT])[0]->isParUpdatable())
        {
          // compute auxiliary terms for the perturbed gradient reconstruction from other faces
          //computeCellGradsMinusOtherFaceTerms(LEFT );
          //computeCellGradsMinusOtherFaceTerms(RIGHT);

          computeBothJacobsDiffFaceTerm();
        }
        else if ((*m_states[LEFT])[0]->isParUpdatable())
        {
          // compute auxiliary terms for the perturbed gradient reconstruction from other faces
          //computeCellGradsMinusOtherFaceTerms(RIGHT);
	    
	  //computeCellGradsMinusOtherFaceTerms(LEFT );

          computeOneJacobDiffFaceTerm(LEFT );
          
	  //computeBothJacobsDiffFaceTerm();
        }
        else if ((*m_states[RIGHT])[0]->isParUpdatable())
        {
          // compute auxiliary terms for the perturbed gradient reconstruction from other faces
          //computeCellGradsMinusOtherFaceTerms(LEFT );
	    
	  //computeCellGradsMinusOtherFaceTerms(RIGHT);

          computeOneJacobDiffFaceTerm(RIGHT);
          
	  //computeBothJacobsDiffFaceTerm();
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
  
  //// Loop over the elements to calculate the divergence of the continuous flux
  
//  // loop over element types, for the moment there should only be one
//  const CFuint nbrElemTypes = elemType->size();
//  cf_assert(nbrElemTypes == 1);
//  for (m_iElemType = 0; m_iElemType < nbrElemTypes; ++m_iElemType)
//  {
//    // get start and end indexes for this type of element
//    const CFuint startIdx = (*elemType)[m_iElemType].getStartIdx();
//    const CFuint endIdx   = (*elemType)[m_iElemType].getEndIdx();
//
//    // loop over cells
//    for (CFuint elemIdx = startIdx; elemIdx < endIdx; ++elemIdx)
//    {
//      // build the GeometricEntity
//      geoDataCell.idx = elemIdx;
//      m_cell = m_cellBuilder->buildGE();
//
//      // get the states in this cell
//      m_cellStates = m_cell->getStates();
//      
//      // get the neighbouring faces
//      m_faces[0] = m_cell->getNeighborGeos();
//
//      // if the states in the cell are parallel updatable, compute the resUpdates (-divFC)
//      if ((*m_cellStates)[0]->isParUpdatable())
//      {
//	// set the cell data
//	setCellData();
//
//	// compute the divergence of the discontinuous flux (-divFD+divhFD)
//	computeDivDiscontFlx(m_divContFlx);
//      
//	// update RHS
//        updateRHS();
//        
//	// compute the contribution to the jacobian
//        computeJacobDiffVolTerm();
//      } 
//      
//      // divide by the Jacobian to transform the residuals back to the physical domain
//      //divideByJacobDet();
//      
//      // print out the residual updates for debugging
//      if(m_cell->getID() == 191)
//      {
//	CFLog(VERBOSE, "ID  = " << (*m_cellStates)[0]->getLocalID() << "\n");
//        CFLog(VERBOSE, "Update = \n");
//        // get the datahandle of the rhs
//        DataHandle< CFreal > rhs = socket_rhs.getDataHandle();
//        for (CFuint iState = 0; iState < m_nbrSolPnts; ++iState)
//        {
//          CFuint resID = m_nbrEqs*( (*m_cellStates)[iState]->getLocalID() );
//          for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
//          {
//            CFLog(VERBOSE, "" << rhs[resID+iVar] << " ");
//          }
//          CFLog(VERBOSE,"\n");
//          DataHandle<CFreal> updateCoeff = socket_updateCoeff.getDataHandle();
//          CFLog(VERBOSE, "UpdateCoeff: " << updateCoeff[(*m_cellStates)[iState]->getLocalID()] << "\n");
//        }
//      }
//      //release the GeometricEntity
//      m_cellBuilder->releaseGE();
//    }
//  }
}

//////////////////////////////////////////////////////////////////////////////

void DiffRHSJacobFluxReconstruction::computeJacobDiffVolTerm()
{
  const CFuint nbrFaces = m_cell->nbNeighborGeos();
  for (CFuint iFace = 0; iFace < nbrFaces; ++iFace)
  {
    const CFuint nbrFaceNeighbours = (*(m_faces[0]))[iFace]->nbNeighborGeos();
    for (CFuint iSide = 0; iSide < nbrFaceNeighbours; ++iSide)
    {
      GeometricEntity* cell = (*(m_faces[0]))[iFace]->getNeighborGeo(iSide);

      m_faceNghbrStates[0][iFace][iSide] = cell->getStates();
    }
  }

  // compute solution points Jacobian determinants
  m_solJacobDet[0] = m_cell->computeGeometricShapeFunctionJacobianDeterminant(*m_solPntsLocalCoords);

  // get residual factor
  const CFreal resFactor = getMethodData().getResFactor();

  // dereference accumulator
  BlockAccumulator& acc = *m_accSC;

  // set block row and column indices
  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {
    acc.setRowColIndex(iSol,(*m_cellStates)[iSol]->getLocalID());
  }

  // put the perturbed and unperturbed corrections in the correct format
  for (CFuint iState = 0; iState < m_nbrSolPnts; ++iState)
  {
    for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
    {
      m_resUpdates[0][m_nbrEqs*iState+iVar] = m_divContFlx[iState][iVar];
      m_cellGradsBackUp[0][iState][iVar] = (*m_cellGrads[0][iState])[iVar];
    }
  }

  // loop over the states in this cell to perturb the states
  for (m_pertSol = 0; m_pertSol < m_nbrSolPnts; ++m_pertSol)
  {
    // dereference state
    State& pertState = *(*m_cellStates)[m_pertSol];

    // loop over the variables in the state
    for (m_pertVar = 0; m_pertVar < m_nbrEqs; ++m_pertVar)
    {
      // perturb physical variable in state
      m_numJacob->perturb(m_pertVar,pertState[m_pertVar]);

      if (!m_freezeGrads)
      {
        // compute the perturbed cell gradients after the perturbation
        computePerturbedGradients();
      }

      // compute the perturbed residual updates (-divFD)
      computeDivDiscontFlx(m_pertCorrections);

      // put the perturbed and unperturbed corrections in the correct format
      for (CFuint iState = 0; iState < m_nbrSolPnts; ++iState)
      {
        for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
	{
          m_pertResUpdates[0][m_nbrEqs*iState+iVar] = m_pertCorrections[iState][iVar];
        }
      }

      // compute the finite difference derivative of the volume term
      m_numJacob->computeDerivative(m_pertResUpdates[0],m_resUpdates[0],m_derivResUpdates);

      // multiply residual update derivatives with residual factor
      m_derivResUpdates *= resFactor;
      if (m_cell->getID() == 1944) CFLog(VERBOSE, "deriv3: " << m_derivResUpdates << "\n");

      // add the derivative of the residual updates to the accumulator
      CFuint resUpdIdx = 0;
      for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol, resUpdIdx += m_nbrEqs)
      {
        acc.addValues(iSol,m_pertSol,m_pertVar,&m_derivResUpdates[resUpdIdx]);
      }

      // restore physical variable in state
      m_numJacob->restore(pertState[m_pertVar]);

      // restore the gradients in the sol pnts
      for (CFuint iState = 0; iState < m_nbrSolPnts; ++iState)
      {
        for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
        {
          (*m_cellGrads[0][iState])[iVar] = m_cellGradsBackUp[0][iState][iVar];
        }
      }
    }
  }
   //acc.printToScreen();

  if (getMethodData().doComputeJacobian())
  {
    // add the values to the jacobian matrix
    m_lss->getMatrix()->addValues(acc);
  }

  // reset to zero the entries in the block accumulator
  acc.reset();
}

//////////////////////////////////////////////////////////////////////////////

void DiffRHSJacobFluxReconstruction::computeBothJacobsDiffFaceTerm()
{
  // get residual factor
  const CFreal resFactor = getMethodData().getResFactor();

  // dereference accumulator
  BlockAccumulator& acc = *m_acc;

  // set block row and column indices
  CFuint solIdx = 0;
  for (CFuint iSide = 0; iSide < 2; ++iSide)
  {
    for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol, ++solIdx)
    {
      acc.setRowColIndex(solIdx,(*m_states[iSide])[iSol]->getLocalID());
    }
    
    for (CFuint iDim = 0; iDim < m_dim; ++iDim)
    {
      m_neighbCellFluxProjVects[iSide][iDim] = m_cells[iSide]->computeMappedCoordPlaneNormalAtMappedCoords(m_dimList[iDim],*m_solPntsLocalCoords);
    }
  }

  // loop over left and right cell
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

    // loop over the states to perturb the states
    for (m_pertSol = 0; m_pertSol < m_nbrSolPnts; ++m_pertSol)
    {
      // dereference state
      State& pertState = *(*m_states[m_pertSide])[m_pertSol];
      
      // Add the discontinuous gradient
      *m_cellStates = *(m_states[m_pertSide]);
  
      computeCellGradTerm(m_gradTermBefore);

      // loop over the variables in the state
      for (m_pertVar = 0; m_pertVar < m_nbrEqs; ++m_pertVar)
      {
        // perturb physical variable in state
        m_numJacob->perturb(m_pertVar,pertState[m_pertVar]);

        // compute the perturbed gradients in the current cell
        computePerturbedGradientsAnalytical(m_pertSide);
          
//        computePerturbedGradients(m_pertSide);
//
//        // compute the perturbed gradients in the other cell
//        computePertGradsFromFaceTerm(iOtherSide);         

	// compute the perturbed left and right states in the flx pnts
	computeFlxPntStatesAndGrads();

	// compute perturbed FI
	computeInterfaceFlxCorrection();
	
	// compute the perturbed corrections
	computeCorrection(m_pertSide, m_pertDivContFlx[m_pertSide]);
	computeCorrection(iOtherSide, m_pertDivContFlx[iOtherSide]);

        // put the perturbed and unperturbed corrections in the correct format
        for (CFuint iState = 0; iState < m_nbrSolPnts; ++iState)
        {
          for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
	  {
            m_pertResUpdates[m_pertSide][m_nbrEqs*iState+iVar] = m_pertDivContFlx[m_pertSide][iState][iVar];
	    m_pertResUpdates[iOtherSide][m_nbrEqs*iState+iVar] = m_pertDivContFlx[iOtherSide][iState][iVar];
          }
        }

        // update the perturbed cell residual
        
        // compute the finite difference derivative of the face term
        m_numJacob->computeDerivative(m_pertResUpdates[m_pertSide],m_resUpdates[m_pertSide],m_derivResUpdates);

        // multiply residual update derivatives with residual factor
        m_derivResUpdates *= resFactor;
	if (m_cells[m_pertSide]->getID() == 1) 
	{
	  CFLog(VERBOSE, "pert1: " << m_pertResUpdates[m_pertSide] << "\n");
	  CFLog(VERBOSE, "unpert1: " << m_resUpdates[m_pertSide] << "\n");
	  CFLog(VERBOSE, "deriv1: " << m_derivResUpdates << "\n");
	}

        // add the derivative of the residual updates to the accumulator
        CFuint resUpdIdx = 0;
        for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol, resUpdIdx += m_nbrEqs)
        {
          acc.addValues(iSol+pertSideTerm,m_pertSol+pertSideTerm,m_pertVar,&m_derivResUpdates[resUpdIdx]);
        }

        // compute the perturbed diffusive residual in the other cell
        computePertCellDiffResiduals(iOtherSide);
        
        m_derivResUpdates = m_resUpdates[iOtherSide] + m_unpertAllCellDiffRes[otherCellID];
        
        // update the perturbed cell residual
        // compute the finite difference derivative of the other cell the diffusive residual
        m_numJacob->computeDerivative(m_pertCellDiffRes,m_derivResUpdates,m_derivCellDiffRes);

        // multiply residual update derivatives with residual factor
        m_derivCellDiffRes *= resFactor;
	if (m_cells[m_pertSide]->getID() == 1) CFLog(VERBOSE, "deriv2: " << m_derivCellDiffRes << "\n");

        // add the derivative of the residual updates to the accumulator
        resUpdIdx = 0;
        for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol, resUpdIdx += m_nbrEqs)
        {
          acc.addValues(iSol+otherSideTerm,m_pertSol+pertSideTerm,m_pertVar,&m_derivCellDiffRes[resUpdIdx]);
        }
                
        // Add internal cell contributions if needed
        const CFuint cellID = m_cells[m_pertSide]->getID();
        if (!m_cellFlags[cellID])
        {
          computeDivDiscontFlxNeighb(m_pertCellDiffRes,m_pertSide);
          
          // update the perturbed cell residual
          // compute the finite difference derivative of the other cell the diffusive residual
          m_numJacob->computeDerivative(m_pertCellDiffRes,m_unpertAllCellDiffRes[cellID],m_derivCellDiffRes);

          // multiply residual update derivatives with residual factor
          m_derivCellDiffRes *= resFactor;

          // add the derivative of the residual updates to the accumulator
          resUpdIdx = 0;
          for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol, resUpdIdx += m_nbrEqs)
          {
            acc.addValues(iSol+pertSideTerm,m_pertSol+pertSideTerm,m_pertVar,&m_derivCellDiffRes[resUpdIdx]);
          }
        }

        // restore physical variable in state
        m_numJacob->restore(pertState[m_pertVar]);
	
	// restore the gradients in the sol pnts
        for (CFuint iState = 0; iState < m_nbrSolPnts; ++iState)
        {
          for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
          {
            (*m_cellGrads[m_pertSide][iState])[iVar] = m_cellGradsBackUp[m_pertSide][iState][iVar];
	    (*m_cellGrads[iOtherSide][iState])[iVar] = m_cellGradsBackUp[iOtherSide][iState][iVar];
          }
        }
      }
    }
  }
   //acc.printToScreen();

  if (getMethodData().doComputeJacobian())
  {
    // add the values to the jacobian matrix
    m_lss->getMatrix()->addValues(acc);
  }

  // reset to zero the entries in the block accumulator
  acc.reset();
}

//////////////////////////////////////////////////////////////////////////////

void DiffRHSJacobFluxReconstruction::computeOneJacobDiffFaceTerm(const CFuint side)
{
  // get residual factor
  const CFreal resFactor = getMethodData().getResFactor();

  // dereference accumulator
  BlockAccumulator& acc = *m_acc;

  // set block row and column indices
  CFuint solIdx = 0;
  for (CFuint iSide = 0; iSide < 2; ++iSide)
  {
    for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol, ++solIdx)
    {
      acc.setRowColIndex(solIdx,(*m_states[iSide])[iSol]->getLocalID());
    }
    
    for (CFuint iDim = 0; iDim < m_dim; ++iDim)
    {
      m_neighbCellFluxProjVects[iSide][iDim] = m_cells[iSide]->computeMappedCoordPlaneNormalAtMappedCoords(m_dimList[iDim],*m_solPntsLocalCoords);
    }
  }

  // loop over left and right cell
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

    // loop over the states to perturb the states
    for (m_pertSol = 0; m_pertSol < m_nbrSolPnts; ++m_pertSol)
    {
      // dereference state
      State& pertState = *(*m_states[m_pertSide])[m_pertSol];
      
      // Add the discontinuous gradient
      *m_cellStates = *(m_states[m_pertSide]);
  
      computeCellGradTerm(m_gradTermBefore);

      // loop over the variables in the state
      for (m_pertVar = 0; m_pertVar < m_nbrEqs; ++m_pertVar)
      {
        // perturb physical variable in state
        m_numJacob->perturb(m_pertVar,pertState[m_pertVar]);

        // compute the perturbed gradients in the current cell
        computePerturbedGradientsAnalytical(m_pertSide);
          
//        computePerturbedGradients(m_pertSide);
//
//        // compute the perturbed gradients in the other cell
//        computePertGradsFromFaceTerm(iOtherSide);         

	// compute the perturbed left and right states in the flx pnts
	computeFlxPntStatesAndGrads();

	// compute perturbed FI
	computeInterfaceFlxCorrection();
	
	// compute the perturbed corrections
	computeCorrection(side, m_pertDivContFlx[side]);

        // put the perturbed and unperturbed corrections in the correct format
        for (CFuint iState = 0; iState < m_nbrSolPnts; ++iState)
        {
          for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
	  {
            m_pertResUpdates[side][m_nbrEqs*iState+iVar] = m_pertDivContFlx[side][iState][iVar];
          }
        }

        // update the perturbed cell residual
        
        // compute the finite difference derivative of the face term
        if (m_pertSide == side)
        {
          m_numJacob->computeDerivative(m_pertResUpdates[m_pertSide],m_resUpdates[m_pertSide],m_derivResUpdates);

          // multiply residual update derivatives with residual factor
          m_derivResUpdates *= resFactor;
	  if (m_cells[m_pertSide]->getID() == 1) 
	  {
	    CFLog(VERBOSE, "pert1: " << m_pertResUpdates[m_pertSide] << "\n");
	    CFLog(VERBOSE, "unpert1: " << m_resUpdates[m_pertSide] << "\n");
	    CFLog(VERBOSE, "deriv1: " << m_derivResUpdates << "\n");
	  }

          // add the derivative of the residual updates to the accumulator
          CFuint resUpdIdx = 0;
          for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol, resUpdIdx += m_nbrEqs)
          {
            acc.addValues(iSol+pertSideTerm,m_pertSol+pertSideTerm,m_pertVar,&m_derivResUpdates[resUpdIdx]);
          }
          
          // Add internal cell contributions if needed
          const CFuint cellID = m_cells[m_pertSide]->getID();
          if (!m_cellFlags[cellID])
          {
            computeDivDiscontFlxNeighb(m_pertCellDiffRes,m_pertSide);
          
            // update the perturbed cell residual
            // compute the finite difference derivative of the other cell the diffusive residual
            m_numJacob->computeDerivative(m_pertCellDiffRes,m_unpertAllCellDiffRes[cellID],m_derivCellDiffRes);

            // multiply residual update derivatives with residual factor
            m_derivCellDiffRes *= resFactor;

            // add the derivative of the residual updates to the accumulator
            resUpdIdx = 0;
            for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol, resUpdIdx += m_nbrEqs)
            {
              acc.addValues(iSol+pertSideTerm,m_pertSol+pertSideTerm,m_pertVar,&m_derivCellDiffRes[resUpdIdx]);
            }
          }
        }
        else
        {
          // compute the perturbed diffusive residual in the other cell
          computePertCellDiffResiduals(iOtherSide);
          
          RealVector temp = m_resUpdates[iOtherSide] + m_unpertAllCellDiffRes[otherCellID];

          // update the perturbed cell residual
          // compute the finite difference derivative of the other cell the diffusive residual
          m_numJacob->computeDerivative(m_pertCellDiffRes,temp,m_derivCellDiffRes);

          // multiply residual update derivatives with residual factor
          m_derivCellDiffRes *= resFactor;
	  if (m_cells[m_pertSide]->getID() == 1) CFLog(VERBOSE, "deriv2: " << m_derivCellDiffRes << "\n");

          // add the derivative of the residual updates to the accumulator
          CFuint resUpdIdx = 0;
          for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol, resUpdIdx += m_nbrEqs)
          {
            acc.addValues(iSol+otherSideTerm,m_pertSol+pertSideTerm,m_pertVar,&m_derivCellDiffRes[resUpdIdx]);
          }
        }

        // restore physical variable in state
        m_numJacob->restore(pertState[m_pertVar]);
	
	// restore the gradients in the sol pnts
        for (CFuint iState = 0; iState < m_nbrSolPnts; ++iState)
        {
          for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
          {
            (*m_cellGrads[m_pertSide][iState])[iVar] = m_cellGradsBackUp[m_pertSide][iState][iVar];
	    (*m_cellGrads[iOtherSide][iState])[iVar] = m_cellGradsBackUp[iOtherSide][iState][iVar];
          }
        }
      }
    }
  }
   //acc.printToScreen();

  if (getMethodData().doComputeJacobian())
  {
    // add the values to the jacobian matrix
    m_lss->getMatrix()->addValues(acc);
  }

  // reset to zero the entries in the block accumulator
  acc.reset();
}

//////////////////////////////////////////////////////////////////////////////

void DiffRHSJacobFluxReconstruction::storeBackups()
{
//  for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
//  {
//    m_cellStatesFlxPntBackup[iFlxPnt] = *(m_cellStatesFlxPnt[iFlxPnt]);
//    
//    for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
//    {
//      *(m_cellGradFlxPntBackup[iFlxPnt][iVar]) = *(m_cellGradFlxPnt[iFlxPnt][iVar]);
//    }
//  }
  
  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  { 
    for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
    {
      m_cellGradsBackUp[LEFT][iSol][iVar] = (*m_cellGrads[LEFT][iSol])[iVar];
      m_cellGradsBackUp[RIGHT][iSol][iVar] = (*m_cellGrads[RIGHT][iSol])[iVar];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void DiffRHSJacobFluxReconstruction::restoreFromBackups()
{
//  for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
//  {
//    if (m_influencedFlxPnts[iFlxPnt])
//    {
//      (*(m_cellStatesFlxPnt[iFlxPnt]))[m_pertVar] = m_cellStatesFlxPntBackup[iFlxPnt][m_pertVar];
//    }
//      
//    *(m_cellGradFlxPnt[iFlxPnt][m_pertVar]) = *(m_cellGradFlxPntBackup[iFlxPnt][m_pertVar]);
//    
//  }
  
  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  { 
    for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
    {
      (*m_cellGrads[LEFT][iSol])[iVar] = m_cellGradsBackUp[LEFT][iSol][iVar];
      (*m_cellGrads[RIGHT][iSol])[iVar] = m_cellGradsBackUp[RIGHT][iSol][iVar];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void DiffRHSJacobFluxReconstruction::computePerturbedGradients(const CFuint side)
{ 
  for (CFuint iDim = 0; iDim < m_dim; ++iDim)
  {
    m_cellFluxProjVects[iDim] = m_cells[m_pertSide]->computeMappedCoordPlaneNormalAtMappedCoords(m_dimList[iDim],*m_solPntsLocalCoords);
  }
    
  // Add the discontinuous gradient
  *m_cellStates = *(m_states[side]);
  
  computeCellGradTerm(m_gradTerm);
  
  // Loop over solution pnts to calculate the grad updates
  for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolPnts; ++iSolPnt)
  {
    // Loop over  variables
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      //set the grad updates to 0 
      m_gradUpdates[side][iSolPnt][iEq] = 0.0;
    }
  }
  
  // Loop over solution pnts to calculate the grad updates
  for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolPnts; ++iSolPnt)
  {
    // Loop over  variables
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      // Loop over gradient directions
      for (CFuint iDir = 0; iDir < m_dim; ++iDir)
      {
        m_projectedCorrL = m_gradTerm(iEq,iSolPnt) * m_cellFluxProjVects[iDir][iSolPnt];
        
        // Loop over solution pnts to count factor of all sol pnt polys
        for (CFuint jSolPnt = 0; jSolPnt < m_nbrSolSolDep; ++jSolPnt)
        {
          const CFuint jSolIdx = (*m_solSolDep)[iSolPnt][jSolPnt];
	  
          // compute the grad updates
          m_gradUpdates[side][jSolIdx][iEq] += (*m_solPolyDerivAtSolPnts)[jSolIdx][iDir][iSolPnt]*m_projectedCorrL;
	}
      }
    }
  }
  
//  for (CFuint jSolPnt = 0; jSolPnt < m_nbrSolSolDep; ++jSolPnt)
//  { 
//    const CFuint jSolIdx = (*m_solSolDep)[m_pertSol][jSolPnt];
//
//    // update gradients
//    (*m_cellGrads[side][jSolIdx])[m_pertVar] = m_gradUpdates[side][jSolIdx][m_pertVar];
//  }
  
  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {
    // update gradients
    for (CFuint iGrad = 0; iGrad < m_nbrEqs; ++iGrad)
    {
      (*m_cellGrads[side][iSol])[iGrad] = m_gradUpdates[side][iSol][iGrad];
    }   
  }
  
  // get face Jacobian vector sizes in the flux points
  DataHandle< vector< CFreal > > faceJacobVecSizeFaceFlxPnts = socket_faceJacobVecSizeFaceFlxPnts.getDataHandle();
  
  // Add the contribution of the correction to the gradients for each face
  // compute other face contributions to the gradients
  const CFuint nbrOtherFaces = m_otherFaceLocalIdxs[side].size();
  for (CFuint iFace = 0; iFace < nbrOtherFaces; ++iFace)
  {
    // get local face index
    const CFuint faceIdx = m_otherFaceLocalIdxs[side][iFace];
    const CFuint orient = (*m_faceOrients[side])[faceIdx];
            
    // Loop over flux points to set the normal vectors
    for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
    {
      // get face Jacobian vector size
      m_faceJacobVecAbsSizeFlxPnts[iFlxPnt] = faceJacobVecSizeFaceFlxPnts[(*m_faces[side])[faceIdx]->getID()][iFlxPnt];
      
      // compute face Jacobian vectors
      m_faceJacobVecs = (*m_faces[side])[faceIdx]->computeFaceJacobDetVectorAtMappedCoords(*m_flxLocalCoords);

      // set unit normal vector
      m_unitNormalFlxPnts[iFlxPnt] = m_faceJacobVecs[iFlxPnt]/m_faceJacobVecAbsSizeFlxPnts[iFlxPnt];
      
      m_flxPntCoords[iFlxPnt] = (*m_faces[side])[faceIdx]->computeCoordFromMappedCoord((*m_flxLocalCoords)[iFlxPnt]);
    }

    if ((*m_isFaceOnBoundary[side])[faceIdx])
    {
      // Loop over flux points to extrapolate the states to the flux points
      for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
      {
        *(m_cellStatesFlxPnt[0][iFlxPnt]) = 0.0;

        const CFuint currFlxIdx = (*m_faceFlxPntConn)[faceIdx][iFlxPnt];
        
        for (CFuint iSol = 0; iSol < m_nbrSolDep; ++iSol)
        {
          const CFuint solIdx = (*m_flxSolDep)[currFlxIdx][iSol];
            
          *(m_cellStatesFlxPnt[0][iFlxPnt]) += (*m_solPolyValsAtFlxPnts)[currFlxIdx][solIdx]*(*((*m_states[side])[solIdx]));
        }
      }
  
      // compute ghost states
      (*m_bcStateComputers)[(*m_faceBCIdx[side])[faceIdx]]->computeGhostStates(m_cellStatesFlxPnt[0],m_flxPntGhostSol,m_unitNormalFlxPnts,m_flxPntCoords);
      
      computeBndGradTerms(m_gradTermL,m_gradTermR);
      
      // Loop over solution pnts to reset the grad updates
      for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolPnts; ++iSolPnt)
      {
        // Loop over  variables
        for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
        {
          //set the grad updates to 0 
          m_gradUpdates[side][iSolPnt][iEq] = 0.0;
        }
      }
      
      // compute the face corrections to the gradients
      for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
      {
        const CFuint flxIdx = (*m_faceFlxPntConn)[faceIdx][iFlx];
        
        // Loop over  variables
        for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
        {
          const CFreal avgSol = (m_gradTermL(iEq,iFlx)+m_gradTermR(iEq,iFlx))/2.0;
          ///@todo check if faceLocalDir is ok & faceFlxPntConn
	  m_projectedCorrL = (avgSol-m_gradTermL(iEq,iFlx))*(m_faceJacobVecAbsSizeFlxPnts[iFlx]*(*m_faceLocalDir)[faceIdx])*m_unitNormalFlxPnts[iFlx];

          // Loop over solution pnts to calculate the grad updates
          for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolDep; ++iSolPnt)
          {
            const CFuint iSolIdx = (*m_flxSolDep)[flxIdx][iSolPnt];

	    /// @todo Check if this is also OK for triangles!!
	    m_gradUpdates[side][iSolIdx][iEq] += m_projectedCorrL*m_corrFctDiv[iSolIdx][flxIdx];
          }
        }
      }
      
//      for (CFuint jSolPnt = 0; jSolPnt < m_nbrSolSolDep; ++jSolPnt)
//  { 
//    const CFuint jSolIdx = (*m_solSolDep)[m_pertSol][jSolPnt];
//
//    // update gradients
//    (*m_cellGrads[side][jSolIdx])[m_pertVar] += m_gradUpdates[side][jSolIdx][m_pertVar];
//  }

      // add the contribution to the gradients
      for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
      {
        for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
        {
          (*m_cellGrads[side][iSol])[iEq] += m_gradUpdates[side][iSol][iEq];
          
//          if (m_cells[m_pertSide]->getID() == 1) 
//	  {
//            RealVector temp = m_gradUpdates[side][iSol][iEq]/m_solJacobDet[side][iSol];
//            CFLog(INFO,"Num bndFace: " << iSol << ", "  << temp << "\n");
//          }
        }
      }
    }
    else
    {
      // cell side with respect to this face
      const CFuint cellSide = (*m_currCellSide[side])[faceIdx];
      
      // loop over flx pnts to extrapolate the states to the flux points and get the 
      // discontinuous normal flux in the flux points and the Riemann flux
      for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
      {     
        // local flux point indices in the left and right cell
        const CFuint flxPntIdxL = (*m_faceFlxPntConnPerOrient)[orient][LEFT][iFlxPnt];
        const CFuint flxPntIdxR = (*m_faceFlxPntConnPerOrient)[orient][RIGHT][iFlxPnt];
    
        *(m_cellStatesFlxPnt[LEFT][iFlxPnt]) = 0.0;
        *(m_cellStatesFlxPnt[RIGHT][iFlxPnt]) = 0.0;

        // extrapolate the left and right states to the flx pnts
        for (CFuint iSol = 0; iSol < m_nbrSolDep; ++iSol)
        {
          const CFuint solIdxL = (*m_flxSolDep)[flxPntIdxL][iSol];
          const CFuint solIdxR = (*m_flxSolDep)[flxPntIdxR][iSol];
      
          *(m_cellStatesFlxPnt[LEFT][iFlxPnt]) += (*m_solPolyValsAtFlxPnts)[flxPntIdxL][solIdxL]*(*(*m_faceNghbrStates[side][iFace][LEFT])[solIdxL]);
          *(m_cellStatesFlxPnt[RIGHT][iFlxPnt]) += (*m_solPolyValsAtFlxPnts)[flxPntIdxR][solIdxR]*(*(*m_faceNghbrStates[side][iFace][RIGHT])[solIdxR]);
        }
      }
      
      computeFaceGradTerms(m_gradTermL,m_gradTermR);
      
      // Loop over solution pnts to reset the grad updates
      for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolPnts; ++iSolPnt)
      {
        // Loop over  variables
        for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
        {
          //set the grad updates to 0 
          m_gradUpdates[cellSide][iSolPnt][iEq] = 0.0;
        }
      }
      
      // compute the face corrections to the gradients
      for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
      {
        const CFuint flxIdx = (*m_faceFlxPntConnPerOrient)[orient][cellSide][iFlx];
        
        // Loop over  variables
        for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
        {
          const CFreal avgSol = (m_gradTermL(iEq,iFlx)+m_gradTermR(iEq,iFlx))/2.0;
          ///@todo check if faceLocalDir is ok & faceFlxPntConn
	  if (cellSide == LEFT)
	  {
	    m_projectedCorrL = (avgSol-m_gradTermL(iEq,iFlx))*(m_faceJacobVecAbsSizeFlxPnts[iFlx]*(*m_faceMappedCoordDir)[orient][cellSide])*m_unitNormalFlxPnts[iFlx];
	  }
	  else
	  {
	    m_projectedCorrL = (avgSol-m_gradTermR(iEq,iFlx))*(m_faceJacobVecAbsSizeFlxPnts[iFlx]*(*m_faceMappedCoordDir)[orient][cellSide])*m_unitNormalFlxPnts[iFlx];
	  }

          // Loop over solution pnts to calculate the grad updates
          for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolDep; ++iSolPnt)
          {
            const CFuint iSolIdx = (*m_flxSolDep)[flxIdx][iSolPnt];

	    /// @todo Check if this is also OK for triangles!!
	    m_gradUpdates[cellSide][iSolIdx][iEq] += m_projectedCorrL*m_corrFctDiv[iSolIdx][flxIdx];
          }
        }
      }
      
//      for (CFuint jSolPnt = 0; jSolPnt < m_nbrSolSolDep; ++jSolPnt)
//  { 
//    const CFuint jSolIdx = (*m_solSolDep)[m_pertSol][jSolPnt];
//
//    // update gradients
//    (*m_cellGrads[side][jSolIdx])[m_pertVar] += m_gradUpdates[cellSide][jSolIdx][m_pertVar];
//  }

      // add the contribution to the gradients
      for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
      {
        for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
        {
          (*m_cellGrads[side][iSol])[iEq] += m_gradUpdates[cellSide][iSol][iEq];
          
//          if (m_cells[m_pertSide]->getID() == 5) 
//	  {
//            RealVector temp = m_gradUpdates[cellSide][iSol][iEq]/m_solJacobDet[side][iSol];
//            CFLog(INFO,"Num otherFace: " << iSol << ", "  << temp << "\n");
//          }
        }
      }
    }
  }
  
  // Add the contribution of the correction of the gradients for this face
  
  // compute face Jacobian vectors
  m_faceJacobVecs = m_face->computeFaceJacobDetVectorAtMappedCoords(*m_flxLocalCoords);
      
  // Loop over flux points to set the normal vectors
  for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
  {
    // get face Jacobian vector size
    m_faceJacobVecAbsSizeFlxPnts[iFlxPnt] = faceJacobVecSizeFaceFlxPnts[m_face->getID()][iFlxPnt];

    // set unit normal vector
    m_unitNormalFlxPnts[iFlxPnt] = m_faceJacobVecs[iFlxPnt]/m_faceJacobVecAbsSizeFlxPnts[iFlxPnt];
  }
    
  // loop over flx pnts to extrapolate the states to the flux points
  for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
  {     
    // local flux point indices in the left and right cell
    const CFuint flxPntIdxL = (*m_faceFlxPntConnPerOrient)[m_orient][LEFT][iFlxPnt];
    const CFuint flxPntIdxR = (*m_faceFlxPntConnPerOrient)[m_orient][RIGHT][iFlxPnt];
   
    *(m_cellStatesFlxPnt[LEFT][iFlxPnt]) = 0.0;
    *(m_cellStatesFlxPnt[RIGHT][iFlxPnt]) = 0.0;

    // extrapolate the left and right states to the flx pnts
    for (CFuint iSol = 0; iSol < m_nbrSolDep; ++iSol)
    {
      const CFuint solIdxL = (*m_flxSolDep)[flxPntIdxL][iSol];
      const CFuint solIdxR = (*m_flxSolDep)[flxPntIdxR][iSol];
          
      *(m_cellStatesFlxPnt[LEFT][iFlxPnt]) += (*m_solPolyValsAtFlxPnts)[flxPntIdxL][solIdxL]*(*(*m_states[LEFT])[solIdxL]);
      *(m_cellStatesFlxPnt[RIGHT][iFlxPnt]) += (*m_solPolyValsAtFlxPnts)[flxPntIdxR][solIdxR]*(*(*m_states[RIGHT])[solIdxR]);
    }
  }
      
  computeFaceGradTerms(m_gradTermL,m_gradTermR);
  
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
        
    // Loop over  variables
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      const CFreal avgSol = (m_gradTermL(iEq,iFlx)+m_gradTermR(iEq,iFlx))/2.0;
      ///@todo check if faceLocalDir is ok & faceFlxPntConn
      m_projectedCorrL = (avgSol-m_gradTermL(iEq,iFlx))*(m_faceJacobVecAbsSizeFlxPnts[iFlx]*(*m_faceMappedCoordDir)[m_orient][LEFT])*m_unitNormalFlxPnts[iFlx];
      m_projectedCorrR = (avgSol-m_gradTermR(iEq,iFlx))*(m_faceJacobVecAbsSizeFlxPnts[iFlx]*(*m_faceMappedCoordDir)[m_orient][RIGHT])*m_unitNormalFlxPnts[iFlx];

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

  // get jacobian determinants at solution points
  m_jacobDets[0] = m_cells[side]->computeGeometricShapeFunctionJacobianDeterminant(*m_solPntsLocalCoords);
  
//  for (CFuint jSolPnt = 0; jSolPnt < m_nbrSolSolDep; ++jSolPnt)
//  { 
//    const CFuint jSolIdx = (*m_solSolDep)[m_pertSol][jSolPnt];
//    
//    // inverse Jacobian determinant
//    const CFreal invJacobDet = 1.0/m_jacobDets[0][jSolIdx];
//
//    // update gradients
//    (*m_cellGrads[side][jSolIdx])[m_pertVar] += m_gradUpdates[side][jSolIdx][m_pertVar];
//    (*m_cellGrads[side][jSolIdx])[m_pertVar] *= invJacobDet;
//  }

  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {
    // inverse Jacobian determinant
    const CFreal invJacobDet = 1.0/m_jacobDets[0][iSol];

    // update gradients
    for (CFuint iGrad = 0; iGrad < m_nbrEqs; ++iGrad)
    {
      (*m_cellGrads[side][iSol])[iGrad] += m_gradUpdates[side][iSol][iGrad];
      (*m_cellGrads[side][iSol])[iGrad] *= invJacobDet;
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void DiffRHSJacobFluxReconstruction::computePerturbedGradientsAnalytical(const CFuint side)
{ 
  // Add the discontinuous gradient
  *m_cellStates = *(m_states[side]);

  computeCellGradTerm(m_gradTerm);

  for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
  {
    m_eps[iEq] = m_gradTerm(iEq,m_pertSol) - m_gradTermBefore(iEq,m_pertSol);
  }
  
  // Loop over solution pnts to calculate the grad updates
  for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolSolDep; ++iSolPnt)
  {
    const CFuint iSolIdx = (*m_solSolDep)[m_pertSol][iSolPnt];
    
    // inverse Jacobian determinant
    const CFreal invJacobDet = 1.0/m_solJacobDet[side][iSolIdx];
    
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      // Loop over gradient directions
      for (CFuint iDir = 0; iDir < m_dim; ++iDir)
      {
        m_projectedCorrL = m_eps[iEq] * m_neighbCellFluxProjVects[m_pertSide][iDir][m_pertSol];
	  
        // compute the grad updates
        (*m_cellGrads[side][iSolIdx])[iEq] += (*m_solPolyDerivAtSolPnts)[iSolIdx][iDir][m_pertSol]*m_projectedCorrL*invJacobDet;
      }
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

    if ((*m_isFaceOnBoundary[side])[faceIdx])
    {  
      // compute face Jacobian vectors
      m_faceJacobVecs = (*m_faces[side])[faceIdx]->computeFaceJacobDetVectorAtMappedCoords(*m_flxLocalCoords);
        
      // Loop over flux points to set the normal vectors
      for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
      {
        const CFuint currFlxIdx = (*m_faceFlxPntConn)[faceIdx][iFlxPnt];
      
        for (CFuint jFlxPnt = 0; jFlxPnt < m_nbrFlxDep; ++jFlxPnt)
        {
          if ((*m_solFlxDep)[m_pertSol][jFlxPnt] == currFlxIdx)
          {
            pertFlxPnt = iFlxPnt;
            pertFlxPntIdx = currFlxIdx;
            
            break;
          }
        }
        
        // get face Jacobian vector size
        m_faceJacobVecAbsSizeFlxPnts[iFlxPnt] = faceJacobVecSizeFaceFlxPnts[(*m_faces[side])[faceIdx]->getID()][iFlxPnt];

        // set unit normal vector
        m_unitNormalFlxPnts[iFlxPnt] = m_faceJacobVecs[iFlxPnt]/m_faceJacobVecAbsSizeFlxPnts[iFlxPnt];
        
        m_flxPntCoords[iFlxPnt] = (*m_faces[side])[faceIdx]->computeCoordFromMappedCoord((*m_flxLocalCoords)[iFlxPnt]);
        
        *(m_cellStatesFlxPnt[0][iFlxPnt]) = 0.0;
        
        for (CFuint iSol = 0; iSol < m_nbrSolDep; ++iSol)
        {
          const CFuint solIdx = (*m_flxSolDep)[currFlxIdx][iSol];
            
          *(m_cellStatesFlxPnt[0][iFlxPnt]) += (*m_solPolyValsAtFlxPnts)[currFlxIdx][solIdx]*(*((*m_states[side])[solIdx]));
        }
      }
      
      // compute ghost states with pert
      (*m_bcStateComputers)[(*m_faceBCIdx[side])[faceIdx]]->computeGhostStates(m_cellStatesFlxPnt[0],m_flxPntGhostSol,m_unitNormalFlxPnts,m_flxPntCoords);
      
      computeBndGradTerms(m_gradTermL,m_gradTermR);
      
      (*(m_cellStatesFlxPnt[0][pertFlxPnt]))[m_pertVar] -= m_numJacob->getEps() * (*m_solPolyValsAtFlxPnts)[pertFlxPntIdx][m_pertSol];
      
      // compute ghost states without pert
      (*m_bcStateComputers)[(*m_faceBCIdx[side])[faceIdx]]->computeGhostStates(m_cellStatesFlxPnt[0],m_flxPntGhostSol,m_unitNormalFlxPnts,m_flxPntCoords);
      
      computeBndGradTerms(m_gradTermL,m_gradTermTemp);
      
      for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
      {
        ///@todo check if faceLocalDir is ok & faceFlxPntConn
        m_projectedCorrL = 0.5*((m_gradTermR(iEq,pertFlxPnt)-m_gradTermTemp(iEq,pertFlxPnt))-m_eps[iEq]*(*m_solPolyValsAtFlxPnts)[pertFlxPntIdx][m_pertSol])*(m_faceJacobVecAbsSizeFlxPnts[pertFlxPnt]*(*m_faceLocalDir)[faceIdx])*m_unitNormalFlxPnts[pertFlxPnt];
        
        // Loop over solution pnts to calculate the grad updates
        for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolDep; ++iSolPnt)
        {
          const CFuint iSolIdx = (*m_flxSolDep)[pertFlxPntIdx][iSolPnt];
        
          // inverse Jacobian determinant
          const CFreal invJacobDet = 1.0/m_solJacobDet[side][iSolIdx];

          /// @todo Check if this is also OK for triangles!!
          (*m_cellGrads[side][iSolIdx])[iEq] += m_projectedCorrL*m_corrFctDiv[iSolIdx][pertFlxPntIdx]*invJacobDet;
//          if (m_cells[m_pertSide]->getID() == 1) 
//	  {
//          RealVector temp = m_projectedCorrL*m_corrFctDiv[iSolIdx][pertFlxPntIdx]*invJacobDet;
//            CFLog(INFO,"Ana Bnd: " << iSolIdx << ", "  << temp << "\n");
//          }
        }
      }
    }
    else
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
      
      for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
      {
        ///@todo check if faceLocalDir is ok & faceFlxPntConn
        m_projectedCorrL = -0.5*m_eps[iEq]*(*m_solPolyValsAtFlxPnts)[pertFlxPntIdx][m_pertSol]*(m_faceJacobVecAbsSizeFlxPnts[pertFlxPnt]*(*m_faceMappedCoordDir)[orient][cellSide])*m_unitNormalFlxPnts[pertFlxPnt];

        // Loop over solution pnts to calculate the grad updates
        for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolDep; ++iSolPnt)
        {
          const CFuint iSolIdx = (*m_flxSolDep)[pertFlxPntIdx][iSolPnt];
        
          // inverse Jacobian determinant
          const CFreal invJacobDet = 1.0/m_solJacobDet[side][iSolIdx];

          /// @todo Check if this is also OK for triangles!!
          (*m_cellGrads[side][iSolIdx])[iEq] += m_projectedCorrL*m_corrFctDiv[iSolIdx][pertFlxPntIdx]*invJacobDet;
//          if (m_cells[m_pertSide]->getID() == 5) 
//	  {
//            RealVector temp = m_projectedCorrL*m_corrFctDiv[iSolIdx][pertFlxPntIdx]*invJacobDet;
//              CFLog(INFO,"Ana otherFace: " << iSolIdx << ", "  << temp << "\n");
//          }
        }
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
      
  for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
  {
    ///@todo check if faceLocalDir is ok & faceFlxPntConn
    m_projectedCorrL = -0.5*m_eps[iEq]*(*m_solPolyValsAtFlxPnts)[pertFlxPntIdx][m_pertSol]*(m_faceJacobVecAbsSizeFlxPnts[pertFlxPnt]*(*m_faceMappedCoordDir)[m_orient][side])*m_unitNormalFlxPnts[pertFlxPnt];
    m_projectedCorrR = 0.5*m_eps[iEq]*(*m_solPolyValsAtFlxPnts)[pertFlxPntIdx][m_pertSol]*(m_faceJacobVecAbsSizeFlxPnts[pertFlxPnt]*(*m_faceMappedCoordDir)[m_orient][otherSide])*m_unitNormalFlxPnts[pertFlxPnt];

    // Loop over solution pnts to calculate the grad updates
    for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolDep; ++iSolPnt)
    {
      const CFuint iSolIdx = (*m_flxSolDep)[pertFlxPntIdx][iSolPnt];
      const CFuint iSolIdxOtherSide = (*m_flxSolDep)[pertFlxPntIdxOtherSide][iSolPnt];
    
      // inverse Jacobian determinant
      const CFreal invJacobDet = 1.0/m_solJacobDet[side][iSolIdx];
      const CFreal invJacobDetOtherSide = 1.0/m_solJacobDet[otherSide][iSolIdxOtherSide];

      /// @todo Check if this is also OK for triangles!!
      (*m_cellGrads[side][iSolIdx])[iEq] += m_projectedCorrL*m_corrFctDiv[iSolIdx][pertFlxPntIdx]*invJacobDet;
      (*m_cellGrads[otherSide][iSolIdxOtherSide])[iEq] += m_projectedCorrR*m_corrFctDiv[iSolIdxOtherSide][pertFlxPntIdxOtherSide]*invJacobDetOtherSide;
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void DiffRHSJacobFluxReconstruction::computePerturbedGradients()
{ 
  computeCellGradTerm(m_gradTerm);
  
  // Loop over solution pnts to calculate the grad updates
  for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolPnts; ++iSolPnt)
  {
    // Loop over  variables
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      //set the grad updates to 0 
      m_gradUpdates[0][iSolPnt][iEq] = 0.0;
    }
  }
  
  // Loop over solution pnts to calculate the grad updates
  for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolPnts; ++iSolPnt)
  {
    // Loop over  variables
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      // Loop over gradient directions
      for (CFuint iDir = 0; iDir < m_dim; ++iDir)
      {
        m_projectedCorrL = m_gradTerm(iEq,iSolPnt) * m_cellFluxProjVects[iDir][iSolPnt];
        
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
  
  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {
    // update gradients
    for (CFuint iGrad = 0; iGrad < m_nbrEqs; ++iGrad)
    {
      (*m_cellGrads[0][iSol])[iGrad] = m_gradUpdates[0][iSol][iGrad];
    }
  }
  
  // get face Jacobian vector sizes in the flux points
  DataHandle< vector< CFreal > > faceJacobVecSizeFaceFlxPnts = socket_faceJacobVecSizeFaceFlxPnts.getDataHandle();
  
  // Add the contribution of the correction to the gradients for each face
  // compute other face contributions to the gradients
  const CFuint nbrFaces = m_cell->nbNeighborGeos();
  for (CFuint iFace = 0; iFace < nbrFaces; ++iFace)
  {
    // get local face index
    const CFuint orient = (*m_faceOrientsCell)[iFace];
    
    // compute face Jacobian vectors
    m_faceJacobVecs = (*m_faces[0])[iFace]->computeFaceJacobDetVectorAtMappedCoords(*m_flxLocalCoords);

    // Loop over flux points to set the normal vectors
    for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
    {
      // get face Jacobian vector size
      m_faceJacobVecAbsSizeFlxPnts[iFlxPnt] = faceJacobVecSizeFaceFlxPnts[(*m_faces[0])[iFace]->getID()][iFlxPnt];

      // set unit normal vector
      m_unitNormalFlxPnts[iFlxPnt] = m_faceJacobVecs[iFlxPnt]/m_faceJacobVecAbsSizeFlxPnts[iFlxPnt];
      
      m_flxPntCoords[iFlxPnt] = (*m_faces[0])[iFace]->computeCoordFromMappedCoord((*m_flxLocalCoords)[iFlxPnt]);
    }

    if ((*m_isFaceOnBoundaryCell)[iFace])
    {
      // Loop over flux points to extrapolate the states to the flux points
      for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
      {
        *(m_cellStatesFlxPnt[0][iFlxPnt]) = 0.0;

        const CFuint currFlxIdx = (*m_faceFlxPntConn)[iFace][iFlxPnt];
        for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
        {
          *(m_cellStatesFlxPnt[0][iFlxPnt]) += (*m_solPolyValsAtFlxPnts)[currFlxIdx][iSol]*(*((*m_cellStates)[iSol]));
        }
      }
  
      // compute ghost states
      (*m_bcStateComputers)[(*m_faceBCIdxCell)[iFace]]->computeGhostStates(m_cellStatesFlxPnt[0],m_flxPntGhostSol,m_unitNormalFlxPnts,m_flxPntCoords);
      
      computeBndGradTerms(m_gradTermL,m_gradTermR);
      
      // Loop over solution pnts to reset the grad updates
      for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolPnts; ++iSolPnt)
      {
        // Loop over  variables
        for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
        {
          //set the grad updates to 0 
          m_gradUpdates[0][iSolPnt][iEq] = 0.0;
        }
      }
      
      // compute the face corrections to the gradients
      for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
      {
        const CFuint flxIdx = (*m_faceFlxPntConn)[iFace][iFlx];
        
        // Loop over  variables
        for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
        {
          const CFreal avgSol = (m_gradTermL(iEq,iFlx)+m_gradTermR(iEq,iFlx))/2.0;
          ///@todo check if faceLocalDir is ok & faceFlxPntConn
	  m_projectedCorrL = (avgSol-m_gradTermL(iEq,iFlx))*(m_faceJacobVecAbsSizeFlxPnts[iFlx]*(*m_faceLocalDir)[iFace])*m_unitNormalFlxPnts[iFlx];

          // Loop over solution pnts to calculate the grad updates
          for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolDep; ++iSolPnt)
          {
            const CFuint iSolIdx = (*m_flxSolDep)[flxIdx][iSolPnt];

	    /// @todo Check if this is also OK for triangles!!
	    m_gradUpdates[0][iSolIdx][iEq] += m_projectedCorrL*m_corrFctDiv[iSolIdx][flxIdx];
          }
        }
      }

      // add the contribution to the gradients
      for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
      {
        for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
        {
          (*m_cellGrads[0][iSol])[iEq] += m_gradUpdates[0][iSol][iEq];
        }
      }
    }
    else
    {
      // cell side with respect to this face
      const CFuint cellSide = (*m_currCellSideCell)[iFace];
      
      // loop over flx pnts to extrapolate the states to the flux points
      for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
      {     
        // local flux point indices in the left and right cell
        const CFuint flxPntIdxL = (*m_faceFlxPntConnPerOrient)[orient][LEFT][iFlxPnt];
        const CFuint flxPntIdxR = (*m_faceFlxPntConnPerOrient)[orient][RIGHT][iFlxPnt];
   
        *(m_cellStatesFlxPnt[LEFT][iFlxPnt]) = 0.0;
        *(m_cellStatesFlxPnt[RIGHT][iFlxPnt]) = 0.0;

        // extrapolate the left and right states to the flx pnts
        for (CFuint iSol = 0; iSol < m_nbrSolDep; ++iSol)
        {
          const CFuint solIdxL = (*m_flxSolDep)[flxPntIdxL][iSol];
          const CFuint solIdxR = (*m_flxSolDep)[flxPntIdxR][iSol];
          
          *(m_cellStatesFlxPnt[LEFT][iFlxPnt]) += (*m_solPolyValsAtFlxPnts)[flxPntIdxL][solIdxL]*(*(*m_faceNghbrStates[0][iFace][LEFT])[solIdxL]);
          *(m_cellStatesFlxPnt[RIGHT][iFlxPnt]) += (*m_solPolyValsAtFlxPnts)[flxPntIdxR][solIdxR]*(*(*m_faceNghbrStates[0][iFace][RIGHT])[solIdxR]);
        }
      }
      
      computeFaceGradTerms(m_gradTermL,m_gradTermR);
      
      // Loop over solution pnts to reset the grad updates
      for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolPnts; ++iSolPnt)
      {
        // Loop over  variables
        for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
        {
          //set the grad updates to 0 
          m_gradUpdates[cellSide][iSolPnt][iEq] = 0.0;
        }
      }
      
      // compute the face corrections to the gradients
      for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
      {
        const CFuint flxIdx = (*m_faceFlxPntConnPerOrient)[orient][cellSide][iFlx];
        
        // Loop over  variables
        for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
        {
          const CFreal avgSol = (m_gradTermL(iEq,iFlx)+m_gradTermR(iEq,iFlx))/2.0;
          ///@todo check if faceLocalDir is ok & faceFlxPntConn
          if (cellSide == LEFT)
	  {
	    m_projectedCorrL = (avgSol-m_gradTermL(iEq,iFlx))*(m_faceJacobVecAbsSizeFlxPnts[iFlx]*(*m_faceMappedCoordDir)[orient][cellSide])*m_unitNormalFlxPnts[iFlx];
	  }
	  else
	  {
	    m_projectedCorrL = (avgSol-m_gradTermR(iEq,iFlx))*(m_faceJacobVecAbsSizeFlxPnts[iFlx]*(*m_faceMappedCoordDir)[orient][cellSide])*m_unitNormalFlxPnts[iFlx];
	  }

          // Loop over solution pnts to calculate the grad updates
          for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolDep; ++iSolPnt)
          {
            const CFuint iSolIdx = (*m_flxSolDep)[flxIdx][iSolPnt];

	    /// @todo Check if this is also OK for triangles!!
	    m_gradUpdates[cellSide][iSolIdx][iEq] += m_projectedCorrL*m_corrFctDiv[iSolIdx][flxIdx];
          }
        }
      }

      // add the contribution to the gradients
      for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
      {
        for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
        {
          (*m_cellGrads[0][iSol])[iEq] += m_gradUpdates[cellSide][iSol][iEq];
        }
      }
    }
  }

  // get jacobian determinants at solution points
  m_jacobDets[0] = m_cell->computeGeometricShapeFunctionJacobianDeterminant(*m_solPntsLocalCoords);

  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {
    // inverse Jacobian determinant
    const CFreal invJacobDet = 1.0/m_jacobDets[0][iSol];

    // update gradients
    for (CFuint iGrad = 0; iGrad < m_nbrEqs; ++iGrad)
    {
      (*m_cellGrads[0][iSol])[iGrad] *= invJacobDet;
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void DiffRHSJacobFluxReconstruction::computePertCellDiffResiduals(const CFuint side)
{  
  // compute the volume term
  computeDivDiscontFlxNeighb(m_pertCellDiffRes,side);

  // add current face diffusive fluxes (m_pertResUpdates is set outside this function)
  m_pertCellDiffRes += m_pertResUpdates[side];
}

//////////////////////////////////////////////////////////////////////////////

void DiffRHSJacobFluxReconstruction::computeDivDiscontFlxNeighb(RealVector& residuals, const CFuint side)
{
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

    m_avgSol = *((*(m_states[side]))[iSolPnt]->getData());

    prepareFluxComputation();

    // calculate the discontinuous flux projected on x, y, z-directions
    for (CFuint iDim = 0; iDim < m_dim; ++iDim)
    {
//       m_contFlx[iSolPnt][iDim] = m_diffusiveVarSet->getFlux(m_avgSol,grad,m_cellFluxProjVects[iDim][iSolPnt],0);
       computeFlux(m_avgSol,m_tempGrad,m_neighbCellFluxProjVects[side][iDim][iSolPnt],0,m_contFlx[iSolPnt][iDim]);
    }

    for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFlxDep; ++iFlxPnt)
    {
      const CFuint flxIdx = (*m_solFlxDep)[iSolPnt][iFlxPnt];
      const CFuint dim = (*m_flxPntFlxDim)[flxIdx];

      m_extrapolatedFluxes[flxIdx] += (*m_solPolyValsAtFlxPnts)[flxIdx][iSolPnt]*(m_contFlx[iSolPnt][dim]);
    }
  }

  // reset the divergence of FC
  residuals = 0.0;
    
  // Loop over solution pnts to calculate the divergence of the discontinuous flux
  for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolPnts; ++iSolPnt)
  {
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
          residuals[m_nbrEqs*iSolPnt+iEq] += polyCoef*(m_contFlx[jSolIdx][iDir][iEq]);
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
        residuals[m_nbrEqs*iSolPnt+iVar] += -m_extrapolatedFluxes[flxIdx][iVar] * divh; 
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void DiffRHSJacobFluxReconstruction::computePertGradsFromFaceTerm(const CFuint side)
{
  // current face contribution to the gradients should have been computed in computePerturbedGradients()!

  // compute perturbed gradients
  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {
    const CFreal invJacobDet = 1.0/m_solJacobDet[side][iSol];
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      (*m_cellGrads[side][iSol])[iEq] =
          m_cellGradsMinusFaceTerm[side][iSol][iEq] +
          m_gradUpdates[side][iSol][iEq]*invJacobDet;
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void DiffRHSJacobFluxReconstruction::computePertGradsFromOtherFaceTerm(const CFuint side, const CFuint iFace)
{
  // get local face index
  const CFuint faceIdx = m_otherFaceLocalIdxs[side][iFace];
  cf_assert(!(*m_isFaceOnBoundary[side])[faceIdx]);
    
  // compute face Jacobian vectors
  m_faceJacobVecs = (*m_faces[side])[iFace]->computeFaceJacobDetVectorAtMappedCoords(*m_flxLocalCoords);
  
  // get face Jacobian vector sizes in the flux points
  DataHandle< vector< CFreal > > faceJacobVecSizeFaceFlxPnts = socket_faceJacobVecSizeFaceFlxPnts.getDataHandle();
  
  vector< RealVector > unitNormalFlxPnts;
  
  vector< vector< CFreal > > faceJacobVecSizeFlxPnts;
  faceJacobVecSizeFlxPnts.resize(m_nbrFaceFlxPnts);
  
  const CFuint orient = (*m_faceOrients[side])[faceIdx];
  
  // cell side with respect to this face
  const CFuint cellSide = (*m_currCellSide[side])[faceIdx];

  for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
  {
    // get face Jacobian vector size
    CFreal faceJacobVecAbsSizeFlxPnts = faceJacobVecSizeFaceFlxPnts[(*m_faces[side])[iFace]->getID()][iFlxPnt];

    // set unit normal vector
    unitNormalFlxPnts.push_back(m_faceJacobVecs[iFlxPnt]/faceJacobVecAbsSizeFlxPnts);
    
    faceJacobVecSizeFlxPnts[iFlxPnt].resize(2);
    // set face Jacobian vector size with sign depending on mapped coordinate direction
    faceJacobVecSizeFlxPnts[iFlxPnt][LEFT] = faceJacobVecAbsSizeFlxPnts*(*m_faceMappedCoordDir)[orient][LEFT];
    faceJacobVecSizeFlxPnts[iFlxPnt][RIGHT] = faceJacobVecAbsSizeFlxPnts*(*m_faceMappedCoordDir)[orient][RIGHT];
  }
    
  // loop over flx pnts to extrapolate the states to the flux points
  for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
  {     
    // local flux point indices in the left and right cell
    const CFuint flxPntIdxL = (*m_faceFlxPntConnPerOrient)[orient][LEFT][iFlxPnt];
    const CFuint flxPntIdxR = (*m_faceFlxPntConnPerOrient)[orient][RIGHT][iFlxPnt];
   
    *(m_cellStatesFlxPnt[LEFT][iFlxPnt]) = 0.0;
    *(m_cellStatesFlxPnt[RIGHT][iFlxPnt]) = 0.0;

    // extrapolate the left and right states to the flx pnts
    for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
    {
      *(m_cellStatesFlxPnt[LEFT][iFlxPnt]) += (*m_solPolyValsAtFlxPnts)[flxPntIdxL][iSol]*(*(*m_faceNghbrStates[side][iFace][LEFT])[iSol]);
      *(m_cellStatesFlxPnt[RIGHT][iFlxPnt]) += (*m_solPolyValsAtFlxPnts)[flxPntIdxR][iSol]*(*(*m_faceNghbrStates[side][iFace][RIGHT])[iSol]);
    }
  }
  
  computeFaceGradTerms(m_gradTermL,m_gradTermR);
      
  // compute the boundary face contribution to the gradients
  for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolPnts; ++iSolPnt)
  {
    // Loop over  variables
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      //set the grad updates to 0 
      m_gradUpdates[cellSide][iSolPnt][iEq] = 0.0;
      
      // compute the face corrections to the gradients
      for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
      {
        const CFreal avgSol = (m_gradTermL(iEq,iFlx)+m_gradTermR(iEq,iFlx))/2.0;

	if (cellSide == LEFT)
	{
          m_projectedCorrL = (avgSol-m_gradTermL(iEq,iFlx))*faceJacobVecSizeFlxPnts[iFlx][cellSide]*m_unitNormalFlxPnts[iFlx];
	}
	else
	{
	  m_projectedCorrL = (avgSol-m_gradTermR(iEq,iFlx))*faceJacobVecSizeFlxPnts[iFlx][cellSide]*m_unitNormalFlxPnts[iFlx];
	}
	/// @todo Check if this is also OK for triangles!!
	m_gradUpdates[cellSide][iSolPnt][iEq] += m_projectedCorrL*m_corrFctDiv[iSolPnt][(*m_faceFlxPntConnPerOrient)[orient][cellSide][iFlx]];
      }
    }
  }

  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {
    const CFreal invJacobDet = 1.0/m_solJacobDet[side][iSol];
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      (*m_cellGrads[side][iSol])[iEq] =
          m_cellGradsMinusOtherFaceTerm[side][iFace][iSol][iEq] +
          m_gradUpdates[cellSide][iSol][iEq]*invJacobDet;
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void DiffRHSJacobFluxReconstruction::setOtherFacesLocalIdxs()
{
  // get face ID of current face
  const CFuint currFaceID = m_face->getID();

  // loop over the faces of the cells
  for (CFuint iSide = 0; iSide < 2; ++iSide)
  {
    const CFuint nbrFaces = m_cells[iSide]->nbNeighborGeos();
    cf_assert(m_faces[iSide]->size() == nbrFaces);
    CFuint iFace = 0;
    for (CFuint faceIdx = 0; faceIdx < nbrFaces; ++faceIdx)
    {
      if ((*m_faces[iSide])[faceIdx]->getID() != currFaceID)
      {
        m_otherFaceLocalIdxs[iSide][iFace] = faceIdx;
        ++iFace;
      }
    }
    cf_assert(iFace == m_otherFaceLocalIdxs[iSide].size());
    cf_assert(iFace == nbrFaces-1);
  }
}

//////////////////////////////////////////////////////////////////////////////

void DiffRHSJacobFluxReconstruction::setFaceNeighbourStates()
{
  for (CFuint iSide = 0; iSide < 2; ++iSide)
  {
    // get neighbour states of other faces
    const CFuint nbrOtherFaces = m_otherFaceLocalIdxs[iSide].size();
    for (CFuint iFace = 0; iFace < nbrOtherFaces; ++iFace)
    {
      // get local face index
      const CFuint faceIdx = m_otherFaceLocalIdxs[iSide][iFace];

      // get neigbouring states
      const CFuint nbrFaceNeighbours = (*m_isFaceOnBoundary[iSide])[faceIdx] ? 1 : 2;
      for (CFuint iSide2 = 0; iSide2 < nbrFaceNeighbours; ++iSide2)
      {
        GeometricEntity* cell = (*m_faces[iSide])[faceIdx]->getNeighborGeo(iSide2);
        m_faceNghbrStates[iSide][iFace][iSide2] = cell->getStates();
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void DiffRHSJacobFluxReconstruction::setFaceNeighbourGradients()
{
  // get the gradients datahandle
  DataHandle< vector< RealVector > > gradients = socket_gradients.getDataHandle();

  for (CFuint iSide = 0; iSide < 2; ++iSide)
  {
    // get neighbour states of other faces
    const CFuint nbrOtherFaces = m_otherFaceLocalIdxs[iSide].size();
    for (CFuint iFace = 0; iFace < nbrOtherFaces; ++iFace)
    {
      // get local face index
      const CFuint faceIdx = m_otherFaceLocalIdxs[iSide][iFace];

      // get neigbouring states
      const CFuint nbrFaceNeighbours = (*m_isFaceOnBoundary[iSide])[faceIdx] ? 1 : 2;
      for (CFuint iSide2 = 0; iSide2 < nbrFaceNeighbours; ++iSide2)
      {
        // get number of states
        const CFuint nbrStates = m_faceNghbrStates[iSide][iFace][iSide2]->size();

        // resize m_faceNghbrGrads[iSide][iFace][iSide2]
        m_faceNghbrGrads[iSide][iFace][iSide2].resize(nbrStates);

        // set the gradients
        for (CFuint iState = 0; iState < nbrStates; ++iState)
        {
          const CFuint stateID = (*m_faceNghbrStates[iSide][iFace][iSide2])[iState]->getLocalID();
          m_faceNghbrGrads[iSide][iFace][iSide2][iState] = &gradients[stateID];
        }
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void DiffRHSJacobFluxReconstruction::computeCellGradsMinusFaceTerm()
{ 
  computeFaceGradTerms(m_gradTermL,m_gradTermR);
  
  // compute the face contribution to the gradients
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
    
    // Loop over  variables
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      const CFreal avgSol = (m_gradTermL(iEq,iFlx)+m_gradTermR(iEq,iFlx))/2.0;
      m_projectedCorrL = (avgSol-m_gradTermL(iEq,iFlx))*(m_faceJacobVecAbsSizeFlxPnts[iFlx]*(*m_faceMappedCoordDir)[m_orient][LEFT])*m_unitNormalFlxPnts[iFlx];
      m_projectedCorrR = (avgSol-m_gradTermR(iEq,iFlx))*(m_faceJacobVecAbsSizeFlxPnts[iFlx]*(*m_faceMappedCoordDir)[m_orient][RIGHT])*m_unitNormalFlxPnts[iFlx];
      
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

  // copy cell gradients and subtract face term
  for (CFuint iSide = 0; iSide < 2; ++iSide)
  {
    for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
    {
      const CFreal invJacobDet = 1.0/m_solJacobDet[iSide][iSol];
      for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
      {
        m_cellGradsMinusFaceTerm[iSide][iSol][iEq] = (*m_cellGrads[iSide][iSol])[iEq];
        m_cellGradsMinusFaceTerm[iSide][iSol][iEq] -= m_gradUpdates[iSide][iSol][iEq]*invJacobDet;
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void DiffRHSJacobFluxReconstruction::computeCellGradsMinusOtherFaceTerms(const CFuint side)
{
  const CFuint nbrOtherFaces = m_otherFaceLocalIdxs[side].size();
  for (CFuint iFace = 0; iFace < nbrOtherFaces; ++iFace)
  {
    // get local face index
    const CFuint faceIdx = m_otherFaceLocalIdxs[side][iFace];
    if(!(*m_isFaceOnBoundary[side])[faceIdx])
    {
      // compute face Jacobian vectors
      m_faceJacobVecs = (*m_faces[side])[iFace]->computeFaceJacobDetVectorAtMappedCoords(*m_flxLocalCoords);
  
      // get face Jacobian vector sizes in the flux points
      DataHandle< vector< CFreal > > faceJacobVecSizeFaceFlxPnts = socket_faceJacobVecSizeFaceFlxPnts.getDataHandle();
  
      const CFuint orient = (*m_faceOrients[side])[faceIdx];
  
      // cell side with respect to this face
      const CFuint cellSide = (*m_currCellSide[side])[faceIdx];

      for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
      {
        // get face Jacobian vector size
        m_faceJacobVecAbsSizeFlxPnts[iFlxPnt] = faceJacobVecSizeFaceFlxPnts[(*m_faces[side])[iFace]->getID()][iFlxPnt];

        // set unit normal vector
        m_unitNormalFlxPnts[iFlxPnt] = m_faceJacobVecs[iFlxPnt]/m_faceJacobVecAbsSizeFlxPnts[iFlxPnt];
      }
    
      // loop over flx pnts to extrapolate the states to the flux points
      for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
      {      
        // local flux point indices in the left and right cell
        const CFuint flxPntIdxL = (*m_faceFlxPntConnPerOrient)[orient][LEFT][iFlxPnt];
        const CFuint flxPntIdxR = (*m_faceFlxPntConnPerOrient)[orient][RIGHT][iFlxPnt];
   
        *(m_cellStatesFlxPnt[LEFT][iFlxPnt]) = 0.0;
        *(m_cellStatesFlxPnt[RIGHT][iFlxPnt]) = 0.0;

        // extrapolate the left and right states to the flx pnts
        for (CFuint iSol = 0; iSol < m_nbrSolDep; ++iSol)
        {
          const CFuint iSolIdxL = (*m_flxSolDep)[flxPntIdxL][iSol];
          const CFuint iSolIdxR = (*m_flxSolDep)[flxPntIdxR][iSol];
            
          *(m_cellStatesFlxPnt[LEFT][iFlxPnt]) += (*m_solPolyValsAtFlxPnts)[flxPntIdxL][iSolIdxL]*(*(*m_faceNghbrStates[side][iFace][LEFT])[iSolIdxL]);
          *(m_cellStatesFlxPnt[RIGHT][iFlxPnt]) += (*m_solPolyValsAtFlxPnts)[flxPntIdxR][iSolIdxR]*(*(*m_faceNghbrStates[side][iFace][RIGHT])[iSolIdxR]);
        }
      }
      
      computeFaceGradTerms(m_gradTermL,m_gradTermR);
      
      for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolPnts; ++iSolPnt)
      {
        for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
        {
          //set the grad updates to 0 
          m_gradUpdates[cellSide][iSolPnt][iEq] = 0.0;
        }
      }
      
      // compute the boundary face contribution to the gradients
      for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
      {
        const CFuint flxIdx = (*m_faceFlxPntConnPerOrient)[orient][cellSide][iFlx];
        
        // Loop over  variables
        for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
        {
          const CFreal avgSol = (m_gradTermL(iEq,iFlx)+m_gradTermR(iEq,iFlx))/2.0;
	  if (cellSide == LEFT)
	  {
            m_projectedCorrL = (avgSol-m_gradTermL(iEq,iFlx))*(m_faceJacobVecAbsSizeFlxPnts[iFlx]*(*m_faceMappedCoordDir)[orient][cellSide])*m_unitNormalFlxPnts[iFlx];
	  }
	  else
	  {
	    m_projectedCorrL = (avgSol-m_gradTermR(iEq,iFlx))*(m_faceJacobVecAbsSizeFlxPnts[iFlx]*(*m_faceMappedCoordDir)[orient][cellSide])*m_unitNormalFlxPnts[iFlx];
	  }
          for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolDep; ++iSolPnt)
          {
              const CFuint iSolIdx = (*m_flxSolDep)[flxIdx][iSolPnt];
              
            // compute the face corrections to the gradient
	    /// @todo Check if this is also OK for triangles!!
	    m_gradUpdates[cellSide][iSolIdx][iEq] += m_projectedCorrL*m_corrFctDiv[iSolIdx][flxIdx];
          }
        }
      }

      // copy cell gradients and subtract face term
      for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
      {
        const CFreal invJacobDet = 1.0/m_solJacobDet[side][iSol];

        for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
        {
          m_cellGradsMinusOtherFaceTerm[side][iFace][iSol][iEq] = (*m_cellGrads[side][iSol])[iEq];
          m_cellGradsMinusOtherFaceTerm[side][iFace][iSol][iEq] -= m_gradUpdates[cellSide][iSol][iEq]*invJacobDet;
        }
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void DiffRHSJacobFluxReconstruction::computeUnpertCellDiffResiduals()
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
    // create a list of the dimensions in which the deriv will be calculated
    for (CFuint iDim = 0; iDim < m_dim; ++iDim)
    {
      m_cellFluxProjVects[iDim] = m_cells[iSide]->computeMappedCoordPlaneNormalAtMappedCoords(m_dimList[iDim],*m_solPntsLocalCoords);
    }

    // set the states
    *m_cellStates = *(m_states[iSide]);

    // make a backup of the grads if necessary
    if (iSide == RIGHT)
    {
      for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
      {
        m_cellGradsBackUp[0][iSol] = *(m_cellGrads[0][iSol]);
        *(m_cellGrads[0][iSol]) = *(m_cellGrads[1][iSol]);
      }
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
        *(m_cellGrads[0][iSol]) = m_cellGradsBackUp[0][iSol];
      }
    }

    // current face term
    m_unpertCellDiffRes[iSide] += m_resUpdates[iSide];
//
//    // other face terms
//    const CFuint nbrOtherFaces = m_otherFaceLocalIdxs[iSide].size();
//    for (CFuint iFace = 0; iFace < nbrOtherFaces; ++iFace)
//    {
//      // get local face index
//      const CFuint faceIdx = m_otherFaceLocalIdxs[iSide][iFace];
//
//      if ((*m_isFaceOnBoundary[iSide])[faceIdx])
//      {
//        // compute the boundary face contribution to the diffusive residuals
//	computeBndRes(iSide, faceIdx, m_pertDivContFlx[0]);
//
//	// put the perturbed and unperturbed corrections in the correct format
//        // using m_pertResUpdates because the values stored in m_resUpdates should be preserved
//        for (CFuint iState = 0; iState < m_nbrSolPnts; ++iState)
//        {
//          for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
//          {
//            m_pertResUpdates[0][m_nbrEqs*iState+iVar] = m_pertDivContFlx[0][iState][iVar];
//          }
//        }
//
//        // add boundary face term
//        m_unpertCellDiffRes[iSide] += m_pertResUpdates[0];
//      }
//      else
//      {
//
//        // compute the internal face contribution to the diffusive residuals
//        // using m_pertResUpdates because the values stored in m_resUpdates should be preserved
//        // cell side with respect to this face
//        computeFaceRes(iSide, faceIdx, iFace, m_pertDivContFlx[0]);
//
//        // put the perturbed and unperturbed corrections in the correct format
//        for (CFuint iState = 0; iState < m_nbrSolPnts; ++iState)
//        {
//          for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
//          {
//            m_pertResUpdates[0][m_nbrEqs*iState+iVar] = m_pertDivContFlx[0][iState][iVar];
//          }
//        }
//
//        // add internal face term
//        m_unpertCellDiffRes[iSide] += m_pertResUpdates[0];
//      }
//    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void DiffRHSJacobFluxReconstruction::computeUnpertCellDiffResiduals(const CFuint side)
{
  // create a list of the dimensions in which the deriv will be calculated
  for (CFuint iDim = 0; iDim < m_dim; ++iDim)
  {
    m_cellFluxProjVects[iDim] = m_cells[side]->computeMappedCoordPlaneNormalAtMappedCoords(m_dimList[iDim],*m_solPntsLocalCoords);
  }

  // set the states
  *m_cellStates = *(m_states[side]);

  // make a backup of the grads if necessary
  if (side == RIGHT)
  {
    for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
    {
      m_cellGradsBackUp[0][iSol] = *(m_cellGrads[0][iSol]);
      *(m_cellGrads[0][iSol]) = *(m_cellGrads[1][iSol]);
    }
  }

  // compute the volume term
  computeDivDiscontFlx(m_pertDivContFlx[0]);

  // put the unpert discontinuous diff residual in the correct format
  for (CFuint iState = 0; iState < m_nbrSolPnts; ++iState)
  {
    for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
    {
      m_unpertCellDiffRes[side][m_nbrEqs*iState+iVar] = m_pertDivContFlx[0][iState][iVar];
    }
  }

  // restore grads if necessary
  if (side == RIGHT)
  {
    for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
    {
      *(m_cellGrads[0][iSol]) = m_cellGradsBackUp[0][iSol];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void DiffRHSJacobFluxReconstruction::computeBndRes(CFuint side, CFuint faceIdx, vector< RealVector >& corrections)
{
  vector< RealVector > unitNormalFlxPnts;
  
  vector< CFreal > faceJacobVecSizeFlxPnts;
  faceJacobVecSizeFlxPnts.resize(m_nbrFaceFlxPnts);
  
  // compute flux point coordinates
  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    m_flxPntCoords[iFlx] = (*m_faces[side])[faceIdx]->computeCoordFromMappedCoord((*m_flxLocalCoords)[iFlx]);	
  }
  
  // compute face Jacobian vectors
  m_faceJacobVecs = (*m_faces[side])[faceIdx]->computeFaceJacobDetVectorAtMappedCoords(*m_flxLocalCoords);
  
  // get face Jacobian vector sizes in the flux points
  DataHandle< vector< CFreal > >
  faceJacobVecSizeFaceFlxPnts = socket_faceJacobVecSizeFaceFlxPnts.getDataHandle();

  // Loop over flux points to extrapolate the states to the flux points and get the 
  // discontinuous normal flux in the flux points and the Riemann flux
  for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
  {
    // get face Jacobian vector size
    CFreal faceJacobVecAbsSizeFlxPnts = faceJacobVecSizeFaceFlxPnts[(*m_faces[side])[faceIdx]->getID()][iFlxPnt];

    // set face Jacobian vector size with sign depending on mapped coordinate direction
    faceJacobVecSizeFlxPnts[iFlxPnt] = faceJacobVecAbsSizeFlxPnts*((*m_faceLocalDir)[faceIdx]);

    // set unit normal vector
    unitNormalFlxPnts.push_back(m_faceJacobVecs[iFlxPnt]/faceJacobVecAbsSizeFlxPnts);
  }

  // Loop over flux points to extrapolate the states and gradients to the flux points
  for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
  {
    *(m_cellStatesFlxPnt[0][iFlxPnt]) = 0.0;
    
    for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
    {
      *(m_cellGradFlxPnt[0][iFlxPnt][iVar]) = 0.0;
    }
  
    const CFuint currFlxIdx = (*m_faceFlxPntConn)[faceIdx][iFlxPnt];
    for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
    {
      *(m_cellStatesFlxPnt[0][iFlxPnt]) += (*m_solPolyValsAtFlxPnts)[currFlxIdx][iSol]*(*((*m_states[side])[iSol]));
      
      for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
      {
        *(m_cellGradFlxPnt[0][iFlxPnt][iVar]) += (*m_solPolyValsAtFlxPnts)[currFlxIdx][iSol]*((*(m_cellGrads[side][iSol]))[iVar]);
      }
    }
  }

  // compute ghost states
  (*m_bcStateComputers)[(*m_faceBCIdx[side])[faceIdx]]->computeGhostStates(m_cellStatesFlxPnt[0],m_flxPntGhostSol,unitNormalFlxPnts,m_flxPntCoords);
   
  // compute ghost gradients
  (*m_bcStateComputers)[(*m_faceBCIdx[side])[faceIdx]]->computeGhostGradients(m_cellGradFlxPnt[0],m_flxPntGhostGrads,unitNormalFlxPnts,m_flxPntCoords);

  // compute the riemann flux minus the discontinuous flx in the flx pnts
  for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
  {
    //RealVector avgSol(m_nbrEqs);
    vector< RealVector* > avgGrad;
    avgGrad.resize(m_nbrEqs);

    for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
    {
      avgGrad[iVar] = new RealVector(m_dim);
      *(avgGrad[iVar]) = (*(m_flxPntGhostGrads[iFlxPnt][iVar]) + *(m_cellGradFlxPnt[0][iFlxPnt][iVar]))/2.0;

      m_avgSol[iVar] = ((*(m_flxPntGhostSol[iFlxPnt]))[iVar] + (*(m_cellStatesFlxPnt[0][iFlxPnt]))[iVar])/2.0; 
    }
    
    m_currFlx = iFlxPnt;
    
    prepareFluxComputation();

    computeFlux(m_avgSol,avgGrad,unitNormalFlxPnts[iFlxPnt],0,m_flxPntRiemannFlux[iFlxPnt]);
    
    m_cellFlx[0][iFlxPnt] = (m_flxPntRiemannFlux[iFlxPnt])*faceJacobVecSizeFlxPnts[iFlxPnt];
    
    for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
    {
      deletePtr(avgGrad[iVar]); 
    }
    avgGrad.clear();
  }

  cf_assert(corrections.size() == m_nbrSolPnts);
  
  for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolPnts; ++iSolPnt)
  {
    // reset the corrections
    corrections[iSolPnt] = 0.0;
    
    cf_assert(corrections[iSolPnt].size() == m_nbrEqs);

    // compute the term due to each flx pnt
    for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
    {
      // divergence of the correctionfct
      const CFreal divh = m_corrFctDiv[iSolPnt][(*m_faceFlxPntConn)[faceIdx][iFlxPnt]];
      
      if (divh != 0)
      {
        // the current correction factor
        RealVector currentCorrFactor = m_cellFlx[0][iFlxPnt];
        cf_assert(currentCorrFactor.size() == m_nbrEqs);
    
        // Fill in the corrections
        for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
        {
          corrections[iSolPnt][iVar] += currentCorrFactor[iVar] * divh;
        }
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void DiffRHSJacobFluxReconstruction::computeFaceRes(CFuint side, CFuint faceIdx, CFuint iFace, vector< RealVector >& corrections)
{
  const CFuint cellSide = (*m_currCellSide[side])[faceIdx];
  
  const CFuint orient = (*m_faceOrients[side])[faceIdx];
  
  vector< RealVector > unitNormalFlxPnts;
  
  vector< vector< CFreal > > faceJacobVecSizeFlxPnts;
  faceJacobVecSizeFlxPnts.resize(m_nbrFaceFlxPnts);
    
  // compute face Jacobian vectors
  m_faceJacobVecs = (*m_faces[side])[faceIdx]->computeFaceJacobDetVectorAtMappedCoords(*m_flxLocalCoords);

  // Loop over flux points to set the normal vectors
  for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
  {
    // get face Jacobian vector sizes in the flux points
    DataHandle< vector< CFreal > >
    faceJacobVecSizeFaceFlxPnts = socket_faceJacobVecSizeFaceFlxPnts.getDataHandle();
  
    // get face Jacobian vector size
    CFreal faceJacobVecAbsSizeFlxPnts = faceJacobVecSizeFaceFlxPnts[(*m_faces[side])[faceIdx]->getID()][iFlxPnt];

    faceJacobVecSizeFlxPnts[iFlxPnt].resize(2);
    // set face Jacobian vector size with sign depending on mapped coordinate direction
    faceJacobVecSizeFlxPnts[iFlxPnt][LEFT] = faceJacobVecAbsSizeFlxPnts*(*m_faceMappedCoordDir)[orient][LEFT];
    faceJacobVecSizeFlxPnts[iFlxPnt][RIGHT] = faceJacobVecAbsSizeFlxPnts*(*m_faceMappedCoordDir)[orient][RIGHT];

    // set unit normal vector
    unitNormalFlxPnts.push_back(m_faceJacobVecs[iFlxPnt]/faceJacobVecAbsSizeFlxPnts);
  }

  // loop over flx pnts to extrapolate the states to the flux points and get the 
  // discontinuous normal flux in the flux points and the Riemann flux
  for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
  {     
    // local flux point indices in the left and right cell
    const CFuint flxPntIdxL = (*m_faceFlxPntConnPerOrient)[orient][LEFT][iFlxPnt];
    const CFuint flxPntIdxR = (*m_faceFlxPntConnPerOrient)[orient][RIGHT][iFlxPnt]; 
    
    *(m_cellStatesFlxPnt[LEFT][iFlxPnt]) = 0.0;
    *(m_cellStatesFlxPnt[RIGHT][iFlxPnt]) = 0.0;
    
    for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
    {
      *(m_cellGradFlxPnt[LEFT][iFlxPnt][iVar]) = 0.0;
      *(m_cellGradFlxPnt[RIGHT][iFlxPnt][iVar]) = 0.0;
    }

    for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
    {
      *(m_cellStatesFlxPnt[LEFT][iFlxPnt]) += (*m_solPolyValsAtFlxPnts)[flxPntIdxL][iSol]*(*((*(m_faceNghbrStates[side][iFace][LEFT]))[iSol]));
      *(m_cellStatesFlxPnt[RIGHT][iFlxPnt]) += (*m_solPolyValsAtFlxPnts)[flxPntIdxR][iSol]*(*((*(m_faceNghbrStates[side][iFace][RIGHT]))[iSol]));
      
      for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
      {
        *(m_cellGradFlxPnt[LEFT][iFlxPnt][iVar]) += (*m_solPolyValsAtFlxPnts)[flxPntIdxL][iSol]*((*(m_faceNghbrGrads[side][iFace][LEFT][iSol]))[iVar]);
	*(m_cellGradFlxPnt[RIGHT][iFlxPnt][iVar]) += (*m_solPolyValsAtFlxPnts)[flxPntIdxR][iSol]*((*(m_faceNghbrGrads[side][iFace][RIGHT][iSol]))[iVar]);
      }
    }
  }

  // Loop over the flux points to calculate FI-FD
  for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
  {
    //RealVector avgSol(m_nbrEqs);
    vector< RealVector* > avgGrad;
    avgGrad.resize(m_nbrEqs);
     
    for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
    {
      avgGrad[iVar] = new RealVector(m_dim);
      *(avgGrad[iVar]) = (*(m_cellGradFlxPnt[LEFT][iFlxPnt][iVar]) + *(m_cellGradFlxPnt[RIGHT][iFlxPnt][iVar]))/2.0;
       
      m_avgSol[iVar] = ((*(m_cellStatesFlxPnt[LEFT][iFlxPnt]))[iVar] + (*(m_cellStatesFlxPnt[RIGHT][iFlxPnt]))[iVar])/2.0; 
    }

    m_currFlx = iFlxPnt;

    prepareFluxComputation();

    computeFlux(m_avgSol,avgGrad,unitNormalFlxPnts[iFlxPnt],0,m_flxPntRiemannFlux[iFlxPnt]);

    // compute FI-FD in the mapped coord frame
    m_cellFlx[cellSide][iFlxPnt] = (m_flxPntRiemannFlux[iFlxPnt])*faceJacobVecSizeFlxPnts[iFlxPnt][cellSide];

    for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
    {
      deletePtr(avgGrad[iVar]); 
    }
    avgGrad.clear();
  }

  // loop over sol pnts to compute the corrections
  for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolPnts; ++iSolPnt)
  {
    // reset the corrections which will be stored in divContFlx in order to be able to reuse updateRHS() 
    corrections[iSolPnt] = 0.0;
    
    cf_assert(corrections[iSolPnt].size() == m_nbrEqs);

    // compute the term due to each flx pnt
    for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
    {
      // divergence of the correction function
      const CFreal divh = m_corrFctDiv[iSolPnt][(*m_faceFlxPntConnPerOrient)[orient][cellSide][iFlxPnt]];

      if (divh != 0)
      {
        // the current correction factor (stored in cellFlx)
        const RealVector currentCorrFactor = m_cellFlx[cellSide][iFlxPnt];
        cf_assert(currentCorrFactor.size() == m_nbrEqs);
    
        // Fill in the corrections
        for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
        {
          corrections[iSolPnt][iVar] += currentCorrFactor[iVar] * divh;
        }
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void DiffRHSJacobFluxReconstruction::computeBndGradTerms(RealMatrix& gradTerm, RealMatrix& ghostGradTerm)
{ 
  for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
  {
    for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
    {
      gradTerm(iEq,iFlx) = (*(m_cellStatesFlxPnt[0][iFlx]->getData()))[iEq];
      ghostGradTerm(iEq,iFlx) = (*(m_flxPntGhostSol[iFlx]->getData()))[iEq];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void DiffRHSJacobFluxReconstruction::computeCellGradTerm(RealMatrix& gradTerm)
{
  for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
  {
    for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
    {
      gradTerm(iEq,iSol) = (*((*m_cellStates)[iSol]->getData()))[iEq];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void DiffRHSJacobFluxReconstruction::computeFaceGradTerms(RealMatrix& gradTermL, RealMatrix& gradTermR)
{
  for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
  {
    for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
    {
      gradTermL(iEq,iFlx) = (*(m_cellStatesFlxPnt[LEFT][iFlx]->getData()))[iEq];
      gradTermR(iEq,iFlx) = (*(m_cellStatesFlxPnt[RIGHT][iFlx]->getData()))[iEq];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void DiffRHSJacobFluxReconstruction::updateRHSUnpertCell(const CFuint side)
{
  // get the datahandle of the rhs
  DataHandle< CFreal > rhs = socket_rhs.getDataHandle();

  // get residual factor
  const CFreal resFactor = getMethodData().getResFactor();

  // update rhs
  for (CFuint iState = 0; iState < m_nbrSolPnts; ++iState)
  {
    CFuint resID = m_nbrEqs*( (*m_states[side])[iState]->getLocalID() );
    for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
    {
      rhs[resID+iVar] += resFactor*m_unpertCellDiffRes[side][m_nbrEqs*iState+iVar];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void DiffRHSJacobFluxReconstruction::setup()
{
  CFAUTOTRACE;
  CFLog(VERBOSE, "DiffJacob setup\n");
  // setup parent class
  DiffRHSFluxReconstruction::setup();
  
  // get CellToFaceGeBuilders
  m_cellBuilders.resize(2);
  m_cellBuilders[LEFT ] = getMethodData().getCellBuilder();
  m_cellBuilders[RIGHT] = getMethodData().getSecondCellBuilder();

  // get some additional data for cell building
  m_isFaceOnBoundary.resize(2);
  m_nghbrCellSide   .resize(2);
  m_currCellSide    .resize(2);
  m_faceOrients     .resize(2);
  m_faceBCIdx       .resize(2);
  for (CFuint iSide = 0; iSide < 2; ++iSide)
  {
    m_isFaceOnBoundary[iSide] = m_cellBuilders[iSide]->getGeoBuilder()->getIsFaceOnBoundary();
    m_nghbrCellSide   [iSide] = m_cellBuilders[iSide]->getGeoBuilder()->getNeighbrCellSide ();
    m_currCellSide    [iSide] = m_cellBuilders[iSide]->getGeoBuilder()->getCurrentCellSide ();
    m_faceOrients     [iSide] = m_cellBuilders[iSide]->getGeoBuilder()->getFaceOrient      ();
    m_faceBCIdx       [iSide] = m_cellBuilders[iSide]->getGeoBuilder()->getFaceBCIdx       ();
  }

  
  // get CellToFaceGeBuilder
  m_cellBuilder      = getMethodData().getCellBuilder();
  m_isFaceOnBoundaryCell = m_cellBuilder->getGeoBuilder()->getIsFaceOnBoundary();
  m_nghbrCellSideCell    = m_cellBuilder->getGeoBuilder()->getNeighbrCellSide ();
  m_currCellSideCell     = m_cellBuilder->getGeoBuilder()->getCurrentCellSide ();
  m_faceOrientsCell      = m_cellBuilder->getGeoBuilder()->getFaceOrient      ();
  m_faceBCIdxCell        = m_cellBuilder->getGeoBuilder()->getFaceBCIdx       ();
  
  // get BC state computers
  m_bcStateComputers = getMethodData().getBCStateComputers();

  // get the linear system solver
  m_lss = getMethodData().getLinearSystemSolver()[0];

  // get the numerical Jacobian computer
  m_numJacob = getMethodData().getNumericalJacobian();

  // get number of equations and dimensionality
  const CFuint dim    = PhysicalModelStack::getActive()->getDim();

  // get the local spectral FD data
  vector< FluxReconstructionElementData* >& frLocalData = getMethodData().getFRLocalData();
  cf_assert(frLocalData.size() > 0);

  // get the number of faces in a cell
  const CFuint nbrFaces = frLocalData[0]->getNbrCellFaces();
  const CFuint nbrFacesM1 = nbrFaces - 1;

  // create blockaccumulator
  m_acc.reset(m_lss->createBlockAccumulator(2*m_nbrSolPnts,2*m_nbrSolPnts,m_nbrEqs));

  // create single cell blockaccumulator
  m_accSC.reset(m_lss->createBlockAccumulator(m_nbrSolPnts,m_nbrSolPnts,m_nbrEqs));

  // resize m_faces
  m_faces.resize(2);

  // resize m_otherFaceLocalIdxs
  m_otherFaceLocalIdxs.resize(2);
  m_otherFaceLocalIdxs[LEFT ].resize(nbrFacesM1);
  m_otherFaceLocalIdxs[RIGHT].resize(nbrFacesM1);

  // resize m_faceNghbrStates
  m_faceNghbrStates.resize(2);
  for (CFuint iSide = 0; iSide < 2; ++iSide)
  {
    m_faceNghbrStates[iSide].resize(nbrFaces);
    for (CFuint iFace = 0; iFace < nbrFaces; ++iFace)
    {
      m_faceNghbrStates[iSide][iFace].resize(2);
    }
  }

  // resize m_faceNghbrGrads
  m_faceNghbrGrads.resize(2);
  for (CFuint iSide = 0; iSide < 2; ++iSide)
  {
    m_faceNghbrGrads[iSide].resize(nbrFacesM1);
    for (CFuint iFace = 0; iFace < nbrFacesM1; ++iFace)
    {
      m_faceNghbrGrads[iSide][iFace].resize(2);
    }
  }

  // resize variables
  const CFuint nbrCellResiduals = m_nbrSolPnts*m_nbrEqs;
  m_pertResUpdates.resize(2);
  m_pertResUpdates[LEFT ].resize(nbrCellResiduals);
  m_pertResUpdates[RIGHT].resize(nbrCellResiduals);
  m_derivResUpdates.resize(nbrCellResiduals);
  m_unpertCellDiffRes       .resize(2);
  m_unpertCellDiffRes[LEFT ].resize(nbrCellResiduals);
  m_unpertCellDiffRes[RIGHT].resize(nbrCellResiduals);
  m_pertCellDiffRes         .resize(nbrCellResiduals);
  m_derivCellDiffRes        .resize(nbrCellResiduals);
  m_divContFlxL.resize(m_nbrSolPnts);
  m_divContFlxR.resize(m_nbrSolPnts);
  m_resUpdates .resize(2);
  m_resUpdates [LEFT ].resize(nbrCellResiduals);
  m_resUpdates [RIGHT].resize(nbrCellResiduals);
  
  m_cellGradsBackUp.resize(2);
  m_cellGradsBackUp[LEFT].resize(m_nbrSolPnts);
  m_cellGradsBackUp[RIGHT].resize(m_nbrSolPnts);
  m_pertDivContFlx .resize(2);
  m_pertDivContFlx [LEFT ].resize(m_nbrSolPnts);
  m_pertDivContFlx [RIGHT].resize(m_nbrSolPnts);
  m_pertCorrections.resize(m_nbrSolPnts);
  m_cellGradFlxPntBackup.resize(m_nbrFaceFlxPnts);
  
  m_dimList.resize(m_dim);

  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {
    m_cellGradsBackUp[LEFT][iSol].resize(m_nbrEqs);
    m_cellGradsBackUp[RIGHT][iSol].resize(m_nbrEqs);
    m_pertDivContFlx [LEFT ][iSol].resize(m_nbrEqs);
    m_pertDivContFlx [RIGHT][iSol].resize(m_nbrEqs);
    m_pertCorrections[iSol].resize(m_nbrEqs);
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      m_cellGradsBackUp[LEFT][iSol][iEq].resize(m_dim);
      m_cellGradsBackUp[RIGHT][iSol][iEq].resize(m_dim);
    }
  }
  
  for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolPnts; ++iSolPnt)
  {
    m_divContFlxL[iSolPnt].resize(m_nbrEqs);
    m_divContFlxR[iSolPnt].resize(m_nbrEqs);
  }

  // allocate memory for perturbed gradients
  m_pertGrads.resize(2);
  for (CFuint iSide = 0; iSide < 2; ++iSide)
  {
    m_pertGrads[iSide].resize(m_nbrSolPnts);
    for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
    {
      m_pertGrads[iSide][iSol] = new vector<RealVector>(m_nbrEqs);
      for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
      {
        (*m_pertGrads[iSide][iSol])[iEq].resize(dim);
      }
    }
  }

  // resize m_cellGradsMinusFaceTerm
  m_cellGradsMinusFaceTerm.resize(2);
  for (CFuint iSide = 0; iSide < 2; ++iSide)
  {
    m_cellGradsMinusFaceTerm[iSide].resize(m_nbrSolPnts);
    for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
    {
      m_cellGradsMinusFaceTerm[iSide][iSol].resize(m_nbrEqs);
      for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
      {
        m_cellGradsMinusFaceTerm[iSide][iSol][iEq].resize(dim);
      }
    }
  }

  // resize m_cellGradsMinusOtherFaceTerm
  m_cellGradsMinusOtherFaceTerm.resize(2);
  for (CFuint iSide = 0; iSide < 2; ++iSide)
  {
    m_cellGradsMinusOtherFaceTerm[iSide].resize(nbrFacesM1);
    for (CFuint iFace = 0; iFace < nbrFacesM1; ++iFace)
    {
      m_cellGradsMinusOtherFaceTerm[iSide][iFace].resize(m_nbrSolPnts);
      for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
      {
        m_cellGradsMinusOtherFaceTerm[iSide][iFace][iSol].resize(m_nbrEqs);
        for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
        {
          m_cellGradsMinusOtherFaceTerm[iSide][iFace][iSol][iEq].resize(dim);
        }
      }
    }
  }

  // resize gradient updates
  m_gradUpdates.resize(2);
  for (CFuint iSide = 0; iSide < 2; ++iSide)
  {
    m_gradUpdates[iSide].resize(m_nbrSolPnts);
    for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
    {
      m_gradUpdates[iSide][iSol].resize(m_nbrEqs);
      for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
      {
        m_gradUpdates[iSide][iSol][iEq].resize(dim);
      }
    }
  }

  // resize neighbouring cells solution points Jacobian determinants
  m_solJacobDet.resize(2);
  for (CFuint iSide = 0; iSide < 2; ++iSide)
  {
    m_solJacobDet[iSide].resize(m_nbrSolPnts);
  }
  
  m_flxPntGhostGrads.resize(m_nbrFaceFlxPnts);
  
  // create internal and ghost states
  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    m_flxPntGhostSol.push_back(new State());
    for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
    {
      m_flxPntGhostGrads[iFlx].push_back(new RealVector(m_dim));
    }
  }

  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    m_flxPntGhostSol[iFlx]->setLocalID(iFlx);
    
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      m_cellGradFlxPntBackup[iFlx].push_back(new RealVector(m_dim));
    }
  }
  
  for (CFuint iDim = 0; iDim < m_dim; ++iDim)
  {
    m_dimList[iDim].resize(m_nbrSolPnts);
    for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolPnts; ++iSolPnt)
    {
      m_dimList[iDim][iSolPnt] = iDim;
    }
  }
  
  m_gradTermL.resize(m_nbrEqs,m_nbrFaceFlxPnts);
  m_gradTermR.resize(m_nbrEqs,m_nbrFaceFlxPnts);
  m_gradTermTemp.resize(m_nbrEqs,m_nbrFaceFlxPnts);
  m_gradTerm.resize(m_nbrEqs,m_nbrSolPnts);
  m_gradTermBefore.resize(m_nbrEqs,m_nbrSolPnts);
  m_projectedCorrL.resize(m_dim);
  m_projectedCorrR.resize(m_dim);
  m_eps.resize(m_nbrEqs);
  m_neighbCellFluxProjVects.resize(2);
  m_neighbCellFluxProjVects[LEFT].resize(m_dim);
  m_neighbCellFluxProjVects[RIGHT].resize(m_dim);
  
  // get the elementTypeData
  SafePtr< vector<ElementTypeData> > elemType = MeshDataStack::getActive()->getElementTypeData();
  
  // get the number of cells in the mesh
  const CFuint nbrCells = (*elemType)[0].getEndIdx();
  
  m_cellFlags.resize(nbrCells);
  m_unpertAllCellDiffRes.resize(nbrCells);
  
  for (CFuint iCell = 0; iCell < nbrCells; ++iCell)
  {
    m_unpertAllCellDiffRes[iCell].resize(nbrCellResiduals);
  }
  
  for (CFuint iDim = 0; iDim < m_dim; ++iDim)
  {
    m_neighbCellFluxProjVects[iDim].resize(m_nbrSolPnts);
    for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolPnts; ++iSolPnt)
    {
      m_neighbCellFluxProjVects[iDim][iSolPnt].resize(m_dim);
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void DiffRHSJacobFluxReconstruction::unsetup()
{
  CFAUTOTRACE;
  
  for (CFuint iSide = 0; iSide < m_pertGrads.size(); ++iSide)
  {
    for (CFuint iSol = 0; iSol < m_pertGrads[iSide].size(); ++iSol)
    {
      deletePtr(m_pertGrads[iSide][iSol]);
    }
    m_pertGrads[iSide].resize(0);
  }
  m_pertGrads.resize(0);
  
  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    deletePtr(m_flxPntGhostSol[iFlx]);
    for (CFuint iGrad = 0; iGrad < m_nbrEqs; ++iGrad)
    {
      deletePtr(m_flxPntGhostGrads[iFlx][iGrad]); 
    }
    m_flxPntGhostGrads[iFlx].clear();
  }
  m_flxPntGhostSol.clear();
  m_flxPntGhostGrads.clear();
  
  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    for (CFuint iGrad = 0; iGrad < m_nbrEqs; ++iGrad)
    {
      deletePtr(m_cellGradFlxPntBackup[iFlx][iGrad]); 
    }
    m_cellGradFlxPntBackup[iFlx].clear();
  }
  
  // unsetup parent class
  DiffRHSFluxReconstruction::unsetup();
}


//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD

