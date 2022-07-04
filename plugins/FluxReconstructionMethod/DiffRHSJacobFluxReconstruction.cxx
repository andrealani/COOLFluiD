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
  m_neighbCellFluxProjVects(),
  m_affectedSolPnts(),
  m_contFlxBackup(),
  m_contFlxNeighb()
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
      }
      
      // release the GeometricEntity
      m_faceBuilder->releaseGE();
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void DiffRHSJacobFluxReconstruction::computeBothJacobsDiffFaceTerm()
{
  // get residual factor
  const CFreal resFactor = getMethodData().getResFactor();

  // dereference accumulator
  BlockAccumulator& acc = *m_acc;

  // set block row and column indices, proj vectors and make a backup of discontinuous fluxes
  CFuint solIdx = 0;
  for (m_pertSide = 0; m_pertSide < 2; ++m_pertSide)
  {
    for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol, ++solIdx)
    {
      acc.setRowColIndex(solIdx,(*m_states[m_pertSide])[iSol]->getLocalID());
    }
    
    for (CFuint iDim = 0; iDim < m_dim+m_ndimplus; ++iDim)
    {
      m_neighbCellFluxProjVects[m_pertSide][iDim] = m_cells[m_pertSide]->computeMappedCoordPlaneNormalAtMappedCoords(m_dimList[iDim],*m_solPntsLocalCoords);
    }
    
    // Loop over solution points to calculate the discontinuous flux.
    for (m_pertSol = 0; m_pertSol < m_nbrSolPnts; ++m_pertSol)
    { 
      for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
      {
        *(m_tempGrad[iVar]) = (*(m_cellGrads[m_pertSide][m_pertSol]))[iVar];
      }

      m_avgSol = *((*(m_states[m_pertSide]))[m_pertSol]->getData());

      prepareFluxComputation();

      // calculate the discontinuous flux projected on x, y, z-directions
      for (CFuint iDim = 0; iDim < m_dim+m_ndimplus; ++iDim)
      {
        computeFlux(m_avgSol,m_tempGrad,m_neighbCellFluxProjVects[m_pertSide][iDim][m_pertSol],0,m_contFlxNeighb[m_pertSide][m_pertSol][iDim]);
        m_contFlxBackup[m_pertSide][m_pertSol][iDim] = m_contFlxNeighb[m_pertSide][m_pertSol][iDim];
      }
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
    
    // Add the discontinuous gradient
    *m_cellStates = *(m_states[m_pertSide]);
  
    computeCellGradTerm(m_gradTermBefore);

    // loop over the states to perturb the states
    for (m_pertSol = 0; m_pertSol < m_nbrSolPnts; ++m_pertSol)
    {
      // dereference state
      State& pertState = *(*m_states[m_pertSide])[m_pertSol];

      // reset affected sol pnts
      for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
      {
        m_affectedSolPnts[LEFT][iSol] = false;
        m_affectedSolPnts[RIGHT][iSol] = false;
      }

      // loop over the variables in the state
      for (m_pertVar = 0; m_pertVar < m_nbrEqs; ++m_pertVar)
      {
        // perturb physical variable in state
        m_numJacob->perturb(m_pertVar,pertState[m_pertVar]);

        // compute the perturbed gradients in the current cell
        computePerturbedGradientsAnalytical(m_pertSide); 
        
//            for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
//    {
//        RealVector temp(m_nbrEqs);
//        RealVector temp2(m_nbrEqs);
//        RealVector dq(m_nbrEqs);
//        for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
//      {
//temp[iEq] =m_cellGradsBackUp[LEFT][iSol][iEq][XX];
//        }
//        
//        for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
//      {
//temp2[iEq] = (*(m_cellGrads[LEFT][iSol]))[iEq][XX];
//        }
//        
// m_numJacob->computeDerivative(temp,temp2,dq);
// if(m_cells[LEFT]->getID()==1&&m_pertVar==3) CFLog(INFO,"side: " << m_pertSide << ", sol: " << m_pertSol << ", to side: 0, sol: " << iSol << ": " << dq[3]/0.0013935);
// 
// for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
//      {
//temp[iEq] =m_cellGradsBackUp[LEFT][iSol][iEq][YY];
//        }
//        
//        for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
//      {
//temp2[iEq] = (*(m_cellGrads[LEFT][iSol]))[iEq][YY];
//        }
//        
// m_numJacob->computeDerivative(temp,temp2,dq);
// if(m_cells[LEFT]->getID()==1&&m_pertVar==3) CFLog(INFO," " << dq[3]/0.0013935 << "\n");
// 
// for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
//      {
//temp[iEq] =m_cellGradsBackUp[RIGHT][iSol][iEq][XX];
//        }
//        
//        for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
//      {
//temp2[iEq] = (*(m_cellGrads[RIGHT][iSol]))[iEq][XX];
//        }
//        
// m_numJacob->computeDerivative(temp,temp2,dq);
// if(m_cells[LEFT]->getID()==1&&m_pertVar==3) CFLog(INFO,"side: " << m_pertSide << ", sol: " << m_pertSol << ", to side: 1, sol: " << iSol << ": " << dq[3]/0.0013935);
// 
// for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
//      {
//temp[iEq] =m_cellGradsBackUp[RIGHT][iSol][iEq][YY];
//        }
//        
//        for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
//      {
//temp2[iEq] = (*(m_cellGrads[RIGHT][iSol]))[iEq][YY];
//        }
        
// m_numJacob->computeDerivative(temp,temp2,dq);
// if(m_cells[LEFT]->getID()==1&&m_pertVar==3) CFLog(INFO," " << dq[3]/0.0013935 << "\n");
//        }

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
	if (m_cells[m_pertSide]->getID() == 94) 
	{
	  CFLog(VERBOSE, "pert1: " << m_pertResUpdates[m_pertSide] << "\n");
	  CFLog(VERBOSE, "unpert1: " << m_resUpdates[m_pertSide] << "\n");
	  CFLog(VERBOSE, "deriv1: " << m_derivResUpdates << "\n");
          CFLog(VERBOSE, "otherID: " << m_cells[iOtherSide]->getID() << "\n");
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
	//if (m_cells[0]->getID() == 1) CFLog(INFO, "pertSol: " << m_pertSol << ", pertVar: " << m_pertVar << ", Jllav: " << m_derivCellDiffRes << "\n");

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
//            if (m_cells[LEFT]->getID() == 1 && iSol+pertSideTerm==0 && m_pertSol+pertSideTerm==0 && m_pertVar==3) CFLog(INFO, "adding: " << m_derivCellDiffRes << "\n"); 
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
          
          for (CFuint iSide = 0; iSide < 2; ++iSide)
          {
            if (m_affectedSolPnts[iSide][iState])
            {
              // calculate the discontinuous flux projected on x, y, z-directions
              for (CFuint iDim = 0; iDim < m_dim+m_ndimplus; ++iDim)
              {
                m_contFlxNeighb[iSide][iState][iDim] = m_contFlxBackup[iSide][iState][iDim];
              }
            }
          }
        }
      }
    }
  }
  
  if (m_cells[LEFT]->getID() == 1) 
  {
      //CFLog(INFO, "ACCDiff: " << acc.getValue(0,4,3,3) << "\n");
      //acc.printToScreen();
  }
  if (m_cells[RIGHT]->getID() == 1) 
  {
      //CFLog(INFO, "ACCDiff: " << acc.getValue(4,4,3,3) << "\n");
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

void DiffRHSJacobFluxReconstruction::computeOneJacobDiffFaceTerm(const CFuint side)
{
  // get residual factor
  const CFreal resFactor = getMethodData().getResFactor();

  // dereference accumulator
  BlockAccumulator& acc = *m_acc;

  // set block row and column indices
  CFuint solIdx = 0;
  for (m_pertSide = 0; m_pertSide < 2; ++m_pertSide)
  {
    for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol, ++solIdx)
    {
      acc.setRowColIndex(solIdx,(*m_states[m_pertSide])[iSol]->getLocalID());
    }
    
    for (CFuint iDim = 0; iDim < m_dim+m_ndimplus; ++iDim)
    {
      m_neighbCellFluxProjVects[m_pertSide][iDim] = m_cells[m_pertSide]->computeMappedCoordPlaneNormalAtMappedCoords(m_dimList[iDim],*m_solPntsLocalCoords);
    }
    
    // Loop over solution points to calculate the discontinuous flux.
    for (m_pertSol = 0; m_pertSol < m_nbrSolPnts; ++m_pertSol)
    { 
      for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
      {
        *(m_tempGrad[iVar]) = (*(m_cellGrads[m_pertSide][m_pertSol]))[iVar];
      }

      m_avgSol = *((*(m_states[m_pertSide]))[m_pertSol]->getData());

      prepareFluxComputation();

      // calculate the discontinuous flux projected on x, y, z-directions
      for (CFuint iDim = 0; iDim < m_dim+m_ndimplus; ++iDim)
      {
        computeFlux(m_avgSol,m_tempGrad,m_neighbCellFluxProjVects[m_pertSide][iDim][m_pertSol],0,m_contFlxNeighb[m_pertSide][m_pertSol][iDim]);
        m_contFlxBackup[m_pertSide][m_pertSol][iDim] = m_contFlxNeighb[m_pertSide][m_pertSol][iDim];
      }
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
    
    // Add the discontinuous gradient
    *m_cellStates = *(m_states[m_pertSide]);
  
    computeCellGradTerm(m_gradTermBefore);

    // loop over the states to perturb the states
    for (m_pertSol = 0; m_pertSol < m_nbrSolPnts; ++m_pertSol)
    {
      // dereference state
      State& pertState = *(*m_states[m_pertSide])[m_pertSol];

      // reset affected sol pnts
      for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
      {
        m_affectedSolPnts[LEFT][iSol] = false;
        m_affectedSolPnts[RIGHT][iSol] = false;
      }

      // loop over the variables in the state
      for (m_pertVar = 0; m_pertVar < m_nbrEqs; ++m_pertVar)
      {
        // perturb physical variable in state
        m_numJacob->perturb(m_pertVar,pertState[m_pertVar]);

        // compute the perturbed gradients in the current cell
        computePerturbedGradientsAnalytical(m_pertSide);        

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
          
          for (CFuint iSide = 0; iSide < 2; ++iSide)
          {
            if (m_affectedSolPnts[iSide][iState])
            {
              // calculate the discontinuous flux projected on x, y, z-directions
              for (CFuint iDim = 0; iDim < m_dim+m_ndimplus; ++iDim)
              {
                m_contFlxNeighb[iSide][iState][iDim] = m_contFlxBackup[iSide][iState][iDim];
              }
            }
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
    
    m_affectedSolPnts[side][iSolIdx] = true;
    
    // inverse Jacobian determinant
    const CFreal invJacobDet = 1.0/m_solJacobDet[side][iSolIdx];
    
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      // Loop over gradient directions
      for (CFuint iDir = 0; iDir < m_dim; ++iDir)
      {
        m_projectedCorrL = m_eps[iEq] * (m_neighbCellFluxProjVects[m_pertSide][iDir+m_ndimplus][m_pertSol]);
	  
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
      
      m_affectedSolPnts[otherSide][iSolIdxOtherSide] = true;
    
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
    if (m_affectedSolPnts[side][iSolPnt])
    {
      for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
      {
        *(m_tempGrad[iVar]) = (*(m_cellGrads[side][iSolPnt]))[iVar];
      }

      m_avgSol = *((*(m_states[side]))[iSolPnt]->getData());

      prepareFluxComputation();

      // calculate the discontinuous flux projected on x, y, z-directions
      for (CFuint iDim = 0; iDim < m_dim+m_ndimplus; ++iDim)
      {
//       m_contFlx[iSolPnt][iDim] = m_diffusiveVarSet->getFlux(m_avgSol,grad,m_cellFluxProjVects[iDim][iSolPnt],0);
        computeFlux(m_avgSol,m_tempGrad,m_neighbCellFluxProjVects[side][iDim][iSolPnt],0,m_contFlxNeighb[side][iSolPnt][iDim]);
      }
    }

    for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFlxDep; ++iFlxPnt)
    {
      const CFuint flxIdx = (*m_solFlxDep)[iSolPnt][iFlxPnt];
      const CFuint dim = (*m_flxPntFlxDim)[flxIdx];

      m_extrapolatedFluxes[flxIdx] += (*m_solPolyValsAtFlxPnts)[flxIdx][iSolPnt]*(m_contFlxNeighb[side][iSolPnt][dim]);
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
          residuals[m_nbrEqs*iSolPnt+iEq] += polyCoef*(m_contFlxNeighb[side][jSolIdx][iDir+m_ndimplus][iEq]);
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

void DiffRHSJacobFluxReconstruction::computeUnpertCellDiffResiduals(const CFuint side)
{
  // create a list of the dimensions in which the deriv will be calculated
  for (CFuint iDim = 0; iDim < m_dim+m_ndimplus; ++iDim)
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

  m_dimList.resize(m_dim+m_ndimplus);

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
  
  for (CFuint iDim = 0; iDim < m_dim+m_ndimplus; ++iDim)
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
  m_neighbCellFluxProjVects[LEFT].resize(m_dim+m_ndimplus);
  m_neighbCellFluxProjVects[RIGHT].resize(m_dim+m_ndimplus);
  m_affectedSolPnts.resize(2);
  m_affectedSolPnts[LEFT].resize(m_nbrSolPnts);
  m_affectedSolPnts[RIGHT].resize(m_nbrSolPnts);
  m_contFlxBackup.resize(2);
  m_contFlxNeighb.resize(2);
  m_contFlxBackup[LEFT].resize(m_nbrSolPnts);
  m_contFlxBackup[RIGHT].resize(m_nbrSolPnts);
  m_contFlxNeighb[LEFT].resize(m_nbrSolPnts);
  m_contFlxNeighb[RIGHT].resize(m_nbrSolPnts);

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

  for (CFuint iDim = 0; iDim < m_dim+m_ndimplus; ++iDim)
  {
    m_neighbCellFluxProjVects[LEFT][iDim].resize(m_nbrSolPnts);
    m_neighbCellFluxProjVects[RIGHT][iDim].resize(m_nbrSolPnts);
    for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolPnts; ++iSolPnt)
    {
      m_neighbCellFluxProjVects[LEFT][iDim][iSolPnt].resize(m_dim);
      m_neighbCellFluxProjVects[RIGHT][iDim][iSolPnt].resize(m_dim);
    }
  }
  
  for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolPnts; ++iSolPnt)
  {
    m_contFlxBackup[LEFT][iSolPnt].resize(m_dim+m_ndimplus);
    m_contFlxBackup[RIGHT][iSolPnt].resize(m_dim+m_ndimplus);
    m_contFlxNeighb[LEFT][iSolPnt].resize(m_dim+m_ndimplus);
    m_contFlxNeighb[RIGHT][iSolPnt].resize(m_dim+m_ndimplus);
    
    for (CFuint iDim = 0; iDim < m_dim+m_ndimplus; ++iDim)
    {
      m_contFlxBackup[LEFT][iSolPnt][iDim].resize(m_nbrEqs);
      m_contFlxBackup[RIGHT][iSolPnt][iDim].resize(m_nbrEqs);
      m_contFlxNeighb[LEFT][iSolPnt][iDim].resize(m_nbrEqs);
      m_contFlxNeighb[RIGHT][iSolPnt][iDim].resize(m_nbrEqs);
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

