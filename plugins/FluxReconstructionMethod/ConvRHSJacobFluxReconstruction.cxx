// Copyright (C) 2016 KU Leuven, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/BlockAccumulator.hh"
#include "Framework/LSSMatrix.hh"
#include "Framework/MethodCommandProvider.hh"

#include "Framework/CFSide.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"
#include "Framework/BaseTerm.hh"

#include "MathTools/MathFunctions.hh"

#include "FluxReconstructionMethod/ConvRHSJacobFluxReconstruction.hh"
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

MethodCommandProvider< ConvRHSJacobFluxReconstruction,FluxReconstructionSolverData,FluxReconstructionModule >
  ConvRHSJacobFluxReconstructionProvider("ConvRHSJacob");
  
//////////////////////////////////////////////////////////////////////////////
  
ConvRHSJacobFluxReconstruction::ConvRHSJacobFluxReconstruction(const std::string& name) :
  ConvRHSFluxReconstruction(name),
  m_lss(CFNULL),
  m_numJacob(CFNULL),
  m_acc(CFNULL),
  m_accFace(CFNULL),
  m_pertResUpdates(),
  m_resUpdates(),
  m_derivResUpdates(),
  m_pertCorrections(),
  m_pertDivContFlx(),
  m_cellStatesFlxPntBackup(),
  m_influencedFlxPnt(),
  m_flxPntRiemannFluxBackup(),
  m_pertSol(),
  m_pertVar(),
  m_pertSide(),
  m_extrapolatedFluxesBackup(),
  m_contFlxBackup()
  {
  }
  
//////////////////////////////////////////////////////////////////////////////
  
ConvRHSJacobFluxReconstruction::~ConvRHSJacobFluxReconstruction()
  {
  }

//////////////////////////////////////////////////////////////////////////////

void ConvRHSJacobFluxReconstruction::configure ( Config::ConfigArgs& args )
{
  ConvRHSFluxReconstruction::configure(args);
}  

//////////////////////////////////////////////////////////////////////////////

void ConvRHSJacobFluxReconstruction::execute()
{
  CFAUTOTRACE;
  
  CFLog(VERBOSE, "ConvRHSJacobFluxReconstruction::execute()\n");
  
  // boolean telling whether there is a diffusive term
  const bool hasDiffTerm = getMethodData().hasDiffTerm() || getMethodData().hasArtificialViscosity();
  
  // get the elementTypeData
  SafePtr< vector<ElementTypeData> > elemType = MeshDataStack::getActive()->getElementTypeData();

  // get InnerCells TopologicalRegionSet
  SafePtr<TopologicalRegionSet> cells = MeshDataStack::getActive()->getTrs("InnerCells");

  // get the geodata of the geometric entity builder and set the TRS
  StdTrsGeoBuilder::GeoData& geoDataCell = m_cellBuilder->getDataGE();
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
      
      // if one of the neighbouring cells is parallel updatable or if the gradients have to be computed, set the bnd face data and compute the discontinuous flx
      if ((*m_states[LEFT ])[0]->isParUpdatable() || (*m_states[RIGHT])[0]->isParUpdatable() || hasDiffTerm)
      {
	// set the bnd face data
        setFaceData(m_face->getID());//faceID

	// compute the left and right states in the flx pnts
        computeFlxPntStates();
      }

      // if one of the neighbouring cells is parallel updatable, compute the correction flux
      if ((*m_states[LEFT ])[0]->isParUpdatable() || (*m_states[RIGHT])[0]->isParUpdatable())
      {
	// compute the interface flux
	computeInterfaceFlxCorrection();
	
	// compute the wave speed updates
        computeWaveSpeedUpdates(m_waveSpeedUpd);
	
        // update the wave speed
        updateWaveSpeed();
	
	// compute the correction for the left neighbour
	computeCorrection(LEFT, m_divContFlxL);

	// compute the correction for the right neighbour
	computeCorrection(RIGHT, m_divContFlxR);

	// update RHS
	updateRHSBothSides();
      }
      
      // if there is a diffusive term, compute the gradients
      if (hasDiffTerm)
      {
        computeGradientFaceCorrections();
      }

      // compute the contribution to the numerical jacobian
      if ((*m_states[LEFT])[0]->isParUpdatable() && (*m_states[RIGHT])[0]->isParUpdatable())
      {
        computeBothJacobs();
      }
      else if ((*m_states[LEFT])[0]->isParUpdatable())
      {
        computeOneJacob(LEFT);
      }
      else if ((*m_states[RIGHT])[0]->isParUpdatable())
      {
        computeOneJacob(RIGHT);
      }

      // release the GeometricEntity
      m_faceBuilder->releaseGE();
    }
  }
  
  //// Loop over the elements to calculate the divergence of the continuous flux
  
  // loop over element types, for the moment there should only be one
  const CFuint nbrElemTypes = elemType->size();
  cf_assert(nbrElemTypes == 1);
  for (m_iElemType = 0; m_iElemType < nbrElemTypes; ++m_iElemType)
  {
    // get start and end indexes for this type of element
    const CFuint startIdx = (*elemType)[m_iElemType].getStartIdx();
    const CFuint endIdx   = (*elemType)[m_iElemType].getEndIdx();

    // create blockaccumulator
    m_acc.reset(m_lss->createBlockAccumulator(m_nbrSolPnts,m_nbrSolPnts,m_nbrEqs));

    // loop over cells
    for (CFuint elemIdx = startIdx; elemIdx < endIdx; ++elemIdx)
    {
      // build the GeometricEntity
      geoDataCell.idx = elemIdx;
      m_cell = m_cellBuilder->buildGE();

      // get the states in this cell
      m_cellStates = m_cell->getStates();

      // if the states in the cell are parallel updatable or the gradients need to be computed, set the cell data
      if ((*m_cellStates)[0]->isParUpdatable() || hasDiffTerm)
      {
	// set the cell data
	setCellData();
      }
      
      // if the states in the cell are parallel updatable, compute the divergence of the discontinuous flx (-divFD+divhFD)
      if ((*m_cellStates)[0]->isParUpdatable())
      {
	// compute the residual updates (-divFC)
	computeDivDiscontFlx(m_divContFlx);

	// update RHS
        updateRHS();
      }
      
      // if there is a diffusive term, compute the gradients
      if (hasDiffTerm)
      {
	computeGradients();
      }

      // if the states in the cell are parallel updatable, compute the contribution to the numerical jacobian
      if ((*m_cellStates)[0]->isParUpdatable())
      {
	// add the contributions to the Jacobian
	computeJacobConvCorrection();
      }

      // divide by the Jacobian to transform the residuals back to the physical domain
      //divideByJacobDet();
      
      // print out the residual updates for debugging
      if(m_cell->getID() == 1944)
      {
	CFLog(VERBOSE, "ID  = " << m_cell->getID() << "\n");
        CFLog(VERBOSE, "ConvUpdate = \n");
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

void ConvRHSJacobFluxReconstruction::computeBothJacobs()
{
  // get residual factor
  const CFreal resFactor = getMethodData().getResFactor();

  // dereference accumulator
  BlockAccumulator& acc = *m_accFace;

  // set block row and column indices
  CFuint solIdx = 0;
  for (CFuint iSide = 0; iSide < 2; ++iSide)
  {
    for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol, ++solIdx)
    {
      acc.setRowColIndex(solIdx,(*m_states[iSide])[iSol]->getLocalID());
    }
  }
  
  // put the perturbed and unperturbed corrections in the correct format
  for (CFuint iState = 0; iState < m_nbrSolPnts; ++iState)
  {
    for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
    {
      m_resUpdates[LEFT][m_nbrEqs*iState+iVar] = m_divContFlxL[iState][iVar];
      m_resUpdates[RIGHT][m_nbrEqs*iState+iVar] = m_divContFlxR[iState][iVar];
    }
  }
  
  // store backups of values that will be overwritten
  storeBackups();

  // loop over left and right cell
  for (m_pertSide = 0; m_pertSide < 2; ++m_pertSide)
  {
    // term depending on iSide
    const CFuint pertSideTerm = m_pertSide*m_nbrSolPnts;

    // loop over the states in the left and right cell to perturb the states
    for (m_pertSol = 0; m_pertSol < m_nbrSolPnts; ++m_pertSol)
    {
      // dereference state
      State& pertState = *(*m_states[m_pertSide])[m_pertSol];
      
    // Loop over flux points to determine which flx pnts are influenced by the pert
    if ((m_ndimplus == 3) || (m_ndimplus == 4)) //(elemShape == CFGeoShape::TRIAG or TETRA)
    {
      m_influencedFlxPnt = 0;
      m_NbInfluencedFlxPnts = m_nbrFaceFlxPnts;
    }
    else
    {
      for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
      {
        // get current flx pnt idx
        const CFuint currFlxIdx = (*m_faceFlxPntConnPerOrient)[m_orient][m_pertSide][iFlxPnt];
    
        for (CFuint jFlxPnt = 0; jFlxPnt < m_nbrFlxDep; ++jFlxPnt)
        {
          if (currFlxIdx == (*m_solFlxDep)[m_pertSol][jFlxPnt])
          {
            m_influencedFlxPnt = iFlxPnt;
            m_NbInfluencedFlxPnts= iFlxPnt+1;
          }
        }
      }
    }

      // loop over the variables in the state
      for (m_pertVar = 0; m_pertVar < m_nbrEqs; ++m_pertVar)
      {
        // perturb physical variable in state
        m_numJacob->perturb(m_pertVar,pertState[m_pertVar]);
	
	// compute the perturbed left and right states in the flx pnts
	extrapolatePerturbedState();
	
	// compute perturbed FI
	computePertInterfaceFlxCorrection();
	
	// compute the perturbed corrections
	computePertCorrection(LEFT, m_pertResUpdates[LEFT]);
	computePertCorrection(RIGHT, m_pertResUpdates[RIGHT]);

        // add contributions to the Jacobian
        for (CFuint iSide2 = 0; iSide2 < 2; ++iSide2)
        {
          // factor depending on iSide2
          const CFuint sideTerm = iSide2*m_nbrSolPnts;

          // compute the finite difference derivative of the face term
          m_numJacob->computeDerivative(m_pertResUpdates[iSide2],m_resUpdates[iSide2],m_derivResUpdates);

          // multiply residual update derivatives with residual factor
          m_derivResUpdates *= resFactor;
          
          const CFuint flxIdx = (*m_faceFlxPntConnPerOrient)[m_orient][iSide2][m_influencedFlxPnt];
        
          for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolDep; ++iSolPnt)
          {
            const CFuint solIdx = (*m_flxSolDep)[flxIdx][iSolPnt];
            acc.addValues(solIdx+sideTerm,m_pertSol+pertSideTerm,m_pertVar,&m_derivResUpdates[m_nbrEqs*solIdx]);
          }  
        }

        // restore physical variable in state
        m_numJacob->restore(pertState[m_pertVar]);
        
        // restore values that were overwritten
        restoreFromBackups();
      }
    }
  }
  
//   if (m_cells[LEFT]->getID() == 49 || m_cells[RIGHT]->getID() == 49)
//   {
//   CFLog(VERBOSE,"accFace:\n");
//    acc.printToScreen();
//   }

  if (getMethodData().doComputeJacobian())
  {
    // add the values to the jacobian matrix
    m_lss->getMatrix()->addValues(acc);
  }

  // reset to zero the entries in the block accumulator
  acc.reset();
}

//////////////////////////////////////////////////////////////////////////////

void ConvRHSJacobFluxReconstruction::computeOneJacob(const CFuint side)
{
  // get residual factor
  const CFreal resFactor = getMethodData().getResFactor();

  // dereference accumulator
  BlockAccumulator& acc = *m_accFace;

  // set block row and column indices
  CFuint solIdx = 0;
  for (CFuint iSide = 0; iSide < 2; ++iSide)
  {
    const CFuint nbrSolPnts = m_states[iSide]->size();
    for (CFuint iSol = 0; iSol < nbrSolPnts; ++iSol, ++solIdx)
    {
      acc.setRowColIndex(solIdx,(*m_states[iSide])[iSol]->getLocalID());
    }
  }

  // term depending on the side
  const CFuint sideTerm = side*m_nbrSolPnts;
  
  // put the perturbed and unperturbed corrections in the correct format
  if (side == LEFT)
  {
    for (CFuint iState = 0; iState < m_nbrSolPnts; ++iState)
    {
      for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
      {
        m_resUpdates[side][m_nbrEqs*iState+iVar] = m_divContFlxL[iState][iVar];
      }
    }
  }
  else
  {
    for (CFuint iState = 0; iState < m_nbrSolPnts; ++iState)
    {
      for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
      {
        m_resUpdates[side][m_nbrEqs*iState+iVar] = m_divContFlxR[iState][iVar];
      }
    }
  }
  
  // store backups of values that will be overwritten
  storeBackups();

  // loop over left and right cell
  for (m_pertSide = 0; m_pertSide < 2; ++m_pertSide)
  {
    // term depending on iSide
    const CFuint pertSideTerm = m_pertSide*m_nbrSolPnts;

    // loop over the states to perturb the states
    for (m_pertSol = 0; m_pertSol < m_nbrSolPnts; ++m_pertSol)
    {
      // dereference states
      State& pertState = *(*m_states[m_pertSide])[m_pertSol];
      
      // Loop over flux points to determine which flx pnts are influenced by the pert
      if ((m_ndimplus == 3) || (m_ndimplus == 4)) //(elemShape == CFGeoShape::TRIAG or TETRA)
      {
        m_influencedFlxPnt = 0;
        m_NbInfluencedFlxPnts = m_nbrFaceFlxPnts;
      }
      else
      {
        for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
        {
          // get current flx pnt idx
          const CFuint currFlxIdx = (*m_faceFlxPntConnPerOrient)[m_orient][m_pertSide][iFlxPnt];
      
          for (CFuint jFlxPnt = 0; jFlxPnt < m_nbrFlxDep; ++jFlxPnt)
          {
            if (currFlxIdx == (*m_solFlxDep)[m_pertSol][jFlxPnt])
            {
              m_influencedFlxPnt = iFlxPnt;
              m_NbInfluencedFlxPnts= iFlxPnt+1;
            }
          }
        }
      }

      // loop over the variables in the state
      for (m_pertVar = 0; m_pertVar < m_nbrEqs; ++m_pertVar)
      {
        // perturb physical variable in state
        m_numJacob->perturb(m_pertVar,pertState[m_pertVar]);
	
	// compute the left and right perturbed states in the flx pnts
        extrapolatePerturbedState();
	
	// compute perturbed FI
	computePertInterfaceFlxCorrection();
	
	// compute the perturbed corrections
	computePertCorrection(side, m_pertResUpdates[side]);

        // add contributions to the Jacobian
        // compute the finite difference derivative of the face term
        m_numJacob->computeDerivative(m_pertResUpdates[side],m_resUpdates[side],m_derivResUpdates);

        // multiply residual update derivatives with residual factor
        m_derivResUpdates *= resFactor;
        
        const CFuint flxIdx = (*m_faceFlxPntConnPerOrient)[m_orient][side][m_influencedFlxPnt];
        
        for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolDep; ++iSolPnt)
        {
          const CFuint solIdx = (*m_flxSolDep)[flxIdx][iSolPnt];
          acc.addValues(solIdx+sideTerm,m_pertSol+pertSideTerm,m_pertVar,&m_derivResUpdates[m_nbrEqs*solIdx]);
        }  

        // restore physical variable in state
        m_numJacob->restore(pertState[m_pertVar]);
        
        // restore values that were overwritten
        restoreFromBackups();
      }
    }
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

void ConvRHSJacobFluxReconstruction::extrapolatePerturbedState()
{ 
  for (CFuint iFlxPnt = m_influencedFlxPnt; iFlxPnt < m_NbInfluencedFlxPnts; ++iFlxPnt)
  {
    // local flux point indices in the left and right cell
    const CFuint flxPntIdx = (*m_faceFlxPntConnPerOrient)[m_orient][m_pertSide][iFlxPnt];
      
    // reset states in flx pnt
    (*(m_cellStatesFlxPnt[m_pertSide][iFlxPnt]))[m_pertVar] = 0.0;

    // extrapolate the left and right states to the flx pnts
    for (CFuint iSol = 0; iSol < m_nbrSolDep; ++iSol)
    {
      const CFuint solIdx = (*m_flxSolDep)[flxPntIdx][iSol];
  
      // add the contributions of the current sol pnt
      (*(m_cellStatesFlxPnt[m_pertSide][iFlxPnt]))[m_pertVar] += (*m_solPolyValsAtFlxPnts)[flxPntIdx][solIdx]*(*((*(m_states[m_pertSide]))[solIdx]))[m_pertVar];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvRHSJacobFluxReconstruction::computePertInterfaceFlxCorrection()
{ 
  for (CFuint iFlxPnt = m_influencedFlxPnt; iFlxPnt < m_NbInfluencedFlxPnts; ++iFlxPnt)
  {
    // compute the riemann flux
    m_cellFlx[RIGHT][iFlxPnt] = m_riemannFluxComputer->computeFlux(*(m_cellStatesFlxPnt[LEFT][iFlxPnt]),*(m_cellStatesFlxPnt[RIGHT][iFlxPnt]),m_unitNormalFlxPnts[iFlxPnt]);
    // compute the interface flux in the mapped coord frame and store
    m_cellFlx[LEFT][iFlxPnt] = (m_cellFlx[RIGHT][iFlxPnt])*m_faceJacobVecSizeFlxPnts[iFlxPnt][LEFT];
    m_cellFlx[RIGHT][iFlxPnt] = (m_cellFlx[RIGHT][iFlxPnt])*m_faceJacobVecSizeFlxPnts[iFlxPnt][RIGHT];
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvRHSJacobFluxReconstruction::computePertCorrection(CFuint side, RealVector& corrections)
{ 
  cf_assert(corrections.size() == m_nbrSolPnts*m_nbrEqs);

  // reset corrections
  corrections = 0.0;

  for (CFuint iFlxPnt = m_influencedFlxPnt; iFlxPnt < m_NbInfluencedFlxPnts; ++iFlxPnt)
  {
    const CFuint flxIdx = (*m_faceFlxPntConnPerOrient)[m_orient][side][iFlxPnt];  

    // the current correction factor corresponding to the interface flux (stored in cellFlx)
    const RealVector& currentCorrFactor = m_cellFlx[side][iFlxPnt];

    // loop over sol pnts to compute the corrections
    for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolDep; ++iSolPnt)
    {
      const CFuint solIdx = (*m_flxSolDep)[flxIdx][iSolPnt];

      // divergence of the correction function
      const CFreal divh = m_corrFctDiv[solIdx][flxIdx];

      // Fill in the corrections
      for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
      {
        corrections[m_nbrEqs*solIdx+iVar] -= currentCorrFactor[iVar] * divh; 
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvRHSJacobFluxReconstruction::storeBackups()
{
  for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
  {
    m_cellStatesFlxPntBackup[LEFT][iFlxPnt] = *(m_cellStatesFlxPnt[LEFT][iFlxPnt]);
    m_cellStatesFlxPntBackup[RIGHT][iFlxPnt] = *(m_cellStatesFlxPnt[RIGHT][iFlxPnt]);
    m_flxPntRiemannFluxBackup[LEFT][iFlxPnt] = m_cellFlx[LEFT][iFlxPnt];
    m_flxPntRiemannFluxBackup[RIGHT][iFlxPnt] = m_cellFlx[RIGHT][iFlxPnt];
  }
  
  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {
    m_pertDivContFlx[LEFT][iSol] = m_divContFlxL[iSol];
    m_pertDivContFlx[RIGHT][iSol] = m_divContFlxR[iSol];
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvRHSJacobFluxReconstruction::restoreFromBackups()
{
  for (CFuint iFlxPnt = m_influencedFlxPnt; iFlxPnt < m_NbInfluencedFlxPnts; ++iFlxPnt)
  {
    (*(m_cellStatesFlxPnt[m_pertSide][iFlxPnt]))[m_pertVar] = m_cellStatesFlxPntBackup[m_pertSide][iFlxPnt][m_pertVar];
    m_cellFlx[LEFT][iFlxPnt] = m_flxPntRiemannFluxBackup[LEFT][iFlxPnt];
    m_cellFlx[RIGHT][iFlxPnt] = m_flxPntRiemannFluxBackup[RIGHT][iFlxPnt];
  }    
  const CFuint flxIdxL = (*m_faceFlxPntConnPerOrient)[m_orient][LEFT][m_influencedFlxPnt];
  const CFuint flxIdxR = (*m_faceFlxPntConnPerOrient)[m_orient][RIGHT][m_influencedFlxPnt];
        
  for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolDep; ++iSolPnt)
  {
    const CFuint solIdxL = (*m_flxSolDep)[flxIdxL][iSolPnt];
    const CFuint solIdxR = (*m_flxSolDep)[flxIdxR][iSolPnt];
    m_pertDivContFlx[LEFT][solIdxL] = m_divContFlxL[solIdxL];
    m_pertDivContFlx[RIGHT][solIdxR] = m_divContFlxR[solIdxR];
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvRHSJacobFluxReconstruction::storeBackupsCell()
{
  for (CFuint iFlxPnt = 0; iFlxPnt < m_flxPntsLocalCoords->size(); ++iFlxPnt)
  {
    m_extrapolatedFluxesBackup[iFlxPnt] = m_extrapolatedFluxes[iFlxPnt];
  }
  
  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {
//    m_pertCorrections[iSol] = m_divContFlx[iSol];
    for (CFuint iDim = 0; iDim < m_dim+m_ndimplus; ++iDim)
    {
      m_contFlxBackup[iSol][iDim] = m_contFlx[iSol][iDim];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvRHSJacobFluxReconstruction::restoreFromBackupsCell()
{
  for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFlxDep; ++iFlxPnt)
  {
    const CFuint flxIdx = (*m_solFlxDep)[m_pertSol][iFlxPnt];
    
    m_extrapolatedFluxes[flxIdx] = m_extrapolatedFluxesBackup[flxIdx];
  }
  
//  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
//  { 
//    m_pertCorrections[iSol] = m_divContFlx[iSol];
//  }
  
  for (CFuint iDim = 0; iDim < m_dim+m_ndimplus; ++iDim)
  {
    m_contFlx[m_pertSol][iDim] = m_contFlxBackup[m_pertSol][iDim];
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvRHSJacobFluxReconstruction::computeJacobConvCorrection()
{
  // get residual factor
  const CFreal resFactor = getMethodData().getResFactor();

  // dereference accumulator
  BlockAccumulator& acc = *m_acc;

  // set block row and column indices
  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {
    acc.setRowColIndex(iSol,(*m_cellStates)[iSol]->getLocalID());
    
    for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
    {
      m_resUpdates[0][m_nbrEqs*iSol+iVar] = m_divContFlx[iSol][iVar];
    }
  }
  
  // store backups of values that will be overwritten
  storeBackupsCell();

  // loop over the states/solpnts in this cell to perturb the states
  for (m_pertSol = 0; m_pertSol < m_nbrSolPnts; ++m_pertSol)
  {
    // dereference state
    State& pertState = *(*m_cellStates)[m_pertSol];

    // loop over the variables in the state
    for (m_pertVar = 0; m_pertVar < m_nbrEqs; ++m_pertVar)
    {
      // perturb physical variable in state
      m_numJacob->perturb(m_pertVar,pertState[m_pertVar]);

      // compute the perturbed residual updates (-divFD+divhFD)
      computePertDivDiscontFlx(m_pertResUpdates[0]);

      // compute the finite difference derivative of the volume term
      m_numJacob->computeDerivative(m_pertResUpdates[0],m_resUpdates[0],m_derivResUpdates);

      // multiply residual update derivatives with residual factor
      m_derivResUpdates *= resFactor;
      
      // Loop over affected solution pnts
      for (CFuint jSolPnt = 0; jSolPnt < m_nbrSolSolDep; ++jSolPnt)
      {
        const CFuint jSolIdx = (*m_solSolDep)[m_pertSol][jSolPnt];

        acc.addValues(jSolIdx,m_pertSol,m_pertVar,&m_derivResUpdates[m_nbrEqs*jSolIdx]);
      }

      // restore physical variable in state
      m_numJacob->restore(pertState[m_pertVar]);
      
      // restore overwritten values
      restoreFromBackupsCell();
    }
  }
//   if (m_cell->getID() == 49)
//   {
//   CFLog(VERBOSE,"accVol:\n");
//    acc.printToScreen();
//   }

  if (getMethodData().doComputeJacobian())
  {
    // add the values to the jacobian matrix
    m_lss->getMatrix()->addValues(acc);
  }

  // reset to zero the entries in the block accumulator
  acc.reset();
}

//////////////////////////////////////////////////////////////////////////////

void ConvRHSJacobFluxReconstruction::computePertDivDiscontFlx(RealVector& residuals)
{  
  // reset the extrapolated fluxes
  for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFlxDep; ++iFlxPnt)
  {
    const CFuint flxIdx = (*m_solFlxDep)[m_pertSol][iFlxPnt];
    m_extrapolatedFluxes[flxIdx] = 0.0;
  }

  m_updateVarSet->computePhysicalData(*(*m_cellStates)[m_pertSol], m_pData);

  // calculate the discontinuous flux projected on x, y, z-directions
  for (CFuint iDim = 0; iDim < m_dim+m_ndimplus; ++iDim)
  {
    m_contFlx[m_pertSol][iDim] = m_updateVarSet->getFlux()(m_pData,m_cellFluxProjVects[iDim][m_pertSol]);
  }

  // extrapolate the fluxes to the flux points
  for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFlxDep; ++iFlxPnt)
  {
    const CFuint flxIdx = (*m_solFlxDep)[m_pertSol][iFlxPnt];
    const CFuint dim = (*m_flxPntFlxDim)[flxIdx];
    
    for (CFuint iSol = 0; iSol < m_nbrSolDep; ++iSol)
    {
      const CFuint solIdx = (*m_flxSolDep)[flxIdx][iSol];
      m_extrapolatedFluxes[flxIdx] += (*m_solPolyValsAtFlxPnts)[flxIdx][solIdx]*(m_contFlx[solIdx][dim]);
    }
  } 
  
  // reset the residual updates
  residuals = 0.0;

  // Loop over affected solution pnts
  for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolSolDep; ++iSolPnt)
  {
    const CFuint iSolIdx = (*m_solSolDep)[m_pertSol][iSolPnt];
    
    // Loop over solution pnts to count the factor of all sol pnt polys
    for (CFuint jSolPnt = 0; jSolPnt < m_nbrSolSolDep; ++jSolPnt)
    {
      const CFuint jSolIdx = (*m_solSolDep)[iSolIdx][jSolPnt];

      // Loop over deriv directions and sum them to compute divergence
      for (CFuint iDir = 0; iDir < m_dim; ++iDir)
      {
        const CFreal polyCoef = (*m_solPolyDerivAtSolPnts)[iSolIdx][iDir][jSolIdx]; 

        // Loop over conservative fluxes 
        for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
        {
          // Store divFD in the vector that will be divFC
          residuals[m_nbrEqs*iSolIdx+iEq] -= polyCoef*(m_contFlx[jSolIdx][iDir+m_ndimplus][iEq]);
	}
      }
    }

    // add divhFD to the residual updates
    for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFlxDep; ++iFlxPnt)
    {
      const CFuint flxIdx = (*m_solFlxDep)[iSolIdx][iFlxPnt];

      // get the divergence of the correction function
      const CFreal divh = m_corrFctDiv[iSolIdx][flxIdx];
 
      // Fill in the corrections
      for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
      {
        residuals[m_nbrEqs*iSolIdx+iVar] -= -m_extrapolatedFluxes[flxIdx][iVar] * divh; 
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvRHSJacobFluxReconstruction::setup()
{
  CFAUTOTRACE;
  
  // setup parent class
  ConvRHSFluxReconstruction::setup();
  
  // get the linear system solver
  m_lss = getMethodData().getLinearSystemSolver()[0];

  // get the numerical Jacobian computer
  m_numJacob = getMethodData().getNumericalJacobian();

  // create blockaccumulator
  m_accFace.reset(m_lss->createBlockAccumulator(2*m_nbrSolPnts,2*m_nbrSolPnts,m_nbrEqs));

  // resize variables
  const CFuint resSize = m_nbrSolPnts*m_nbrEqs;
  m_pertResUpdates .resize(2);
  m_pertDivContFlx .resize(2);
  m_cellStatesFlxPntBackup.resize(2);
  m_flxPntRiemannFluxBackup.resize(2);
  m_pertResUpdates [LEFT ].resize(resSize);
  m_pertResUpdates [RIGHT].resize(resSize);
  m_pertDivContFlx [LEFT ].resize(m_nbrSolPnts);
  m_pertDivContFlx [RIGHT].resize(m_nbrSolPnts);
  m_cellStatesFlxPntBackup[LEFT].resize(m_nbrFaceFlxPnts);
  m_cellStatesFlxPntBackup[RIGHT].resize(m_nbrFaceFlxPnts);
  m_flxPntRiemannFluxBackup[LEFT].resize(m_nbrFaceFlxPnts);
  m_flxPntRiemannFluxBackup[RIGHT].resize(m_nbrFaceFlxPnts);
  m_resUpdates .resize(2);
  m_resUpdates [LEFT ].resize(resSize);
  m_resUpdates [RIGHT].resize(resSize);
  m_derivResUpdates.resize(resSize);
  m_pertCorrections.resize(m_nbrSolPnts);
  m_extrapolatedFluxesBackup.resize(m_flxPntsLocalCoords->size());
  m_contFlxBackup.resize(m_nbrSolPnts);
  
  for (CFuint iFlx = 0; iFlx < m_flxPntsLocalCoords->size(); ++iFlx)
  {
    m_extrapolatedFluxesBackup[iFlx].resize(m_nbrEqs);
  }
  
  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    m_cellStatesFlxPntBackup[LEFT][iFlx].resize(m_nbrEqs); 
    m_cellStatesFlxPntBackup[RIGHT][iFlx].resize(m_nbrEqs);
    m_flxPntRiemannFluxBackup[LEFT][iFlx].resize(m_nbrEqs); 
    m_flxPntRiemannFluxBackup[RIGHT][iFlx].resize(m_nbrEqs);
  }
  
  for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolPnts; ++iSolPnt)
  {
    m_pertCorrections[iSolPnt].resize(m_nbrEqs);
    m_pertDivContFlx [LEFT ][iSolPnt].resize(m_nbrEqs);
    m_pertDivContFlx [RIGHT][iSolPnt].resize(m_nbrEqs);
    m_contFlxBackup[iSolPnt].resize(m_dim+m_ndimplus);

    for (CFuint iDim = 0; iDim < m_dim+m_ndimplus; ++iDim)
    {
      m_contFlxBackup[iSolPnt][iDim].resize(m_nbrEqs);
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvRHSJacobFluxReconstruction::unsetup()
{
  CFAUTOTRACE;
  
  // unsetup parent class
  ConvRHSFluxReconstruction::unsetup();
}


//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD

