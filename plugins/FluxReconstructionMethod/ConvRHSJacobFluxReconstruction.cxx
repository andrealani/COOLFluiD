// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
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
  m_divContFlxL(),
  m_divContFlxR(),
  m_pertDivContFlx(),
  m_pertCorrections()
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
  
  CFLog(INFO, "ConvRHSJacobFluxReconstruction::execute()\n");
  
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
      
      computeDiscontinuousFluxes();

      // if one of the neighbouring cells is parallel updatable, compute the correction flux
      if ((*m_states[LEFT ])[0]->isParUpdatable() || (*m_states[RIGHT])[0]->isParUpdatable())
      {
	// compute FI-FD
	computeInterfaceFlxCorrection();
	CFLog(VERBOSE, "real faceID: " << m_face->getID() << ", faceID: " << faceID << "\n");
	
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

	// compute the convective face term contribution to the jacobian
        if ((*m_states[LEFT ])[0]->isParUpdatable() && (*m_states[RIGHT])[0]->isParUpdatable())
        {
          computeBothJacobs();
        }
        else if ((*m_states[LEFT ])[0]->isParUpdatable())
        {
          computeOneJacob(LEFT);
        }
        else if ((*m_states[RIGHT])[0]->isParUpdatable())
        {
          computeOneJacob(RIGHT);
        }

      }
      
      if(m_cells[LEFT]->getID() == 150)
      {
        for (CFuint iState = 0; iState < m_nbrSolPnts; ++iState)
        {
	  
	  CFLog(VERBOSE, "updateFace " << iState << ": " << m_divContFlxL[iState] << "\n");

        }
      }
      if(m_cells[RIGHT]->getID() == 150)
      {
        for (CFuint iState = 0; iState < m_nbrSolPnts; ++iState)
        {
	  
	  CFLog(VERBOSE, "updateFace " << iState << ": " << m_divContFlxR[iState] << "\n");

        }
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

    // create blockaccumulator
    m_acc.reset(m_lss->createBlockAccumulator(m_nbrSolPnts,m_nbrSolPnts,m_nbrEqs));

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

	// add the contributions to the Jacobian
	computeJacobConvCorrection();
      }

      // divide by the Jacobian to transform the residuals back to the physical domain
      //divideByJacobDet();
      
      if(m_cell->getID() == 150)
      {
        for (CFuint iState = 0; iState < m_nbrSolPnts; ++iState)
        {
	  
	  CFLog(VERBOSE, "updateVol " << iState << ": " << m_divContFlx[iState] << "\n");

        }
      }
      
      // print out the residual updates for debugging
      if(m_cell->getID() == 150)
      {
	CFLog(VERBOSE, "ID  = " << (*m_cellStates)[0]->getLocalID() << "\n");
        CFLog(VERBOSE, "UpdateTotal = \n");
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
	  
	  CFLog(VERBOSE, "state " << iState << ": " << *(((*m_cellStates)[iState])->getData()) << "\n");
	  //CFLog(VERBOSE, "coord: " << ((*m_cellStates)[iState])->getCoordinates() << "\n");

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

  // number of solution points in left cell
  const CFuint nbrLeftSolPnts = m_states[LEFT]->size();

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
  
  // put the perturbed and unperturbed corrections in the correct format
  for (CFuint iState = 0; iState < nbrLeftSolPnts; ++iState)
  {
    for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
    {
      m_resUpdates[LEFT][m_nbrEqs*iState+iVar] = m_divContFlxL[iState][iVar];
      m_resUpdates[RIGHT][m_nbrEqs*iState+iVar] = m_divContFlxR[iState][iVar];
    }
  }

  // loop over left and right cell
  for (CFuint iSide = 0; iSide < 2; ++iSide)
  {
    // term depending on iSide
    const CFuint pertSideTerm = iSide*nbrLeftSolPnts;

    // loop over the states in the left and right cell to perturb the states
    const CFuint nbrSolPnts = m_states[iSide]->size();
    for (CFuint iSolPert = 0; iSolPert < nbrSolPnts; ++iSolPert)
    {
      // dereference state
      State& pertState = *(*m_states[iSide])[iSolPert];

      // loop over the variables in the state
      for (CFuint iEqPert = 0; iEqPert < m_nbrEqs; ++iEqPert)
      {
        // perturb physical variable in state
        m_numJacob->perturb(iEqPert,pertState[iEqPert]);

        // backup and reconstruct physical variable in the flux points
        backupAndReconstructPhysVar(iEqPert,*m_states[iSide]);// had arg iSide
	
	// compute the perturbed discontinuous fluxes in the flx pnts
	computeDiscontinuousFluxes();
	
	// compute perturbed FI-FD
	computeInterfaceFlxCorrection();
	
	// compute the perturbed corrections
	computeCorrection(LEFT, m_pertDivContFlx[LEFT]);
	computeCorrection(RIGHT, m_pertDivContFlx[RIGHT]);
      
        // put the perturbed and unperturbed corrections in the correct format
        for (CFuint iState = 0; iState < nbrSolPnts; ++iState)
        {
          for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
	  {
            m_pertResUpdates[LEFT][m_nbrEqs*iState+iVar] = m_pertDivContFlx[LEFT][iState][iVar];
	    m_pertResUpdates[RIGHT][m_nbrEqs*iState+iVar] = m_pertDivContFlx[RIGHT][iState][iVar];
          }
        }

        // add contributions to the Jacobian
        for (CFuint iSide2 = 0; iSide2 < 2; ++iSide2)
        {
          // factor depending on iSide2
          const CFuint sideTerm = iSide2*nbrLeftSolPnts;

          // compute the finite difference derivative of the face term
          m_numJacob->computeDerivative(m_pertResUpdates[iSide2],m_resUpdates[iSide2],m_derivResUpdates);

          // multiply residual update derivatives with residual factor
          m_derivResUpdates *= resFactor;

          // add the derivative of the residual updates to the accumulator
          CFuint resUpdIdx = 0;
          for (CFuint iSol = 0; iSol < nbrSolPnts; ++iSol, resUpdIdx += m_nbrEqs)
          {
            acc.addValues(iSol+sideTerm,iSolPert+pertSideTerm,iEqPert,&m_derivResUpdates[resUpdIdx]);
          }
        }

        // restore physical variable in state
        m_numJacob->restore(pertState[iEqPert]);

        // restore physical variable in the flux points
        restorePhysVar(iEqPert);// had arg iSide
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

  // number of solution points in left cell
  const CFuint nbrLeftSolPnts = m_states[LEFT]->size();

  // term depending on the side
  const CFuint sideTerm = side*nbrLeftSolPnts;
  
  // put the perturbed and unperturbed corrections in the correct format
  for (CFuint iState = 0; iState < nbrLeftSolPnts; ++iState)
  {
    for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
    {
      if (side == LEFT)
      {
        m_resUpdates[side][m_nbrEqs*iState+iVar] = m_divContFlxL[iState][iVar];
      }
      else
      {
	m_resUpdates[side][m_nbrEqs*iState+iVar] = m_divContFlxR[iState][iVar];
      }
    }
  }

  // loop over left and right cell
  for (CFuint iSide = 0; iSide < 2; ++iSide)
  {
    // term depending on iSide
    const CFuint pertSideTerm = iSide*nbrLeftSolPnts;

    // loop over the states to perturb the states
    const CFuint nbrSolPnts = m_states[iSide]->size();
    for (CFuint iSolPert = 0; iSolPert < nbrSolPnts; ++iSolPert)
    {
      // dereference states
      State& pertState = *(*m_states[iSide])[iSolPert];

      // loop over the variables in the state
      for (CFuint iEqPert = 0; iEqPert < m_nbrEqs; ++iEqPert)
      {
        // perturb physical variable in state
        m_numJacob->perturb(iEqPert,pertState[iEqPert]);

        // backup and reconstruct physical variable in the flux points
        backupAndReconstructPhysVar(iEqPert,*m_states[iSide]);// had arg iSide
	
	computeDiscontinuousFluxes();
	
	// compute perturbed FI-FD
	computeInterfaceFlxCorrection();
	
	// compute the perturbed corrections
	computeCorrection(side, m_pertDivContFlx[side]);

        // compute the perturbed face term
        //computeConvFaceTerm(m_pertResUpdates);
	
	// put the perturbed and unperturbed corrections in the correct format
        for (CFuint iState = 0; iState < nbrSolPnts; ++iState)
        {
          for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
	  {
            m_pertResUpdates[side][m_nbrEqs*iState+iVar] = m_pertDivContFlx[side][iState][iVar];
          }
        }

        // add contributions to the Jacobian
        // compute the finite difference derivative of the face term
        m_numJacob->computeDerivative(m_pertResUpdates[side],m_resUpdates[side],m_derivResUpdates);

        // multiply residual update derivatives with residual factor
        m_derivResUpdates *= resFactor;

        // add the derivative of the residual updates to the accumulator
        CFuint resUpdIdx = 0;
        for (CFuint iSol = 0; iSol < nbrSolPnts; ++iSol, resUpdIdx += m_nbrEqs)
        {
          acc.addValues(iSol+sideTerm,iSolPert+pertSideTerm,iEqPert,&m_derivResUpdates[resUpdIdx]);
        }

        // restore physical variable in state
        m_numJacob->restore(pertState[iEqPert]);

        // restore physical variable in the flux points
        restorePhysVar(iEqPert);// had arg iSide
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

void ConvRHSJacobFluxReconstruction::computeJacobConvCorrection()
{ 
  // get number of solution points in this cell
  const CFuint nbrSolPnts = m_cellStates->size();

  // get residual factor
  const CFreal resFactor = getMethodData().getResFactor();

  // dereference accumulator
  BlockAccumulator& acc = *m_acc;


  // set block row and column indices
  for (CFuint iSol = 0; iSol < nbrSolPnts; ++iSol)
  {
    acc.setRowColIndex(iSol,(*m_cellStates)[iSol]->getLocalID());
  }
  
  // put the perturbed and unperturbed corrections in the correct format
  for (CFuint iState = 0; iState < nbrSolPnts; ++iState)
  {
    for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
    {
      m_resUpdates[0][m_nbrEqs*iState+iVar] = m_divContFlx[iState][iVar];
    }
  }

  // loop over the states/solpnts in this cell to perturb the states
  for (CFuint iSolPert = 0; iSolPert < nbrSolPnts; ++iSolPert)
  {
    // dereference state
    State& pertState = *(*m_cellStates)[iSolPert];

    // loop over the variables in the state
    for (CFuint iEqPert = 0; iEqPert < m_nbrEqs; ++iEqPert)
    {
      // perturb physical variable in state
      m_numJacob->perturb(iEqPert,pertState[iEqPert]);

      // backup and reconstruct physical variable in the flux points
      backupAndReconstructPhysVar(iEqPert,*m_cellStates);

      // compute the perturbed residual updates (-divFC)
      computeResUpdates(m_pertCorrections);
      
      // put the perturbed and unperturbed corrections in the correct format
      for (CFuint iState = 0; iState < nbrSolPnts; ++iState)
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

      // add the derivative of the residual updates to the accumulator
      CFuint resUpdIdx = 0;
      for (CFuint iSol = 0; iSol < nbrSolPnts; ++iSol, resUpdIdx += m_nbrEqs)
      {
        acc.addValues(iSol,iSolPert,iEqPert,&m_derivResUpdates[resUpdIdx]);
      }

      // restore physical variable in state
      m_numJacob->restore(pertState[iEqPert]);

      // restore physical variable in the flux points
      restorePhysVar(iEqPert);
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

void ConvRHSJacobFluxReconstruction::backupAndReconstructPhysVar
    (const CFuint iVar, const vector< State* >& cellIntStates)
{
  // backup
  backupPhysVar(iVar);
  

}

//////////////////////////////////////////////////////////////////////////////

void ConvRHSJacobFluxReconstruction::backupPhysVar(const CFuint iVar)
{
  //cf_assert(m_nbrFlxPnts[m_orient] <= m_backupPhysVar.size());
  //cf_assert(m_nbrFlxPnts[m_orient] <= m_backupGhostStates.size());
  //for (CFuint iFlx = 0; iFlx < m_nbrFlxPnts[m_orient]; ++iFlx)
  //{
    //m_backupPhysVar[iFlx] = (*m_flxPntIntSol[iFlx])[iVar];
    //m_backupGhostStates[iFlx] = *m_flxPntGhostSol[iFlx];
  //}
}

//////////////////////////////////////////////////////////////////////////////

void ConvRHSJacobFluxReconstruction::restorePhysVar(const CFuint iVar)
{
  //cf_assert(m_nbrFlxPnts[m_orient] <= m_backupPhysVar.size());
  //cf_assert(m_nbrFlxPnts[m_orient] <= m_backupGhostStates.size());
  //for (CFuint iFlx = 0; iFlx < m_nbrFlxPnts[m_orient]; ++iFlx)
  //{
    //(*m_flxPntIntSol[iFlx])[iVar] = m_backupPhysVar[iFlx];
    //*m_flxPntGhostSol[iFlx] = m_backupGhostStates[iFlx];
  //}
}

//////////////////////////////////////////////////////////////////////////////

void ConvRHSJacobFluxReconstruction::setup()
{
  CFAUTOTRACE;
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
  m_pertResUpdates [LEFT ].resize(resSize);
  m_pertResUpdates [RIGHT].resize(resSize);
  m_pertDivContFlx [LEFT ].resize(m_nbrSolPnts);
  m_pertDivContFlx [RIGHT].resize(m_nbrSolPnts);
  m_resUpdates .resize(2);
  m_resUpdates [LEFT ].resize(resSize);
  m_resUpdates [RIGHT].resize(resSize);
  m_derivResUpdates.resize(resSize);
  
  m_divContFlxL.resize(m_nbrSolPnts);
  m_divContFlxR.resize(m_nbrSolPnts);
  
  m_pertCorrections.resize(m_nbrSolPnts);
  
  for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolPnts; ++iSolPnt)
  {
    m_divContFlxL[iSolPnt].resize(m_nbrEqs);
    m_divContFlxR[iSolPnt].resize(m_nbrEqs);
    m_pertCorrections[iSolPnt].resize(m_nbrEqs);
    m_pertDivContFlx [LEFT ][iSolPnt].resize(m_nbrEqs);
    m_pertDivContFlx [RIGHT][iSolPnt].resize(m_nbrEqs);
  }

}

//////////////////////////////////////////////////////////////////////////////

void ConvRHSJacobFluxReconstruction::unsetup()
{
  CFAUTOTRACE;
  ConvRHSFluxReconstruction::unsetup();
}


//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD

