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

#include "FluxReconstructionMethod/ConvRHSJacobFluxReconstructionBlending.hh"
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

MethodCommandProvider< ConvRHSJacobFluxReconstructionBlending,FluxReconstructionSolverData,FluxReconstructionModule >
  ConvRHSJacobFluxReconstructionBlendingProvider("ConvRHSJacobBlending");
  
//////////////////////////////////////////////////////////////////////////////
  
ConvRHSJacobFluxReconstructionBlending::ConvRHSJacobFluxReconstructionBlending(const std::string& name) :
  ConvRHSFluxReconstructionBlending(name),
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
  m_influencedFlxPnts(),
  m_flxPntRiemannFluxBackup(),
  m_pertSol(),
  m_pertVar(),
  m_pertSide(),
  m_extrapolatedFluxesBackup(),
  m_contFlxBackup(),
  elemShape()
  {
  }
  
//////////////////////////////////////////////////////////////////////////////
  
ConvRHSJacobFluxReconstructionBlending::~ConvRHSJacobFluxReconstructionBlending()
  {
  }

//////////////////////////////////////////////////////////////////////////////

void ConvRHSJacobFluxReconstructionBlending::configure ( Config::ConfigArgs& args )
{
  ConvRHSFluxReconstructionBlending::configure(args);
}  

//////////////////////////////////////////////////////////////////////////////

void ConvRHSJacobFluxReconstructionBlending::execute()
{
  CFAUTOTRACE;
  
  CFLog(VERBOSE, "ConvRHSJacobFluxReconstructionBlending::execute()\n");
  
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

    // Reset the value of m_nbrFaceFlxPnts in case it is not the same for all faces (Prism)
    m_nbrFaceFlxPnts = (*m_faceFlxPntConnPerOrient)[m_orient][0].size();

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

      /////// NEW HERE --- P0 state ///////
      RealMatrix transformationMatrixL;
      RealMatrix filterL;

      filterL.resize(m_nbrSolPnts, m_nbrSolPnts);
      transformationMatrixL.resize(m_nbrSolPnts, m_nbrSolPnts);

      // Construct the filter matrix 
      filterL = 0.0;
      filterL(0, 0) = 1.0;


      // Compute the transformation matrix
    
      transformationMatrixL = m_vdm * filterL * m_vdmInv;
    
        // Applying filter
      for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq) 
      {
        for (CFuint i = 0; i < m_nbrSolPnts; ++i) 
        { 
          m_statesP0[LEFT][i][iEq] = 0.;
          m_statesP0[RIGHT][i][iEq] = 0.;
          for (CFuint j = 0; j < m_nbrSolPnts; ++j) 
          { 
            m_statesP0[LEFT][i][iEq] += transformationMatrixL(i, j) * (*((*(m_states[LEFT]))[j]))[iEq];
            m_statesP0[RIGHT][i][iEq] += transformationMatrixL(i, j) * (*((*(m_states[RIGHT]))[j]))[iEq];
          }
        }
      }

      for (CFuint iState = 0; iState < 1; ++iState)
      {
        for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
        {
          (*((*(m_states_P0[LEFT]))[iState]))[iVar] =m_statesP0[LEFT][iState][iVar];
          (*((*(m_states_P0[RIGHT]))[iState]))[iVar] =m_statesP0[RIGHT][iState][iVar];
        }
      }
    
      //////////////////////////////////////
      
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
        computeWaveSpeedUpdates(m_waveSpeedUpd, m_waveSpeedUpdP0);
	
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

      
        /////// NEW HERE --- P0 state ///////
      
        RealMatrix transformationMatrixL;
        RealMatrix filterL;
      
        filterL.resize(m_nbrSolPnts, m_nbrSolPnts);
        transformationMatrixL.resize(m_nbrSolPnts, m_nbrSolPnts);

        // Construct the filter matrix 
        filterL = 0.0;
        filterL(0, 0) = 1.0;

        // Compute the transformation matrix
      
        transformationMatrixL = m_vdm * filterL * m_vdmInv;
      
          // Applying filter
        for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq) 
        {
          for (CFuint i = 0; i < m_nbrSolPnts; ++i) 
          { 
            m_P0State[i][iEq] = 0.;
            for (CFuint j = 0; j < m_nbrSolPnts; ++j) 
            { 
              m_P0State[i][iEq] += transformationMatrixL(i, j) * (*((*(m_cellStates))[j]))[iEq];
            }
          }
        }

      for (CFuint iState = 0; iState < 1; ++iState)
      {
        for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
        {
          (*((*(m_cellStatesP0))[iState]))[iVar] =m_P0State[iState][iVar];
        }
      }

        //////////////////////////////////////
       
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

void ConvRHSJacobFluxReconstructionBlending::computeBothJacobs()
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

      State& pertStateP0 = *(*m_states_P0[m_pertSide])[0];
      
    // Loop over flux points to determine which flx pnts are influenced by the pert
    m_influencedFlxPnts.resize(0);
    for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
    {
      // get current flx pnt idx
      const CFuint currFlxIdx = (*m_faceFlxPntConnPerOrient)[m_orient][m_pertSide][iFlxPnt];
    
      for (CFuint jFlxPnt = 0; jFlxPnt < m_nbrFlxDep; ++jFlxPnt)
      {
        if (currFlxIdx == (*m_solFlxDep)[m_pertSol][jFlxPnt])
        {
          //m_influencedFlxPnt = iFlxPnt;
          m_influencedFlxPnts.push_back(iFlxPnt);
        }
      }
    }
    m_NbInfluencedFlxPnts= m_influencedFlxPnts.size(); 

      // loop over the variables in the state
      for (m_pertVar = 0; m_pertVar < m_nbrEqs; ++m_pertVar)
      {

        m_numJacob->perturb(m_pertVar,pertStateP0[m_pertVar]);
        for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
        {
          m_statesP0[m_pertSide][0][iVar]= (*((*(m_states_P0[m_pertSide]))[0]))[iVar];
        }
        m_numJacob->restore(pertStateP0[m_pertVar]);


        // perturb physical variable in state
        m_numJacob->perturb(m_pertVar,pertState[m_pertVar]);
	
	// compute the perturbed left and right states in the flx pnts
	extrapolatePerturbedState();
	
	// compute perturbed FI
	computePertInterfaceFlxCorrection();
	
	// compute the perturbed corrections
	computePertCorrection(LEFT, m_pertResUpdates[LEFT]);
	computePertCorrection(RIGHT, m_pertResUpdates[RIGHT]);

//for (CFuint iFlxPnt = 0; iFlxPnt < m_NbInfluencedFlxPnts; ++iFlxPnt)
         // {
        // add contributions to the Jacobian
        for (CFuint iSide2 = 0; iSide2 < 2; ++iSide2)
        {
          // factor depending on iSide2
          const CFuint sideTerm = iSide2*m_nbrSolPnts;

          // compute the finite difference derivative of the face term
          m_numJacob->computeDerivative(m_pertResUpdates[iSide2],m_resUpdates[iSide2],m_derivResUpdates);

          // multiply residual update derivatives with residual factor
          m_derivResUpdates *= resFactor;
          
          const CFuint flxIdx = (*m_faceFlxPntConnPerOrient)[m_orient][iSide2][ m_influencedFlxPnts[0] ];
        m_nbrSolDep = ((*m_flxSolDep)[flxIdx]).size();
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

void ConvRHSJacobFluxReconstructionBlending::computeOneJacob(const CFuint side)
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

      State& pertStateP0 = *(*m_states_P0[m_pertSide])[0];
      
      // Loop over flux points to determine which flx pnts are influenced by the pert
      m_influencedFlxPnts.resize(0);
      for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
      {
        // get current flx pnt idx
        const CFuint currFlxIdx = (*m_faceFlxPntConnPerOrient)[m_orient][m_pertSide][iFlxPnt];
      
        for (CFuint jFlxPnt = 0; jFlxPnt < m_nbrFlxDep; ++jFlxPnt)
        {
          if (currFlxIdx == (*m_solFlxDep)[m_pertSol][jFlxPnt])
          {
            //m_influencedFlxPnt = iFlxPnt;
            m_influencedFlxPnts.push_back(iFlxPnt);
          }
        }
      }
      m_NbInfluencedFlxPnts= m_influencedFlxPnts.size(); 

      // loop over the variables in the state
      for (m_pertVar = 0; m_pertVar < m_nbrEqs; ++m_pertVar)
      {

        m_numJacob->perturb(m_pertVar,pertStateP0[m_pertVar]);
        for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
        {
          m_statesP0[m_pertSide][0][iVar]= (*((*(m_states_P0[m_pertSide]))[0]))[iVar];
        }
        m_numJacob->restore(pertStateP0[m_pertVar]);


        // perturb physical variable in state
        m_numJacob->perturb(m_pertVar,pertState[m_pertVar]);
	
	// compute the left and right perturbed states in the flx pnts
        extrapolatePerturbedState();
	
	// compute perturbed FI
	computePertInterfaceFlxCorrection();
	
	// compute the perturbed corrections
	computePertCorrection(side, m_pertResUpdates[side]);
  
  //for (CFuint iFlxPnt = 0; iFlxPnt < m_NbInfluencedFlxPnts; ++iFlxPnt)
       // {
        // add contributions to the Jacobian
        // compute the finite difference derivative of the face term
        m_numJacob->computeDerivative(m_pertResUpdates[side],m_resUpdates[side],m_derivResUpdates);

        // multiply residual update derivatives with residual factor
        m_derivResUpdates *= resFactor;

        const CFuint flxIdx = (*m_faceFlxPntConnPerOrient)[m_orient][side][ m_influencedFlxPnts[0] ];
        m_nbrSolDep = ((*m_flxSolDep)[flxIdx]).size();
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

void ConvRHSJacobFluxReconstructionBlending::extrapolatePerturbedState()
{ 
  for (CFuint iFlxPnt = 0; iFlxPnt < m_NbInfluencedFlxPnts; ++iFlxPnt)
  {
    // local flux point indices in the left and right cell
    const CFuint flxPntIdx = (*m_faceFlxPntConnPerOrient)[m_orient][m_pertSide][ m_influencedFlxPnts[iFlxPnt] ];
      
    // reset states in flx pnt
    (*(m_cellStatesFlxPnt[m_pertSide][ m_influencedFlxPnts[iFlxPnt] ]))[m_pertVar] = 0.0;
m_nbrSolDep = ((*m_flxSolDep)[flxPntIdx]).size();
    // extrapolate the left and right states to the flx pnts
    for (CFuint iSol = 0; iSol < m_nbrSolDep; ++iSol)
    {
      const CFuint solIdx = (*m_flxSolDep)[flxPntIdx][iSol];
  
      // add the contributions of the current sol pnt
      (*(m_cellStatesFlxPnt[m_pertSide][ m_influencedFlxPnts[iFlxPnt] ]))[m_pertVar] += (*m_solPolyValsAtFlxPnts)[flxPntIdx][solIdx]*(*((*(m_states[m_pertSide]))[solIdx]))[m_pertVar];
    }
  }

  for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPntsP0; ++iFlxPnt)
  {
    // local flux point indices in the left and right cell

    const CFuint flxPntIdx = (*m_faceFlxPntConnPerOrientP0)[m_orient][m_pertSide][ iFlxPnt ];
      
    // reset states in flx pnt
    (*(m_cellStatesFlxPntP0[m_pertSide][iFlxPnt ]))[m_pertVar] = 0.0;
    
    m_nbrSolDep = ((*m_flxSolDepP0)[flxPntIdx]).size();
    // extrapolate the left and right states to the flx pnts
    for (CFuint iSol = 0; iSol < m_nbrSolDep; ++iSol)
    {
      const CFuint solIdx = (*m_flxSolDepP0)[flxPntIdx][iSol];
  
      // add the contributions of the current sol pnt
      (*(m_cellStatesFlxPntP0[m_pertSide][ iFlxPnt ]))[m_pertVar] += (*m_solPolyValsAtFlxPntsP0)[flxPntIdx][solIdx]*((((m_statesP0[m_pertSide]))[solIdx]))[m_pertVar];
    }
  }

}

//////////////////////////////////////////////////////////////////////////////

void ConvRHSJacobFluxReconstructionBlending::computePertInterfaceFlxCorrection()
{ 
  for (CFuint iFlxPnt = 0; iFlxPnt < m_NbInfluencedFlxPnts; ++iFlxPnt)
  {
    // compute the riemann flux
    m_cellFlx[RIGHT][ m_influencedFlxPnts[iFlxPnt] ] = m_riemannFluxComputer->computeFlux(*(m_cellStatesFlxPnt[LEFT][ m_influencedFlxPnts[iFlxPnt] ]),*(m_cellStatesFlxPnt[RIGHT][ m_influencedFlxPnts[iFlxPnt] ]),m_unitNormalFlxPnts[ m_influencedFlxPnts[iFlxPnt] ]);
    // compute the interface flux in the mapped coord frame and store
    m_cellFlx[LEFT][ m_influencedFlxPnts[iFlxPnt] ] = (m_cellFlx[RIGHT][ m_influencedFlxPnts[iFlxPnt] ])*m_faceJacobVecSizeFlxPnts[ m_influencedFlxPnts[iFlxPnt] ][LEFT];
    m_cellFlx[RIGHT][ m_influencedFlxPnts[iFlxPnt] ] = (m_cellFlx[RIGHT][ m_influencedFlxPnts[iFlxPnt] ])*m_faceJacobVecSizeFlxPnts[ m_influencedFlxPnts[iFlxPnt] ][RIGHT];
  }

  for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPntsP0; ++iFlxPnt)
  {
    // compute the riemann flux
    m_cellFlxP0[RIGHT][ iFlxPnt ] = m_riemannFluxComputer->computeFlux(*(m_cellStatesFlxPntP0[LEFT][ iFlxPnt ]),*(m_cellStatesFlxPntP0[RIGHT][ iFlxPnt ]),m_unitNormalFlxPntsP0[ iFlxPnt ]);
    // compute the interface flux in the mapped coord frame and store
    m_cellFlxP0[LEFT][ iFlxPnt ] = (m_cellFlxP0[RIGHT][ iFlxPnt ])*m_faceJacobVecSizeFlxPntsP0[ iFlxPnt ][LEFT];
    m_cellFlxP0[RIGHT][ iFlxPnt ] = (m_cellFlxP0[RIGHT][ iFlxPnt ])*m_faceJacobVecSizeFlxPntsP0[ iFlxPnt ][RIGHT];
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvRHSJacobFluxReconstructionBlending::computePertCorrection(CFuint side, RealVector& corrections)
{ 
  cf_assert(corrections.size() == m_nbrSolPnts*m_nbrEqs);

  DataHandle< CFreal > output = socket_alpha.getDataHandle();
  CFreal alpha = output[(*m_states[side])[0]->getLocalID()];

  // reset corrections
  corrections = 0.0;

  // STEP 1: Add high-order contributions (loop over flux points)
  for (CFuint iFlxPnt = 0; iFlxPnt < m_NbInfluencedFlxPnts; ++iFlxPnt)
  {
    const CFuint flxIdx = (*m_faceFlxPntConnPerOrient)[m_orient][side][ m_influencedFlxPnts[iFlxPnt] ];  

    // the current correction factor corresponding to the interface flux (stored in cellFlx)
    // This is the high-order flux correction
    const RealVector& currentCorrFactor = m_cellFlx[side][ m_influencedFlxPnts[iFlxPnt] ];
    
    m_nbrSolDep = ((*m_flxSolDep)[flxIdx]).size();
    // loop over sol pnts affected by this flux point to compute the corrections
    for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolDep; ++iSolPnt)
    {
      const CFuint solIdx = (*m_flxSolDep)[flxIdx][iSolPnt];

      // divergence of the correction function
      const CFreal divh = m_corrFctDiv[solIdx][flxIdx];

      // Fill in the high-order corrections with blending factor (1-alpha)
      for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
      {
        corrections[m_nbrEqs*solIdx+iVar] -= currentCorrFactor[iVar] * divh * (1.0 - alpha); 
      }
    }
  }
  
  // STEP 2: Add P0 contribution (only once per solution point)
  // This applies to all solution points in the cell
  const CFuint flxIdx = (*m_faceFlxPntConnPerOrient)[m_orient][side][0];
  const CFuint flxIdxP0 = (*m_faceFlxPntConnP0)[m_flxPntFaceConn[flxIdx]][0]; 
  const RealVector& currentCorrFactorP0 = m_cellFlxP0[side][0];
  
  for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolPnts; ++iSolPnt)
  {
    const CFreal divhP0 = m_corrFctDivP0[0][flxIdxP0];
    for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
    {
      corrections[m_nbrEqs*iSolPnt+iVar] -= currentCorrFactorP0[iVar] * divhP0 * alpha;// * (*m_cellAvgSolCoefs)[iSolPnt];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvRHSJacobFluxReconstructionBlending::storeBackups()
{
  for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
  {
    m_cellStatesFlxPntBackup[LEFT][iFlxPnt] = *(m_cellStatesFlxPnt[LEFT][iFlxPnt]);
    m_cellStatesFlxPntBackup[RIGHT][iFlxPnt] = *(m_cellStatesFlxPnt[RIGHT][iFlxPnt]);
    m_flxPntRiemannFluxBackup[LEFT][iFlxPnt] = m_cellFlx[LEFT][iFlxPnt];
    m_flxPntRiemannFluxBackup[RIGHT][iFlxPnt] = m_cellFlx[RIGHT][iFlxPnt];
  }

  for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPntsP0; ++iFlxPnt)
  {
    m_cellStatesFlxPntBackupP0[LEFT][iFlxPnt] = *(m_cellStatesFlxPntP0[LEFT][iFlxPnt]);
    m_cellStatesFlxPntBackupP0[RIGHT][iFlxPnt] = *(m_cellStatesFlxPntP0[RIGHT][iFlxPnt]);
    m_flxPntRiemannFluxBackupP0[LEFT][iFlxPnt] = m_cellFlxP0[LEFT][iFlxPnt];
    m_flxPntRiemannFluxBackupP0[RIGHT][iFlxPnt] = m_cellFlxP0[RIGHT][iFlxPnt];
  }
  
  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {
    m_pertDivContFlx[LEFT][iSol] = m_divContFlxL[iSol];
    m_pertDivContFlx[RIGHT][iSol] = m_divContFlxR[iSol];
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvRHSJacobFluxReconstructionBlending::restoreFromBackups()
{
  for (CFuint iFlxPnt = 0; iFlxPnt < m_NbInfluencedFlxPnts; ++iFlxPnt)
  {
    (*(m_cellStatesFlxPnt[m_pertSide][ m_influencedFlxPnts[iFlxPnt] ]))[m_pertVar] = m_cellStatesFlxPntBackup[m_pertSide][ m_influencedFlxPnts[iFlxPnt] ][m_pertVar];
    m_cellFlx[LEFT][ m_influencedFlxPnts[iFlxPnt] ] = m_flxPntRiemannFluxBackup[LEFT][ m_influencedFlxPnts[iFlxPnt] ];
    m_cellFlx[RIGHT][ m_influencedFlxPnts[iFlxPnt] ] = m_flxPntRiemannFluxBackup[RIGHT][ m_influencedFlxPnts[iFlxPnt] ];
  }

  for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPntsP0; ++iFlxPnt)
  {    
    (*(m_cellStatesFlxPntP0[m_pertSide][ iFlxPnt ]))[m_pertVar] = m_cellStatesFlxPntBackupP0[m_pertSide][ iFlxPnt ][m_pertVar];
    m_cellFlxP0[LEFT][ iFlxPnt ] = m_flxPntRiemannFluxBackupP0[LEFT][ iFlxPnt ];
    m_cellFlxP0[RIGHT][ iFlxPnt ] = m_flxPntRiemannFluxBackupP0[RIGHT][ iFlxPnt ];
  }
  const CFuint flxIdxL = (*m_faceFlxPntConnPerOrient)[m_orient][LEFT][ m_influencedFlxPnts[0] ];
  const CFuint flxIdxR = (*m_faceFlxPntConnPerOrient)[m_orient][RIGHT][ m_influencedFlxPnts[0] ];
  
  m_nbrSolDep = ((*m_flxSolDep)[flxIdxL]).size();
  for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolDep; ++iSolPnt)
  {
    const CFuint solIdxL = (*m_flxSolDep)[flxIdxL][iSolPnt];
    const CFuint solIdxR = (*m_flxSolDep)[flxIdxR][iSolPnt];
    m_pertDivContFlx[LEFT][solIdxL] = m_divContFlxL[solIdxL];
    m_pertDivContFlx[RIGHT][solIdxR] = m_divContFlxR[solIdxR];
  }
  
}

//////////////////////////////////////////////////////////////////////////////

void ConvRHSJacobFluxReconstructionBlending::storeBackupsCell()
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

  /// only one sol pnt in P0
  for (CFuint iFlxPnt = 0; iFlxPnt < m_flxPntsLocalCoordsP0->size(); ++iFlxPnt)
  {
    m_extrapolatedFluxesBackupP0[iFlxPnt] = m_extrapolatedFluxesP0[iFlxPnt];
  }
  
  //for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  //{
//    m_pertCorrections[iSol] = m_divContFlx[iSol];
    for (CFuint iDim = 0; iDim < m_dim+m_ndimplus; ++iDim)
    {
      m_contFlxBackupP0[0][iDim] = m_contFlxP0[0][iDim];
    }
  //}
}

//////////////////////////////////////////////////////////////////////////////

void ConvRHSJacobFluxReconstructionBlending::restoreFromBackupsCell()
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

  /// only one sol pnt in P0
  for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFlxDepP0; ++iFlxPnt)
  {
    const CFuint flxIdx = (*m_solFlxDepP0)[0][iFlxPnt];
    
    m_extrapolatedFluxesP0[flxIdx] = m_extrapolatedFluxesBackupP0[flxIdx];
  }
  
  for (CFuint iDim = 0; iDim < m_dim+m_ndimplus; ++iDim)
  {
    m_contFlxP0[0][iDim] = m_contFlxBackupP0[0][iDim];
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvRHSJacobFluxReconstructionBlending::computeJacobConvCorrection()
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

    State& pertStateP0 = *(*m_cellStatesP0)[0]; //m_P0State[0]; //only one sol pnt per P0 element (sqme contribution to all local sol pnts)

    // loop over the variables in the state
    for (m_pertVar = 0; m_pertVar < m_nbrEqs; ++m_pertVar)
    {

      m_numJacob->perturb(m_pertVar,pertStateP0[m_pertVar]);
      for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
      {
        (*((*(m_PertcellStatesP0))[0]))[iVar] = (*((*(m_cellStatesP0))[0]))[iVar];
      }
      m_numJacob->restore(pertStateP0[m_pertVar]);


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

void ConvRHSJacobFluxReconstructionBlending::computePertDivDiscontFlx(RealVector& residuals)
{  
  DataHandle< CFreal > output = socket_alpha.getDataHandle();
  CFreal alpha = output[((*m_cellStates)[0])->getLocalID()];

    ////P0
    for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFlxDepP0; ++iFlxPnt)
    {
      const CFuint flxIdx = (*m_solFlxDepP0)[0][iFlxPnt];
      m_extrapolatedFluxesP0[flxIdx] = 0.0;
    }
  
    m_updateVarSet->computePhysicalData(*(*m_PertcellStatesP0)[0], m_pData);
  
    // calculate the discontinuous flux projected on x, y, z-directions
    for (CFuint iDim = 0; iDim < m_dim+m_ndimplus; ++iDim)
    {
      m_contFlxP0[0][iDim] = m_updateVarSet->getFlux()(m_pData,m_cellFluxProjVectsP0[iDim][0]);
    }
  
    // extrapolate the fluxes to the flux points
    for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFlxDepP0; ++iFlxPnt)
    {
      const CFuint flxIdx = (*m_solFlxDepP0)[0][iFlxPnt];
      const CFuint dim = (*m_flxPntFlxDimP0)[flxIdx];
      m_extrapolatedFluxesP0[flxIdx] += (*m_solPolyValsAtFlxPntsP0)[flxIdx][0]*(m_contFlxP0[0][dim]);
    }
    ////
  
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
    m_nbrSolDep = ((*m_flxSolDep)[flxIdx]).size();
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

        const CFreal polyCoefP0 = (*m_solPolyDerivAtSolPntsP0)[0][iDir][0];

        // Loop over conservative fluxes 
        for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
        {
          // Store divFD in the vector that will be divFC
          residuals[m_nbrEqs*iSolIdx+iEq] -= polyCoef*(m_contFlx[jSolIdx][iDir][iEq]) * (1.-alpha);// + polyCoefP0*(m_contFlxP0[0][iDir][iEq]) * alpha;
	}
      }
    }
        // add P0 divFD to the residual updates
    for (CFuint iDir = 0; iDir < m_dim; ++iDir)
    {
      const CFreal polyCoefP0 = (*m_solPolyDerivAtSolPntsP0)[0][iDir][0];
      for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
      {
        residuals[m_nbrEqs*iSolIdx+iEq] -= polyCoefP0*(m_contFlxP0[0][iDir][iEq]) * alpha;// * (*m_cellAvgSolCoefs)[iSolPnt];
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
        residuals[m_nbrEqs*iSolIdx+iVar] -= -m_extrapolatedFluxes[flxIdx][iVar] * divh * (1.-alpha);
      }
    }

        // add P0 divhFD to the residual updates
    for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFlxDepP0; ++iFlxPnt)
    {
      const CFuint flxIdx = (*m_solFlxDepP0)[0][iFlxPnt];

      // get the divergence of the correction function
      const CFreal divhP0 = m_corrFctDivP0[0][flxIdx];

      // Fill in the corrections
      for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
      {
        residuals[m_nbrEqs*iSolIdx+iVar] -= -m_extrapolatedFluxesP0[flxIdx][iVar] * divhP0 * alpha;// * (*m_cellAvgSolCoefs)[iSolPnt];
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvRHSJacobFluxReconstructionBlending::setup()
{
  CFAUTOTRACE;
  
  // setup parent class
  ConvRHSFluxReconstructionBlending::setup();
  
  // get the linear system solver
  m_lss = getMethodData().getLinearSystemSolver()[0];

  // get the numerical Jacobian computer
  m_numJacob = getMethodData().getNumericalJacobian();

  m_numJacob_P0 = getMethodData().getNumericalJacobian();

  // create blockaccumulator
  m_accFace.reset(m_lss->createBlockAccumulator(2*m_nbrSolPnts,2*m_nbrSolPnts,m_nbrEqs));

  // get the local FR data
  vector< FluxReconstructionElementData* >& frLocalData = getMethodData().getFRLocalData();

  //get element shape
  elemShape = frLocalData[0]->getShape();

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


  m_cellStatesFlxPntBackupP0.resize(2);
  m_flxPntRiemannFluxBackupP0.resize(2);
  m_cellStatesFlxPntBackupP0[LEFT].resize(m_nbrFaceFlxPntsP0);
  m_cellStatesFlxPntBackupP0[RIGHT].resize(m_nbrFaceFlxPntsP0);
  m_flxPntRiemannFluxBackupP0[LEFT].resize(m_nbrFaceFlxPntsP0);
  m_flxPntRiemannFluxBackupP0[RIGHT].resize(m_nbrFaceFlxPntsP0);
  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPntsP0; ++iFlx)
  {
    m_cellStatesFlxPntBackupP0[LEFT][iFlx].resize(m_nbrEqs); 
    m_cellStatesFlxPntBackupP0[RIGHT][iFlx].resize(m_nbrEqs);
    m_flxPntRiemannFluxBackupP0[LEFT][iFlx].resize(m_nbrEqs); 
    m_flxPntRiemannFluxBackupP0[RIGHT][iFlx].resize(m_nbrEqs);
  }


  m_extrapolatedFluxesBackupP0.resize(m_flxPntsLocalCoordsP0->size());
  m_contFlxBackupP0.resize(1);
  
  for (CFuint iFlx = 0; iFlx < m_flxPntsLocalCoordsP0->size(); ++iFlx)
  {
    m_extrapolatedFluxesBackupP0[iFlx].resize(m_nbrEqs);
  }

  for (CFuint iSolPnt = 0; iSolPnt < 1; ++iSolPnt)
  {
    m_contFlxBackupP0[iSolPnt].resize(m_dim+m_ndimplus);

    for (CFuint iDim = 0; iDim < m_dim+m_ndimplus; ++iDim)
    {
      m_contFlxBackupP0[iSolPnt][iDim].resize(m_nbrEqs);
    }
  }

  m_states_P0.resize(2);

  m_states_P0[LEFT] = new std::vector<Framework::State*>();
  m_states_P0[RIGHT] = new std::vector<Framework::State*>();

  m_states_P0[LEFT]->push_back(new State());
  m_states_P0[RIGHT]->push_back(new State());

  /*m_cellStatesP0 = new std::vector<Framework::State*>();  // allocate

  // later, to push:
  m_cellStatesP0->push_back(new State());  */

  m_PertcellStatesP0 = new std::vector<Framework::State*>();  // allocate

    // start with one (null) slot
    m_PertcellStatesP0->resize(1);

    // prepare dummy coords again
    RealVector dummyCoord;
    dummyCoord.resize(m_dim);
    dummyCoord = 0.0;
  
    // overwrite the nullptr with a real State*
    (*m_PertcellStatesP0)[0] = new State();
    (*m_PertcellStatesP0)[0]->setLocalID(0);
    (*m_PertcellStatesP0)[0]->setSpaceCoordinates( new Node(dummyCoord,false) );

}

//////////////////////////////////////////////////////////////////////////////

void ConvRHSJacobFluxReconstructionBlending::unsetup()
{
  CFAUTOTRACE;
  
  // unsetup parent class
  ConvRHSFluxReconstructionBlending::unsetup();

  deletePtr( (*m_PertcellStatesP0)[0] );

  m_PertcellStatesP0->clear();

  delete m_PertcellStatesP0;
  m_PertcellStatesP0 = nullptr;
}


//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD

