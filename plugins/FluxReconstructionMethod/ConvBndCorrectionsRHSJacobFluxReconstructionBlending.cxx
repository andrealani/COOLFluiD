// Copyright (C) 2016 KU Leuven, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/BlockAccumulator.hh"
#include "Framework/LSSMatrix.hh"
#include "Framework/MethodCommandProvider.hh"

#include "Framework/MeshData.hh"
#include "Framework/BaseTerm.hh"

#include "MathTools/MathFunctions.hh"

#include "FluxReconstructionMethod/ConvBndCorrectionsRHSJacobFluxReconstructionBlending.hh"
#include "FluxReconstructionMethod/FluxReconstruction.hh"
#include "FluxReconstructionMethod/FluxReconstructionElementData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider< ConvBndCorrectionsRHSJacobFluxReconstructionBlending, FluxReconstructionSolverData, FluxReconstructionModule >
  ConvBndCorrectionsRHSJacobFluxReconstructionBlendingProvider("ConvBndCorrectionsRHSJacobBlending");

//////////////////////////////////////////////////////////////////////////////

ConvBndCorrectionsRHSJacobFluxReconstructionBlending::ConvBndCorrectionsRHSJacobFluxReconstructionBlending(const std::string& name) :
  ConvBndCorrectionsRHSFluxReconstructionBlending(name),
  m_lss(CFNULL),
  m_numJacob(CFNULL),
  m_acc(CFNULL),
  m_pertResUpdates(),
  m_resUpdates(),
  m_derivResUpdates(),
  m_pertCorrections(),
  m_pertSol(),
  m_pertVar(),
  m_solFlxDep(CFNULL),
  m_nbrFlxDep(),
  m_cellStatesFlxPntBackup(),
  m_influencedFlxPnt(),
  m_influencedFlxPnts(),
  m_flxPntRiemannFluxBackup(),
  elemShape()
{
}

//////////////////////////////////////////////////////////////////////////////

ConvBndCorrectionsRHSJacobFluxReconstructionBlending::~ConvBndCorrectionsRHSJacobFluxReconstructionBlending()
{
}

//////////////////////////////////////////////////////////////////////////////

void ConvBndCorrectionsRHSJacobFluxReconstructionBlending::configure ( Config::ConfigArgs& args )
{
  ConvBndCorrectionsRHSFluxReconstructionBlending::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void ConvBndCorrectionsRHSJacobFluxReconstructionBlending::executeOnTrs()
{
  CFAUTOTRACE;

  // get InnerCells TopologicalRegionSet
  SafePtr<TopologicalRegionSet> cellTrs = MeshDataStack::getActive()->getTrs("InnerCells");

  // get current bnd face TRS
  SafePtr<TopologicalRegionSet> faceTrs = getCurrentTRS();
  
  CFLog(VERBOSE,"ConvBndCorrectionsRHSJacobFluxReconstructionBlending::executeOnTRS: " << faceTrs->getName() << "\n");

  // get bndFacesStartIdxs from FluxReconstructionMethodData
  map< std::string , vector< vector< CFuint > > >&
    bndFacesStartIdxsPerTRS = getMethodData().getBndFacesStartIdxs();
  vector< vector< CFuint > > bndFacesStartIdxs = bndFacesStartIdxsPerTRS[faceTrs->getName()];

  // number of face orientations (should be the same for all TRs)
  cf_assert(bndFacesStartIdxs.size() != 0);
  CFuint nbOrients = bndFacesStartIdxs[0].size()-1;

  // number of TRs
  const CFuint nbTRs = faceTrs->getNbTRs();
  cf_assert(bndFacesStartIdxs.size() == nbTRs);

  // get the geodata of the face builder and set the TRSs
  FaceToCellGEBuilder::GeoData& geoData = m_faceBuilder->getDataGE();
  geoData.cellsTRS = cellTrs;
  geoData.facesTRS = faceTrs;
  geoData.isBoundary = true;
  
  // boolean telling whether there is a diffusive term
  const bool hasDiffTerm = getMethodData().hasDiffTerm() || getMethodData().hasArtificialViscosity();
  m_bcStateComputer->preProcess();
  // loop over TRs
  for (CFuint iTR = 0; iTR < nbTRs; ++iTR)
  {
    nbOrients = bndFacesStartIdxs[iTR].size()-1;
    // loop over different orientations
    for (m_orient = 0; m_orient < nbOrients; ++m_orient)
    {
      CFLog(VERBOSE,"m_orient: " << m_orient << "\n");
      
      // Reset the value of m_nbrFaceFlxPnts in case it is not the same for all faces (Prism)
      m_nbrFaceFlxPnts=(*m_faceFlxPntConn)[m_orient].size();
      
      // select the correct flx pnts on the face out of all cell flx pnts for the current orient
      for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
      {
        m_flxPntsLocalCoords[iFlx] = (*m_allCellFlxPnts)[(*m_faceFlxPntConn)[m_orient][iFlx]];
      }
      
      // start and stop index of the faces with this orientation
      const CFuint startFaceIdx = bndFacesStartIdxs[iTR][m_orient  ];
      const CFuint stopFaceIdx  = bndFacesStartIdxs[iTR][m_orient+1];

      // loop over faces with this orientation
      for (CFuint faceID = startFaceIdx; faceID < stopFaceIdx; ++faceID)
      {
        // build the face GeometricEntity
        geoData.idx = faceID;
        m_face = m_faceBuilder->buildGE();

        // get the neighbouring cell
        m_intCell = m_face->getNeighborGeo(0);
	
	// get the states in the neighbouring cell
        m_cellStates = m_intCell->getStates();
	
        CFLog(VERBOSE,"cellID: " << m_intCell->getID() << "\n");
	CFLog(VERBOSE,"coord state 0: " << (((*m_cellStates)[0])->getCoordinates()) << "\n");


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


        // if cell is parallel updatable or the gradients have to be computed, compute the necessary data
        if ((*m_cellStates)[0]->isParUpdatable() || hasDiffTerm)
        {  
	  // set the bnd face data
	  setBndFaceData(m_face->getID());//faceID
	  
	  // compute the perturbed states and ghost states in the flx pnts
          computeFlxPntStates();
	}
	
	// if the cell is parallel updatable, compute the flx correction
	if ((*m_cellStates)[0]->isParUpdatable())
	{
	  // compute FI-FD
          computeInterfaceFlxCorrection();
	  
          // compute the wave speed updates
          computeWaveSpeedUpdates(m_waveSpeedUpd, m_waveSpeedUpdP0);
      
          // update the wave speeds
          updateWaveSpeed();
       
	  // compute the correction -(FI-FD)divh of the bnd face for each sol pnt
          computeCorrection(m_corrections);
	  
	  // update the rhs
          updateRHS();
	}
	  
	// if there is a diffusive term, compute the gradients
        if (hasDiffTerm)
        {
          computeGradientBndFaceCorrections();
        }
        
        const CFuint iter = SubSystemStatusStack::getActive()->getNbIter();
    
        const CFuint iterFreeze = getMethodData().getFreezeJacobIter();
    
        const CFuint interval = iter - iterFreeze;
      
        if (!getMethodData().freezeJacob() || iter < iterFreeze || interval % getMethodData().getFreezeJacobInterval() == 0)
        {
	
	  // if the cell is parallel updatable, compute the contribution to the numerical jacobian
	  if ((*m_cellStates)[0]->isParUpdatable())
	  {
	    // compute the convective boundary flux correction contribution to the jacobian
	    computeJacobConvBndCorrection();
          }
        }
        
        // release the face
        m_faceBuilder->releaseGE();
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvBndCorrectionsRHSJacobFluxReconstructionBlending::computeJacobConvBndCorrection()
{
  // get residual factor
  const CFreal resFactor = getMethodData().getResFactor();

  // dereference accumulator
  BlockAccumulator& acc = *m_acc;

  // set block row and column indices
  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {
    acc.setRowColIndex(iSol,(*m_cellStates)[iSol]->getLocalID());
    
    // put the perturbed and unperturbed corrections in the correct format
    for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
    {
      m_resUpdates[m_nbrEqs*iSol+iVar] = m_corrections[iSol][iVar];
    }
  }
  
  // store backups of values that will be overridden
  storeBackups();

  // loop over the states in the internal cell to perturb the states
  for (m_pertSol = 0; m_pertSol < m_nbrSolPnts; ++m_pertSol)
  {
    // dereference state
    State& pertState = *(*m_cellStates)[m_pertSol];
    
    // Loop over flux points to determine which flx pnts are influenced by the pert
    m_influencedFlxPnts.resize(0);
    for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
    {
      // get current flx pnt idx
      const CFuint currFlxIdx = (*m_faceFlxPntConn)[m_orient][iFlxPnt];
    
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
      // perturb physical variable in HO state
      m_numJacob->perturb(m_pertVar,pertState[m_pertVar]);
      
      // Recompute P0 state from perturbed HO states (P0 = cell average depends on ALL sol pnts)
      for (CFuint i = 0; i < m_nbrSolPnts; ++i) 
      { 
        m_P0State[i][m_pertVar] = 0.0;
        for (CFuint j = 0; j < m_nbrSolPnts; ++j) 
        { 
          m_P0State[i][m_pertVar] += m_vdm(i,0) * m_vdmInv(0,j) * (*(*m_cellStates)[j])[m_pertVar];
        }
      }
      // Update P0 state pointer
      (*(*m_cellStatesP0)[0])[m_pertVar] = m_P0State[0][m_pertVar];
      
      // compute the perturbed states and ghost states in the flx pnts
      extrapolatePerturbedState();

      // compute the perturbed interface flx correction
      computePertInterfaceFlxCorrection();
      computePertCorrection(m_pertCorrections);
      
      // put the perturbed and unperturbed corrections in the correct format
      for (CFuint iState = 0; iState < m_nbrSolPnts; ++iState)
      {
        for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
	{
          m_pertResUpdates[m_nbrEqs*iState+iVar] = m_pertCorrections[iState][iVar];
        }
      }

      // add contributions to the Jacobian
      // compute the finite difference derivative of the face term
      m_numJacob->computeDerivative(m_pertResUpdates,m_resUpdates,m_derivResUpdates);

      // multiply residual update derivatives with residual factor
      m_derivResUpdates *= resFactor;
      
      // reset updated flags
      for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
      {
        m_solPntUpdated[iSol] = false;
      }
      
      for (CFuint iFlxPnt = 0; iFlxPnt < m_NbInfluencedFlxPnts; ++iFlxPnt)
      {   
        const CFuint flxIdx = (*m_faceFlxPntConn)[m_orient][ m_influencedFlxPnts[iFlxPnt] ];

        m_nbrSolDep = ((*m_flxSolDep)[flxIdx]).size();
        for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolDep; ++iSolPnt)
        {
          const CFuint solIdx = (*m_flxSolDep)[flxIdx][iSolPnt];
          
          if (!m_solPntUpdated[solIdx])
          {
            m_solPntUpdated[solIdx] = true;
            acc.addValues(solIdx,m_pertSol,m_pertVar,&m_derivResUpdates[m_nbrEqs*solIdx]);
          }
        }  
      } 

      // restore physical variable in state
      m_numJacob->restore(pertState[m_pertVar]);
      
      // restore the overridden values
      restoreFromBackups();
    }
  }
  
//   if (m_intCell->getID() == 56)
//   {
//   CFLog(VERBOSE,"accBndFace:\n");
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

void ConvBndCorrectionsRHSJacobFluxReconstructionBlending::extrapolatePerturbedState()
{       
  for (CFuint iFlxPnt = 0; iFlxPnt < m_NbInfluencedFlxPnts; ++iFlxPnt)
  {   
    // get current flx pnt idx
    const CFuint currFlxIdx = (*m_faceFlxPntConn)[m_orient][ m_influencedFlxPnts[iFlxPnt] ];
    // reset the extrapolated states
    (*(m_cellStatesFlxPnt[ m_influencedFlxPnts[iFlxPnt] ]))[m_pertVar] = 0.0;
    
    // extrapolate the states to current flx pnt
    m_nbrSolDep = ((*m_flxSolDep)[currFlxIdx]).size();
    for (CFuint iSol = 0; iSol < m_nbrSolDep; ++iSol)
    {
      const CFuint solIdx = (*m_flxSolDep)[currFlxIdx][iSol];

    (*(m_cellStatesFlxPnt[ m_influencedFlxPnts[iFlxPnt] ]))[m_pertVar] += (*m_solPolyValsAtFlxPnts)[currFlxIdx][solIdx]*(*((*m_cellStates)[solIdx]))[m_pertVar];
    }
  }  
  // compute ghost states
  m_bcStateComputer->computeGhostStates(m_cellStatesFlxPnt,m_flxPntGhostSol,m_unitNormalFlxPnts,m_flxPntCoords);

  for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPntsP0; ++iFlxPnt)
  {   
    // get current flx pnt idx
    const CFuint currFlxIdx = (*m_faceFlxPntConnP0)[m_orient][iFlxPnt];
    // reset the extrapolated states
    (*(m_cellStatesFlxPntP0[ iFlxPnt ]))[m_pertVar] = 0.0;
    
    // extrapolate the states to current flx pnt
    m_nbrSolDep = ((*m_flxSolDepP0)[currFlxIdx]).size();
    for (CFuint iSol = 0; iSol < m_nbrSolDep; ++iSol)
    {
      const CFuint solIdx = (*m_flxSolDepP0)[currFlxIdx][iSol];

    (*(m_cellStatesFlxPntP0[ iFlxPnt ]))[m_pertVar] += (*m_solPolyValsAtFlxPntsP0)[currFlxIdx][solIdx]*((((m_P0State))[0]))[m_pertVar];
    }
  }  
  // compute ghost states
  m_bcStateComputer->computeGhostStates(m_cellStatesFlxPntP0,m_flxPntGhostSolP0,m_unitNormalFlxPntsP0,m_flxPntCoordsP0);
}

//////////////////////////////////////////////////////////////////////////////

void ConvBndCorrectionsRHSJacobFluxReconstructionBlending::computePertInterfaceFlxCorrection()
{ 
  for (CFuint iFlxPnt = 0; iFlxPnt < m_NbInfluencedFlxPnts; ++iFlxPnt){   
    // compute the riemann flux in the flx pnts
    m_flxPntRiemannFlux[ m_influencedFlxPnts[iFlxPnt] ] = m_riemannFluxComputer->computeFlux(*(m_cellStatesFlxPnt[ m_influencedFlxPnts[iFlxPnt] ]),
                          *(m_flxPntGhostSol[ m_influencedFlxPnts[iFlxPnt] ]),
                          m_unitNormalFlxPnts[ m_influencedFlxPnts[iFlxPnt] ]);
      
    // store the local Riemann flux, scaled with geometrical Jacobian
    m_flxPntRiemannFlux[ m_influencedFlxPnts[iFlxPnt] ] *= m_faceJacobVecSizeFlxPnts[ m_influencedFlxPnts[iFlxPnt] ];
  }

  for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPntsP0; ++iFlxPnt){   
    // compute the riemann flux in the flx pnts
    m_flxPntRiemannFluxP0[ iFlxPnt ] = m_riemannFluxComputer->computeFlux(*(m_cellStatesFlxPntP0[ iFlxPnt ]),
                          *(m_flxPntGhostSolP0[ iFlxPnt ]),
                          m_unitNormalFlxPntsP0[ iFlxPnt ]);
      
    // store the local Riemann flux, scaled with geometrical Jacobian
    m_flxPntRiemannFluxP0[ iFlxPnt ] *= m_faceJacobVecSizeFlxPntsP0[ iFlxPnt ];
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvBndCorrectionsRHSJacobFluxReconstructionBlending::computePertCorrection(vector< RealVector >& corrections)
{ 
  cf_assert(corrections.size() == m_nbrSolPnts); 

  // Get alpha (frozen from socket)
  DataHandle< CFreal > output = socket_alpha.getDataHandle();
  const CFreal alpha = output[((*m_cellStates)[0])->getLocalID()];

  // reset updated flags
  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {
    m_solPntUpdated[iSol] = false;
  }

  // reset corrections for all influenced solution points (following original structure)
  for (CFuint iFlxPnt = 0; iFlxPnt < m_NbInfluencedFlxPnts; ++iFlxPnt)
  {
    const CFuint flxIdx = (*m_faceFlxPntConn)[m_orient][m_influencedFlxPnts[iFlxPnt]];
    m_nbrSolDep = ((*m_flxSolDep)[flxIdx]).size();      
    for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolDep; ++iSolPnt)
    {
      const CFuint solIdx = (*m_flxSolDep)[flxIdx][iSolPnt];
      if (!m_solPntUpdated[solIdx])
      {
        m_solPntUpdated[solIdx] = true;
        corrections[solIdx] = 0.0;
      }
    }
  }  
  
  // STEP 1: Add high-order contributions (loop over flux points)
  for (CFuint iFlxPnt = 0; iFlxPnt < m_NbInfluencedFlxPnts; ++iFlxPnt)
  { 
    const CFuint flxIdx = (*m_faceFlxPntConn)[m_orient][ m_influencedFlxPnts[iFlxPnt] ];
    
    // the current correction factor previously computed
    // This is the high-order flux correction
    const RealVector& currentCorrFactor = m_flxPntRiemannFlux[ m_influencedFlxPnts[iFlxPnt] ];

    cf_assert(currentCorrFactor.size() == m_nbrEqs);
    m_nbrSolDep = ((*m_flxSolDep)[flxIdx]).size();
    // compute the term due to each flx pnt for all affected solution points
    for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolDep; ++iSolPnt)
    {
      const CFuint solIdx = (*m_flxSolDep)[flxIdx][iSolPnt];

      // divergence of the correctionfct
      const CFreal divh = m_corrFctDiv[solIdx][flxIdx];

      // Fill in the high-order corrections with blending factor (1-alpha)
      for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
      {  
        corrections[solIdx][iVar] -= currentCorrFactor[iVar] * divh * (1.0 - alpha); 
      }
    }
  }

  // STEP 2: Add P0 contribution (only once per solution point)
  // This applies to all solution points in the cell
  const CFuint flxIdx0 = (*m_faceFlxPntConn)[m_orient][0];
  const CFuint flxIdxP0 = (*m_faceFlxPntConnP0)[m_flxPntFaceConn[flxIdx0]][0]; 
  const RealVector& currentCorrFactorP0 = m_flxPntRiemannFluxP0[0];

  for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolPnts; ++iSolPnt)
  {
    const CFreal divhP0 = m_corrFctDivP0[0][flxIdxP0];
    for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
    {
      corrections[iSolPnt][iVar] -= currentCorrFactorP0[iVar] * divhP0 * alpha;// * (*m_cellAvgSolCoefs)[iSolPnt];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvBndCorrectionsRHSJacobFluxReconstructionBlending::storeBackups()
{
  for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
  {
    m_cellStatesFlxPntBackup[iFlxPnt] = *(m_cellStatesFlxPnt[iFlxPnt]);
    m_flxPntRiemannFluxBackup[iFlxPnt] = m_flxPntRiemannFlux[iFlxPnt];
  }

  for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPntsP0; ++iFlxPnt)
  {
    m_cellStatesFlxPntBackupP0[iFlxPnt] = *(m_cellStatesFlxPntP0[iFlxPnt]);
    m_flxPntRiemannFluxBackupP0[iFlxPnt] = m_flxPntRiemannFluxP0[iFlxPnt];
  }
  
  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {
    m_pertCorrections[iSol] = m_corrections[iSol];
    // Backup P0 states
    m_P0StateBackup[iSol] = m_P0State[iSol];
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvBndCorrectionsRHSJacobFluxReconstructionBlending::restoreFromBackups()
{
  for (CFuint iFlxPnt = 0; iFlxPnt < m_NbInfluencedFlxPnts; ++iFlxPnt)
  { 
  (*(m_cellStatesFlxPnt[ m_influencedFlxPnts[iFlxPnt] ]))[m_pertVar] = m_cellStatesFlxPntBackup[ m_influencedFlxPnts[iFlxPnt] ][m_pertVar];
  m_flxPntRiemannFlux[ m_influencedFlxPnts[iFlxPnt] ] = m_flxPntRiemannFluxBackup[ m_influencedFlxPnts[iFlxPnt] ];
  }

  for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPntsP0; ++iFlxPnt)
  { 
  (*(m_cellStatesFlxPntP0[ iFlxPnt ]))[m_pertVar] = m_cellStatesFlxPntBackupP0[ iFlxPnt ][m_pertVar];
  m_flxPntRiemannFluxP0[ iFlxPnt ] = m_flxPntRiemannFluxBackupP0[ iFlxPnt ];
  }
  
  // Restore P0 states
  for (CFuint i = 0; i < m_nbrSolPnts; ++i) 
  { 
    m_P0State[i][m_pertVar] = m_P0StateBackup[i][m_pertVar];
  }
  (*(*m_cellStatesP0)[0])[m_pertVar] = m_P0StateBackup[0][m_pertVar];
  
  // reset updated flags
  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {
    m_solPntUpdated[iSol] = false;
  }
  
  // loop over all influenced flux points
  for (CFuint iInfluencedFlx = 0; iInfluencedFlx < m_NbInfluencedFlxPnts; ++iInfluencedFlx)
  {
    const CFuint flxIdx = (*m_faceFlxPntConn)[m_orient][m_influencedFlxPnts[iInfluencedFlx]];
    m_nbrSolDep = ((*m_flxSolDep)[flxIdx]).size();
    for (CFuint iSol = 0; iSol < m_nbrSolDep; ++iSol)
    {
      const CFuint solIdx = (*m_flxSolDep)[flxIdx][iSol];
      if (!m_solPntUpdated[solIdx])
      {
        m_solPntUpdated[solIdx] = true;
        m_pertCorrections[solIdx] = m_corrections[solIdx];
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvBndCorrectionsRHSJacobFluxReconstructionBlending::setup()
{
  CFAUTOTRACE;

  // setup parent class
  ConvBndCorrectionsRHSFluxReconstructionBlending::setup();
  
  // get the linear system solver
  m_lss = getMethodData().getLinearSystemSolver()[0];

  // get the numerical Jacobian computer
  m_numJacob = getMethodData().getNumericalJacobian();

  m_numJacob_P0 = getMethodData().getNumericalJacobian();

  // create blockaccumulator
  m_acc.reset(m_lss->createBlockAccumulator(m_nbrSolPnts,m_nbrSolPnts,m_nbrEqs));
  
  // get the local FR data
  vector< FluxReconstructionElementData* >& frLocalData = getMethodData().getFRLocalData();
  
  //get element shape
  elemShape = frLocalData[0]->getShape();

  m_solFlxDep = frLocalData[0]->getSolPntFlxDependency();

  m_nbrFlxDep = ((*m_solFlxDep)[0]).size();
  
  // resize variables
  const CFuint nbrRes = m_nbrSolPnts*m_nbrEqs;
  m_pertResUpdates .resize(nbrRes);
  m_derivResUpdates.resize(nbrRes);
  m_resUpdates .resize(nbrRes);
  m_pertCorrections.resize(m_nbrSolPnts);
  m_cellStatesFlxPntBackup.resize(m_nbrFaceFlxPnts);
  m_flxPntRiemannFluxBackup.resize(m_nbrFaceFlxPnts);
  
  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    m_cellStatesFlxPntBackup[iFlx].resize(m_nbrEqs); 
    m_flxPntRiemannFluxBackup[iFlx].resize(m_nbrEqs);
  }
  
  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {
    m_pertCorrections[iSol].resize(m_nbrEqs);
  }

  m_cellStatesFlxPntBackupP0.resize(m_nbrFaceFlxPntsP0);
  m_flxPntRiemannFluxBackupP0.resize(m_nbrFaceFlxPntsP0);
  
  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPntsP0; ++iFlx)
  {
    m_cellStatesFlxPntBackupP0[iFlx].resize(m_nbrEqs); 
    m_flxPntRiemannFluxBackupP0[iFlx].resize(m_nbrEqs);
  }
  
  m_cellStatesP0 = new std::vector<Framework::State*>();  // allocate

  // later, to push:
  m_cellStatesP0->push_back(new State()); 
  
  // resize m_solPntUpdated
  m_solPntUpdated.resize(m_nbrSolPnts);
  
  // Initialize P0 state backup
  m_P0StateBackup.resize(m_nbrSolPnts);
  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {
    m_P0StateBackup[iSol].resize(m_nbrEqs);
  }
  
}

//////////////////////////////////////////////////////////////////////////////

void ConvBndCorrectionsRHSJacobFluxReconstructionBlending::unsetup()
{
  CFAUTOTRACE;

  // unsetup parent class
  ConvBndCorrectionsRHSFluxReconstructionBlending::unsetup();
}
//////////////////////////////////////////////////////////////////////////////

    } // namespace FluxReconstructionMethod

} // namespace COOLFluiD
