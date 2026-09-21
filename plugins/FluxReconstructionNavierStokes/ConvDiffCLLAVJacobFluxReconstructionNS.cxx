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

#include "FluxReconstructionNavierStokes/ConvDiffCLLAVJacobFluxReconstructionNS.hh"
#include "FluxReconstructionNavierStokes/FluxReconstructionNavierStokes.hh"
#include "FluxReconstructionMethod/FluxReconstructionElementData.hh"
#include "NavierStokes/NavierStokesVarSet.hh"

#include "NavierStokes/EulerTerm.hh"
#include "NavierStokes/EulerVarSet.hh"

#include "NavierStokes/Euler2DVarSet.hh"


//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider< ConvDiffCLLAVJacobFluxReconstructionNS,
		       FluxReconstructionSolverData,
		       FluxReconstructionNavierStokesModule >
convDiffCLLAVRHSJacobNSFluxReconstructionProvider("ConvDiffCLLAVRHSJacobNS");
  
//////////////////////////////////////////////////////////////////////////////
  
ConvDiffCLLAVJacobFluxReconstructionNS::ConvDiffCLLAVJacobFluxReconstructionNS(const std::string& name) :
  ConvDiffCLLAVJacobFluxReconstruction(name),
  m_tempGradTerm(),
  m_tempGradTermL(),
  m_tempGradTermR(),
  m_diffusiveVarSetNS(CFNULL),
  m_tempStatesL(),
  m_tempStatesR(),
  m_tempStatesL2(),
  m_tempStatesR2(),
  m_tempStatesCell(),
  m_eulerVarSet(CFNULL),
  m_msEulerTerm(CFNULL),
  m_nbrSpecies(),
  m_pData(),
  m_pData2(),
  m_eulerVarSet2(CFNULL),
  m_tempGradTermJacob(),
  m_tempStatesJacob(),
  m_tempGradTermJacob2(),
  m_tempStatesJacob2(),
  m_unpertGradVars(),
  m_pertGradVars()
{
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffCLLAVJacobFluxReconstructionNS::configure ( Config::ConfigArgs& args )
{
  ConvDiffCLLAVJacobFluxReconstruction::configure(args);
} 
  
//////////////////////////////////////////////////////////////////////////////

void ConvDiffCLLAVJacobFluxReconstructionNS::computeWaveSpeedUpdates(vector< CFreal >& waveSpeedUpd)
{
  // compute the wave speed updates for the neighbouring cells
  cf_assert(waveSpeedUpd.size() == 2);
  
  // here convective and artificial parts are added!
  ConvDiffCLLAVJacobFluxReconstruction::computeWaveSpeedUpdates(waveSpeedUpd);
          
  // now add diffusive part
  CFreal visc = 1.0;

  const CFreal dynVisc = m_diffusiveVarSetNS->getCurrDynViscosity();
  
  const CFreal factorPr = min(m_diffusiveVarSetNS->getModel().getPrandtl(),1.0);
  cf_assert(factorPr>0.0);
  
  for (CFuint iSide = 0; iSide < 2; ++iSide)
  {
    //for (CFuint iFlx = 0; iFlx < m_cellFlx[iSide].size(); ++iFlx)
    for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
    {
      const CFreal jacobXJacobXIntCoef = m_faceJacobVecAbsSizeFlxPnts[iFlx]*
                                 m_faceJacobVecAbsSizeFlxPnts[iFlx]*
                                   (*m_faceIntegrationCoefs)[iFlx]*
                                   m_cflConvDiffRatio;
      const CFreal rho = m_diffusiveVarSetNS->getDensity(*(m_cellStatesFlxPnt[iSide][iFlx]));
      visc = dynVisc/rho/factorPr;
      
      // transform update states to physical data to calculate eigenvalues
      waveSpeedUpd[iSide] += visc*jacobXJacobXIntCoef/m_cellVolume[iSide];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffCLLAVJacobFluxReconstructionNS::computeInterfaceFlxCorrection()
{
  ConvDiffLLAVJacobFluxReconstruction::computeInterfaceFlxCorrection();
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffCLLAVJacobFluxReconstructionNS::computeRiemannFluxJacobianNum(const CFreal resFactor)
{
  ConvDiffLLAVJacobFluxReconstruction::computeRiemannFluxJacobianNum(resFactor);
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffCLLAVJacobFluxReconstructionNS::computeRiemannFluxToGradJacobianNum(const CFreal resFactor)
{
  ConvDiffLLAVJacobFluxReconstruction::computeRiemannFluxToGradJacobianNum(resFactor);
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffCLLAVJacobFluxReconstructionNS::computeGradVarsToStateJacobianNum()
{
  CFLog(VERBOSE, "NS computeGradVarsToStateJacobianNum\n");

  for (m_pertSide = 0; m_pertSide < 2; ++m_pertSide)
  { 
    for (m_pertSol = 0; m_pertSol < m_nbrSolPnts; ++m_pertSol)
    {
      m_tempStatesJacob[0] = (*m_states[m_pertSide])[m_pertSol];
    
      m_diffusiveVarSetNS->setGradientVars(m_tempStatesJacob,m_tempGradTermJacob,1);
      
      for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
      {
        m_unpertGradVars[iEq] = m_tempGradTermJacob(iEq,0);
      }

      // dereference state
      State& pertState = *(*m_states[m_pertSide])[m_pertSol];

      for (m_pertVar = 0; m_pertVar < m_nbrEqs; ++m_pertVar)
      {
        // perturb physical variable in state
        m_numJacob->perturb(m_pertVar,pertState[m_pertVar]);
        
        m_tempStatesJacob2[0] = &pertState;
  
        m_diffusiveVarSetNS->setGradientVars(m_tempStatesJacob2,m_tempGradTermJacob2,1);
        
        for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
        {
          m_pertGradVars[iEq] = m_tempGradTermJacob2(iEq,0);
        }
        
        m_numJacob->computeDerivative(m_unpertGradVars,m_pertGradVars,m_gradVarsToStateJacobian[m_pertSide][m_pertSol][m_pertVar]);

        // restore physical variable in state
        m_numJacob->restore(pertState[m_pertVar]);
      }
    }
  }
  ////@TODO find a better way to to this
  
  // get current iteration
  const CFuint iter = SubSystemStatusStack::getActive()->getNbIter();
  
  // get the face start indexes
  vector< CFuint >& innerFacesStartIdxs = getMethodData().getInnerFacesStartIdxs();
  
  if (iter == 1 && m_faceID == innerFacesStartIdxs[0])
  {
      for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
      {
        (*(*m_states[0])[0])[iEq] += 0.1;
      }
      
      m_tempStatesJacob[0] = (*m_states[0])[0];
    
      m_diffusiveVarSetNS->setGradientVars(m_tempStatesJacob,m_tempGradTermJacob,1);
      
      for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
      {
        m_unpertGradVars[iEq] = m_tempGradTermJacob(iEq,0);
      }

      // dereference state
      State& pertState = *(*m_states[0])[0];

      for (m_pertVar = 0; m_pertVar < m_nbrEqs; ++m_pertVar)
      {
        // perturb physical variable in state
        m_numJacob->perturb(m_pertVar,pertState[m_pertVar]);
        
        m_tempStatesJacob2[0] = &pertState;
  
        m_diffusiveVarSetNS->setGradientVars(m_tempStatesJacob2,m_tempGradTermJacob2,1);
        
        for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
        {
          m_pertGradVars[iEq] = m_tempGradTermJacob2(iEq,0);
        }
        
        m_numJacob->computeDerivative(m_unpertGradVars,m_pertGradVars,m_tempFlux);
        
        m_varToGradVarDep[m_pertVar].resize(0);
      
        for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
        {
          if (m_tempFlux[iEq] != 0.0) m_varToGradVarDep[m_pertVar].push_back(iEq);
        }
        
        m_nbrVarToGradVarDep[m_pertVar] = m_varToGradVarDep[m_pertVar].size();       

        // restore physical variable in state
        m_numJacob->restore(pertState[m_pertVar]);
      }
        
      for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
      {
        (*(*m_states[0])[0])[iEq] -= 0.1;
      }
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffCLLAVJacobFluxReconstructionNS::computeBndGradTerms(RealMatrix& gradTerm, RealMatrix& ghostGradTerm)
{ 
  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    m_tempStatesL2[iFlx] = m_cellStatesFlxPnt[0][iFlx]->getData();
    m_tempStatesR2[iFlx] = m_flxPntGhostSol[iFlx]->getData();
  }
  
  m_diffusiveVarSetNS->setGradientVars(m_tempStatesL2,gradTerm,m_nbrFaceFlxPnts);
  m_diffusiveVarSetNS->setGradientVars(m_tempStatesR2,ghostGradTerm,m_nbrFaceFlxPnts);
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffCLLAVJacobFluxReconstructionNS::computeFlxPntGradTerm(const CFuint side, RealMatrix& gradTerm)
{
  vector< RealVector* > tempStates;
  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    tempStates.push_back(m_cellStatesFlxPnt[side][iFlx]->getData());
  }

  m_diffusiveVarSetNS->setGradientVars(tempStates,gradTerm,m_nbrFaceFlxPnts);
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffCLLAVJacobFluxReconstructionNS::computeCellGradTerm(RealMatrix& gradTerm)
{   
  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {
    m_tempStatesCell[iSol] = (*m_cellStates)[iSol]->getData();
  }
  
  m_diffusiveVarSetNS->setGradientVars(m_tempStatesCell,gradTerm,m_nbrSolPnts);
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffCLLAVJacobFluxReconstructionNS::computeFaceGradTerms(RealMatrix& gradTermL, RealMatrix& gradTermR)
{
  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    m_tempStatesL2[iFlx] = m_cellStatesFlxPnt[LEFT][iFlx]->getData();
    m_tempStatesR2[iFlx] = m_cellStatesFlxPnt[RIGHT][iFlx]->getData();
  }
  
  m_diffusiveVarSetNS->setGradientVars(m_tempStatesL2,gradTermL,m_nbrFaceFlxPnts);
  m_diffusiveVarSetNS->setGradientVars(m_tempStatesR2,gradTermR,m_nbrFaceFlxPnts);
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffCLLAVJacobFluxReconstructionNS::prepareFluxComputation()
{
  const bool isPerturb = this->getMethodData().isPerturb();
  const CFuint iPerturbVar = this->getMethodData().iPerturbVar();

  m_diffusiveVarSetNS->setComposition(m_avgSol, isPerturb, iPerturbVar);
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffCLLAVJacobFluxReconstructionNS::computeGradients()
{
  ConvDiffLLAVJacobFluxReconstruction::computeGradients();
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffCLLAVJacobFluxReconstructionNS::computeGradientFaceCorrections()
{
  ConvDiffLLAVJacobFluxReconstruction::computeGradientFaceCorrections();
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffCLLAVJacobFluxReconstructionNS::computeSmoothness()
{ 
  CFreal sNum = 0.0;
  
  CFreal sDenom = 0.0;
  
  // get datahandle
  DataHandle< CFreal > monPhysVar = socket_monPhysVar.getDataHandle();
  
  const bool RhoivtTv = getMethodData().getUpdateVarStr() == "RhoivtTv";
  
  if (RhoivtTv && m_monitoredPhysVar < m_pData.size())
  {
    for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
    {
      CFreal stateP = 0.0;
      CFreal diffStatesPPMinOne = 0.0;

      m_eulerVarSet2->computePhysicalData(*((*m_cellStates)[iSol]),m_pData);
      m_eulerVarSet2->computePhysicalData(m_statesPMinOne[iSol],m_pData2);

      stateP = m_pData[m_monitoredPhysVar];
      diffStatesPPMinOne = stateP - m_pData2[m_monitoredPhysVar];
    
      monPhysVar[(((*m_cellStates)[iSol]))->getLocalID()] = stateP;

      sNum += diffStatesPPMinOne*diffStatesPPMinOne;
      sDenom += stateP*stateP;
    }
  }
  else if (!RhoivtTv && m_monitoredPhysVar < m_pData.size())
  {
    for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
    {
      CFreal stateP = 0.0;
      CFreal diffStatesPPMinOne = 0.0;

      m_eulerVarSet->computePhysicalData(*((*m_cellStates)[iSol]),m_pData);
      m_eulerVarSet->computePhysicalData(m_statesPMinOne[iSol],m_pData2);

      stateP = m_pData[m_monitoredPhysVar];
      diffStatesPPMinOne = stateP - m_pData2[m_monitoredPhysVar];
    
      monPhysVar[(((*m_cellStates)[iSol]))->getLocalID()] = stateP;

      sNum += diffStatesPPMinOne*diffStatesPPMinOne;
      sDenom += stateP*stateP;
    }
  }
  else
  {
    for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
    {
      CFreal stateP = 0.0;
      CFreal diffStatesPPMinOne = 0.0;

      stateP = (*((*m_cellStates)[iSol]))[m_monitoredVar];
      diffStatesPPMinOne = stateP - m_statesPMinOne[iSol][m_monitoredVar];

      sNum += diffStatesPPMinOne*diffStatesPPMinOne;
      sDenom += stateP*stateP;
    }
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

void ConvDiffCLLAVJacobFluxReconstructionNS::computeSmoothness(const CFuint side)
{ 
  CFreal sNum = 0.0;
  
  CFreal sDenom = 0.0;
  
  // get datahandle
  DataHandle< CFreal > monPhysVar = socket_monPhysVar.getDataHandle();
  
  const bool RhoivtTv = getMethodData().getUpdateVarStr() == "RhoivtTv";
  
  if (RhoivtTv && m_monitoredPhysVar < m_pData.size())
  {
    for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
    {
      CFreal stateP = 0.0;
      CFreal diffStatesPPMinOne = 0.0;

      m_eulerVarSet2->computePhysicalData(*((*m_states[side])[iSol]),m_pData);
      m_eulerVarSet2->computePhysicalData(m_statesPMinOne[iSol],m_pData2);

      stateP = m_pData[m_monitoredPhysVar];
      diffStatesPPMinOne = stateP - m_pData2[m_monitoredPhysVar];

      sNum += diffStatesPPMinOne*diffStatesPPMinOne;
      sDenom += stateP*stateP;
    }
  }
  else if (!RhoivtTv && m_monitoredPhysVar < m_pData.size())
  {
    for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
    {
      CFreal stateP = 0.0;
      CFreal diffStatesPPMinOne = 0.0;

      m_eulerVarSet->computePhysicalData(*((*m_states[side])[iSol]),m_pData);
      m_eulerVarSet->computePhysicalData(m_statesPMinOne[iSol],m_pData2);

      stateP = m_pData[m_monitoredPhysVar];
      diffStatesPPMinOne = stateP - m_pData2[m_monitoredPhysVar];

      sNum += diffStatesPPMinOne*diffStatesPPMinOne;
      sDenom += stateP*stateP;
    }
  }
  else
  {
    for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
    {
      CFreal stateP = 0.0;
      CFreal diffStatesPPMinOne = 0.0;

      stateP = (*((*m_cellStates)[iSol]))[m_monitoredVar];
      diffStatesPPMinOne = stateP - m_statesPMinOne[iSol][m_monitoredVar];

      sNum += diffStatesPPMinOne*diffStatesPPMinOne;
      sDenom += stateP*stateP;
    }
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

void ConvDiffCLLAVJacobFluxReconstructionNS::computeEpsilon0()
{ 
  // compute a cell average characteristic flow speed. Note that a straight average is used, not a weighted one, maybe change this
  CFreal wavespeed = 0.0;
  
  const bool RhoivtTv = getMethodData().getUpdateVarStr() == "RhoivtTv";
  
  if (RhoivtTv)
  {
    for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
    {
      m_eulerVarSet2->computePhysicalData(*((*m_cellStates)[iSol]),m_pData);
      
      wavespeed += m_pData[EulerTerm::V] + m_pData[EulerTerm::A];
      cf_assert(m_pData[EulerTerm::V] > 0.0 && m_pData[EulerTerm::A] > 0.0);
    }
  }
  else
  {
    for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
    {
      m_eulerVarSet->computePhysicalData(*((*m_cellStates)[iSol]),m_pData);

      wavespeed += m_pData[EulerTerm::V] + m_pData[EulerTerm::A];
      cf_assert(m_pData[EulerTerm::V] > 0.0 && m_pData[EulerTerm::A] > 0.0);
    }
  }
  
  cf_assert(wavespeed > 0.0);
  
  wavespeed /= m_nbrSolPnts;
    
  const CFreal peclet = computePeclet();
  
  const CFreal oneOverDim = 1./m_dim;
  
  // get the cell volumes
  DataHandle< CFreal > cellVolumes = socket_cellVolumes.getDataHandle();
  
  const CFreal h = pow(cellVolumes[m_cell->getID()],oneOverDim);
  //if(m_cell->getID()==1) CFLog(INFO, "wvspd: " << h*wavespeed << ", P: " << peclet << ", dx: " << m_subcellRes << "\n");
  m_epsilon0 = max(h*wavespeed*(2.0/peclet - m_subcellRes/peclet),0.0);
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffCLLAVJacobFluxReconstructionNS::computeEpsilon0(const CFuint side)
{ 
  // compute a cell average characteristic flow speed. Note that a straight average is used, not a weighted one, maybe change this
  CFreal wavespeed = 0.0;
  
  const bool RhoivtTv = getMethodData().getUpdateVarStr() == "RhoivtTv";
  
  if (RhoivtTv)
  {
    for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
    {
      m_eulerVarSet2->computePhysicalData(*((*m_states[side])[iSol]),m_pData);
      
      wavespeed += m_pData[EulerTerm::V] + m_pData[EulerTerm::A];
      cf_assert(m_pData[EulerTerm::V] > 0.0 && m_pData[EulerTerm::A] > 0.0);
    }
  }
  else
  {
    for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
    {
      m_eulerVarSet->computePhysicalData(*((*m_states[side])[iSol]),m_pData);

      wavespeed += m_pData[EulerTerm::V] + m_pData[EulerTerm::A];
      cf_assert(m_pData[EulerTerm::V] > 0.0 && m_pData[EulerTerm::A] > 0.0);
    }
  }
  
  cf_assert(wavespeed > 0.0);
  
  wavespeed /= m_nbrSolPnts;
    
  const CFreal peclet = computePeclet();
  
  const CFreal oneOverDim = 1./m_dim;
  
  // get the cell volumes
  DataHandle< CFreal > cellVolumes = socket_cellVolumes.getDataHandle();
  
  const CFreal h = pow(cellVolumes[m_cells[side]->getID()],oneOverDim);
  
  m_epsilon0 = max(h*wavespeed*(2.0/peclet - m_subcellRes/peclet),0.0);
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffCLLAVJacobFluxReconstructionNS::computeEpsToStateJacobianAna()
{
  CFLog(VERBOSE, "NS computeEpsToStateJacobianAna\n");
    
  for (m_pertSide = 0; m_pertSide < 2; ++m_pertSide)
  {
    for (m_pertSol = 0; m_pertSol < m_nbrSolPnts; ++m_pertSol)
    {
      m_epsJacobian[m_pertSide][m_pertSol] = 0.0;
    }
  }
  
  const bool RhoivtTv = getMethodData().getUpdateVarStr() == "RhoivtTv";
  
  // get the cell volumes
  DataHandle< CFreal > cellVolumes = socket_cellVolumes.getDataHandle();
  
  const CFreal oneOverDim = 1./m_dim;
  
  const CFreal peclet = computePeclet();
  
  for (m_pertSide = 0; m_pertSide < 2; ++m_pertSide)
  { 
    CFreal h_f_S;
    
    const CFreal h = pow(cellVolumes[m_cells[m_pertSide]->getID()],oneOverDim);
    
    const CFreal h_f = h*(2.0/peclet - m_subcellRes/peclet);
    
    computeProjStates(m_statesPMinOne, m_pertSide);
    
    computeSmoothness(m_pertSide);
    
    CFreal sBefore;
    
    if (m_s < m_s0 - m_kappa)
    {
      h_f_S = 0.0;
      
      sBefore = 0.0;
    }
    else if (m_s > m_s0 + m_kappa)
    { 
      h_f_S = h_f/m_nbrSolPnts;
      
      sBefore = 1.0;
    }
    else
    {   
      sBefore = 0.5*(1.0 + sin(0.5*MathTools::MathConsts::CFrealPi()*(m_s-m_s0)/m_kappa));
              
      h_f_S = h_f/m_nbrSolPnts*sBefore;
    }
    
    CFreal wallFactor = 1.0;
    
    if (m_useWallCutOff)
    {
      // Get the wall distance
      DataHandle< CFreal > wallDist = socket_wallDistance.getDataHandle();
  
      CFreal centroidDistance = 0.0;
      
      for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
      {
        const CFuint stateID = (*m_states[m_pertSide])[iSol]->getLocalID();
        centroidDistance += wallDist[stateID];
      }
    
      centroidDistance /= m_nbrSolPnts;
    
      if (centroidDistance < m_wallCutOff) 
      {
        if (centroidDistance < 0.5*m_wallCutOff)
        {
          wallFactor = 0.0; 
        }
        else
        {
          wallFactor = 0.5*(1.0 + sin(0.5*MathTools::MathConsts::CFrealPi()*(centroidDistance-0.75*m_wallCutOff)/(0.25*m_wallCutOff)));
        }
      }
    }
    
    sBefore *= wallFactor;
      
    for (m_pertSol = 0; m_pertSol < m_nbrSolPnts; ++m_pertSol)
    {
      // dereference state
      State& pertState = *(*m_states[m_pertSide])[m_pertSol];
      
      if (RhoivtTv)
      {
        m_eulerVarSet2->computePhysicalData(*((*m_states[m_pertSide])[m_pertSol]),m_pData);
      }
      else
      {
        m_eulerVarSet->computePhysicalData(*((*m_states[m_pertSide])[m_pertSol]),m_pData);
      }
        
      const CFreal lambda = m_pData[EulerTerm::V] + m_pData[EulerTerm::A];
      cf_assert(m_pData[EulerTerm::V] > 0.0 && m_pData[EulerTerm::A] > 0.0);
          
      const CFreal h_f_lambda = h_f * lambda;
      
      for (m_pertVar = 0; m_pertVar < m_nbrEqs; ++m_pertVar)
      {
        const CFreal uBefore = pertState[m_pertVar];
        
        // perturb physical variable in state
        m_numJacob->perturb(m_pertVar,pertState[m_pertVar]);
        
        const CFreal eps = pertState[m_pertVar] - uBefore;
        
        if (RhoivtTv)
        {
          m_eulerVarSet2->computePhysicalData(*((*m_states[m_pertSide])[m_pertSol]),m_pData);
        }
        else
        {
          m_eulerVarSet->computePhysicalData(*((*m_states[m_pertSide])[m_pertSol]),m_pData);
        }
        
        const CFreal lambdaPert = m_pData[EulerTerm::V] + m_pData[EulerTerm::A];
          
        computeProjStates(m_statesPMinOne, m_pertSide);
    
        computeSmoothness(m_pertSide);
        
        CFreal sPert;
        
        if (m_s < m_s0 - m_kappa)
        { 
          sPert = 0.0;
        }
        else if (m_s > m_s0 + m_kappa)
        { 
          sPert = 1.0;
        }
        else
        {   
          sPert = 0.5*(1.0 + sin(0.5*MathTools::MathConsts::CFrealPi()*(m_s-m_s0)/m_kappa));
        }
        
        sPert *= wallFactor;

        // restore physical variable in state
        m_numJacob->restore(pertState[m_pertVar]);  
          
        m_epsJacobian[m_pertSide][m_pertSol][m_pertVar] += h_f_S*(lambdaPert - lambda)/eps;
          
        m_epsJacobian[m_pertSide][m_pertSol][m_pertVar] += h_f_lambda*(sPert - sBefore)/eps;
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffCLLAVJacobFluxReconstructionNS::setup()
{
  CFAUTOTRACE;
  
  // setup parent class
  ConvDiffCLLAVJacobFluxReconstruction::setup();
  
  // get the diffusive varset
  m_diffusiveVarSetNS = (getMethodData().getDiffusiveVar()).d_castTo< NavierStokesVarSet >();
  
  m_updateToSolutionVecTrans = getMethodData().getUpdateToSolutionVecTrans();
  
  m_updateToSolutionVecTrans->setup(2);
  
  m_tempGradTerm.resize(m_nbrEqs,m_nbrSolPnts);
  m_tempGradTermJacob.resize(m_nbrEqs,1);
  m_tempGradTermJacob2.resize(m_nbrEqs,1);

  
  m_tempGradTermL.resize(m_nbrEqs,m_nbrFaceFlxPnts);
  m_tempGradTermR.resize(m_nbrEqs,m_nbrFaceFlxPnts);
  
  m_tempStatesL.resize(m_nbrFaceFlxPnts);
  m_tempStatesR.resize(m_nbrFaceFlxPnts);
  m_tempStatesL2.resize(m_nbrFaceFlxPnts);
  m_tempStatesR2.resize(m_nbrFaceFlxPnts);
  m_tempStatesJacob.resize(1);
  m_tempStatesJacob2.resize(1);
  m_unpertGradVars.resize(m_nbrEqs);
  m_pertGradVars.resize(m_nbrEqs);
  
  m_tempStatesCell.resize(m_nbrSolPnts);
  
  const bool RhoivtTv = getMethodData().getUpdateVarStr() == "RhoivtTv";
  
  if(!RhoivtTv)
  {
    // get Euler varset
    m_eulerVarSet = getMethodData().getUpdateVar().d_castTo<EulerVarSet>();
    m_eulerVarSet->getModel()->resizePhysicalData(m_pData);
    m_eulerVarSet->getModel()->resizePhysicalData(m_pData2);
  } 
  else
  {
    m_eulerVarSet2 = getMethodData().getUpdateVar().d_castTo< MultiScalarVarSet< Euler2DVarSet > >();  
    
    m_msEulerTerm = PhysicalModelStack::getActive()-> getImplementor()->getConvectiveTerm().d_castTo< MultiScalarTerm< EulerTerm > >();
    if (m_msEulerTerm.isNull())
    {
      throw Common::ShouldNotBeHereException (FromHere(),"Update variable set is not MultiScalar EulerTerm in BCNoSlipWallrvt!");
    }
  
    m_nbrSpecies = m_msEulerTerm->getNbScalarVars(0);
    
    m_eulerVarSet2->getModel()->resizePhysicalData(m_pData);
    m_eulerVarSet2->getModel()->resizePhysicalData(m_pData2);
  } 
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD
