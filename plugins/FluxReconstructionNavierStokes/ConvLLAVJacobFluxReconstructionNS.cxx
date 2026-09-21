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

#include "FluxReconstructionNavierStokes/ConvLLAVJacobFluxReconstructionNS.hh"
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

MethodCommandProvider< ConvLLAVJacobFluxReconstructionNS,
		       FluxReconstructionSolverData,
		       FluxReconstructionNavierStokesModule >
convLLAVRHSJacobNSFluxReconstructionProvider("ConvLLAVRHSJacobNS");
  
//////////////////////////////////////////////////////////////////////////////
  
ConvLLAVJacobFluxReconstructionNS::ConvLLAVJacobFluxReconstructionNS(const std::string& name) :
  ConvLLAVJacobFluxReconstruction(name),
  m_eulerVarSet(CFNULL),
  m_msEulerTerm(CFNULL),
  m_nbrSpecies(),
  m_pData(),
  m_pData2(),
  m_eulerVarSet2(CFNULL)
{
}

//////////////////////////////////////////////////////////////////////////////

void ConvLLAVJacobFluxReconstructionNS::configure ( Config::ConfigArgs& args )
{
  ConvLLAVJacobFluxReconstruction::configure(args);
} 

//////////////////////////////////////////////////////////////////////////////

void ConvLLAVJacobFluxReconstructionNS::computeInterfaceFlxCorrection()
{
  ConvLLAVJacobFluxReconstruction::computeInterfaceFlxCorrection();
}

//////////////////////////////////////////////////////////////////////////////

void ConvLLAVJacobFluxReconstructionNS::computeRiemannFluxJacobianNum(const CFreal resFactor)
{
  ConvLLAVJacobFluxReconstruction::computeRiemannFluxJacobianNum(resFactor);
}

//////////////////////////////////////////////////////////////////////////////

void ConvLLAVJacobFluxReconstructionNS::prepareFluxComputation()
{
  //const bool isPerturb = this->getMethodData().isPerturb();
  //const CFuint iPerturbVar = this->getMethodData().iPerturbVar();

  //m_diffusiveVarSetNS->setComposition(m_avgSol, isPerturb, iPerturbVar);
}

//////////////////////////////////////////////////////////////////////////////

void ConvLLAVJacobFluxReconstructionNS::computeGradients()
{
  ConvLLAVJacobFluxReconstruction::computeGradients();
}

//////////////////////////////////////////////////////////////////////////////

void ConvLLAVJacobFluxReconstructionNS::computeGradientFaceCorrections()
{
  ConvLLAVJacobFluxReconstruction::computeGradientFaceCorrections();
}

//////////////////////////////////////////////////////////////////////////////

void ConvLLAVJacobFluxReconstructionNS::computeSmoothness()
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

void ConvLLAVJacobFluxReconstructionNS::computeSmoothness(const CFuint side)
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

void ConvLLAVJacobFluxReconstructionNS::computeEpsilon0()
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

void ConvLLAVJacobFluxReconstructionNS::computeEpsilon0(const CFuint side)
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

void ConvLLAVJacobFluxReconstructionNS::computeEpsToStateJacobianAna()
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

void ConvLLAVJacobFluxReconstructionNS::setup()
{
  CFAUTOTRACE;
  
  // setup parent class
  ConvLLAVJacobFluxReconstruction::setup();
  
  // get the diffusive varset
  //m_diffusiveVarSetNS = (getMethodData().getDiffusiveVar()).d_castTo< NavierStokesVarSet >();
  
  m_updateToSolutionVecTrans = getMethodData().getUpdateToSolutionVecTrans();
  
  m_updateToSolutionVecTrans->setup(2);
  
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
