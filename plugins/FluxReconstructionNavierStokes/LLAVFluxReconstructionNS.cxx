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

#include "NavierStokes/EulerTerm.hh"

#include "NavierStokes/NavierStokesVarSet.hh"
#include "NavierStokes/EulerVarSet.hh"

#include "FluxReconstructionNavierStokes/LLAVFluxReconstructionNS.hh"
#include "FluxReconstructionMethod/FluxReconstruction.hh"
#include "FluxReconstructionMethod/FluxReconstructionElementData.hh"

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

MethodCommandProvider< LLAVFluxReconstructionNS,
		       FluxReconstructionSolverData,
		       FluxReconstructionModule >
LLAVFluxReconstructionNSFluxReconstructionProvider("LLAVNS");

//////////////////////////////////////////////////////////////////////////////
  
LLAVFluxReconstructionNS::LLAVFluxReconstructionNS(const std::string& name) :
  LLAVFluxReconstruction(name),
  m_diffVarSet(CFNULL),
  m_gradsBackUp(),
  m_eulerVarSet(CFNULL),
  m_msEulerTerm(CFNULL),
  m_nbrSpecies(),
  m_pData(),
  m_pData2(),
  m_eulerVarSet2(CFNULL)
  {
  }

//////////////////////////////////////////////////////////////////////////////

void LLAVFluxReconstructionNS::configure ( Config::ConfigArgs& args )
{
  LLAVFluxReconstruction::configure(args);
}  

//////////////////////////////////////////////////////////////////////////////

void LLAVFluxReconstructionNS::setFaceData(CFuint faceID)
{
  LLAVFluxReconstruction::setFaceData(faceID);

//   m_diffVarSet = getMethodData().getDiffusiveVar();
// 
//   SafePtr< NavierStokesVarSet > navierStokesVarSet = m_diffVarSet.d_castTo< NavierStokesVarSet >();
// 
//   for (CFuint iSide = 0; iSide < 2; ++iSide)
//   {
//     vector< vector< RealVector* > > grads;
//     vector< RealVector* > tempStates;
//     vector< vector< RealVector > > temp;
//     temp.resize(m_nbrSolPnts);
//     grads.resize(m_nbrSolPnts);
// 
//     // make a back up of the grads
//     for (CFuint iState = 0; iState < m_nbrSolPnts; ++iState)
//     {
//       tempStates.push_back((*(m_states[iSide]))[iState]->getData());
// 
//       temp[iState] = *(m_cellGrads[iSide][iState]);
//       
//       grads[iState].resize(m_nbrEqs);
// 
//       for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
//       {
//         cf_assert(temp[iState].size() == m_nbrEqs);
//         grads[iState][iVar] = & (temp[iState][iVar]);
//       }
//     }
// 
//     navierStokesVarSet->setStateGradients(tempStates,grads,m_gradsBackUp[iSide],m_nbrSolPnts);
// 
//     // make store the correct gradients
//     for (CFuint iState = 0; iState < m_nbrSolPnts; ++iState)
//     {    
//       for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
//       { 
//         (*m_cellGrads[iSide][iState])[iVar] = *(m_gradsBackUp[iSide][iState][iVar]);
//       }
//     }
//   }
  if (getMethodData().hasDiffTerm() || getMethodData().getUpdateVarStr() == "Puvt")
  {
    // get the gradients datahandle
    DataHandle< vector< RealVector > > gradientsAV = socket_gradientsAV.getDataHandle();

    for (CFuint iSide = 0; iSide < 2; ++iSide)
    {
      for (CFuint iState = 0; iState < m_nbrSolPnts; ++iState)
      {
        const CFuint stateID = (*(m_states[iSide]))[iState]->getLocalID();
        m_cellGrads[iSide][iState] = &gradientsAV[stateID];
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void LLAVFluxReconstructionNS::setCellData()
{
  LLAVFluxReconstruction::setCellData();
  
//   m_diffVarSet = getMethodData().getDiffusiveVar();
//   
//   SafePtr< NavierStokesVarSet > navierStokesVarSet = m_diffVarSet.d_castTo< NavierStokesVarSet >();
//   
//   vector< vector< RealVector* > > grads;
//   vector< RealVector* > tempStates;
//   vector< vector< RealVector > > temp;
//   temp.resize(m_nbrSolPnts);
//   grads.resize(m_nbrSolPnts);
//   
//   // make a back up of the grads
//   for (CFuint iState = 0; iState < m_nbrSolPnts; ++iState)
//   {
//     tempStates.push_back((*m_cellStates)[iState]->getData());
//     
//     temp[iState] = *(m_cellGrads[0][iState]);
//       
//     grads[iState].resize(m_nbrEqs);
//     
//     for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
//     {
//       cf_assert(temp[iState].size() == m_nbrEqs);
//       grads[iState][iVar] = & (temp[iState][iVar]);
//     }
//   }
//   
//   navierStokesVarSet->setStateGradients(tempStates,grads,m_gradsBackUp[0],m_nbrSolPnts);
//     
//   // make store the correct gradients
//   for (CFuint iState = 0; iState < m_nbrSolPnts; ++iState)
//   {    
//     for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
//     { 
//       (*m_cellGrads[0][iState])[iVar] = *(m_gradsBackUp[0][iState][iVar]);
//     }
//   }
  
  if (getMethodData().hasDiffTerm() || getMethodData().getUpdateVarStr() == "Puvt")
  {
    // get the gradients datahandle
    DataHandle< vector< RealVector > > gradientsAV = socket_gradientsAV.getDataHandle();

    // get the grads in the current cell
    for (CFuint iState = 0; iState < m_nbrSolPnts; ++iState)
    {
      const CFuint stateID = (*(m_cellStates))[iState]->getLocalID();
      m_cellGrads[0][iState] = &gradientsAV[stateID];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

CFreal LLAVFluxReconstructionNS::computeViscCoef(RealVector* state)
{
  CFreal result = 1.0;
  
//   if (getMethodData().getUpdateVarStr() == "Cons")
//   {
//     result = 1./(*state)[0];
//   }
//   else if (getMethodData().getUpdateVarStr() == "Puvt")
//   {
//     const CFreal R = m_eulerVarSet->getModel()->getR();
//     result = (*state)[m_nbrEqs-1]/(*state)[0]*R;
//   }
//   else if (getMethodData().getUpdateVarStr() == "RhoivTv")
//   {
//     for (CFuint iSpecies = 0; iSpecies < m_nbrSpecies; ++iSpecies)
//     {
//       result += (*state)[iSpecies];
//     }
//     result = 1.0/result;
//   }
  
  return result;
}

//////////////////////////////////////////////////////////////////////////////

void LLAVFluxReconstructionNS::computeSmoothness()
{ 
  CFreal sNum = 0.0;
  
  CFreal sDenom = 0.0;
  
  // get datahandle
  DataHandle< CFreal > monPhysVar = socket_monPhysVar.getDataHandle();
  
  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {
    CFreal stateP = 0.0;
    CFreal diffStatesPPMinOne = 0.0;
    
    if (m_monitoredPhysVar < m_pData.size())
    {
      RealVector statePMinOne = *((*m_cellStates)[iSol]->getData()) - m_statesPMinOne[iSol];
      
      const bool RhoivtTv = getMethodData().getUpdateVarStr() == "RhoivtTv";
      
      if (RhoivtTv)
      {
        m_eulerVarSet2->computePhysicalData(*((*m_cellStates)[iSol]),m_pData);
        m_eulerVarSet2->computePhysicalData(m_statesPMinOne[iSol],m_pData2);
      }
      else 
      {
	m_eulerVarSet->computePhysicalData(*((*m_cellStates)[iSol]),m_pData);
        m_eulerVarSet->computePhysicalData(m_statesPMinOne[iSol],m_pData2);
      }
      
      stateP = m_pData[m_monitoredPhysVar];
      diffStatesPPMinOne = stateP - m_pData2[m_monitoredPhysVar];
    
      monPhysVar[(((*m_cellStates)[iSol]))->getLocalID()] = stateP;
    }
    else
    {
      stateP = (*((*m_cellStates)[iSol]))[m_monitoredVar];
      diffStatesPPMinOne = stateP - m_statesPMinOne[iSol][m_monitoredVar];
    }
    
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

CFreal LLAVFluxReconstructionNS::computePeclet()
{
  const CFreal machInf = m_eulerVarSet->getModel()->getMachInf();
  const CFreal velInf = m_eulerVarSet->getModel()->getVelInf();
  const CFreal pressInf = m_eulerVarSet->getModel()->getPressInf();
  const CFreal gamma = m_eulerVarSet->getModel()->getGamma();
  
//   cf_assert(machInf > 1.0);
//   cf_assert(velInf > 0.0);
//   cf_assert(pressInf > 0.0);
  
  const CFreal rhoInf = gamma*pressInf*machInf*machInf/(velInf*velInf);
  
  const CFreal factor = 90.0*4./3.*(m_order+1.)/(m_order+2.);
  
  //CFreal result = factor/(m_peclet*(machInf-1.0)*rhoInf);
  
  CFreal result = m_peclet;
  
  return result;
}

//////////////////////////////////////////////////////////////////////////////

void LLAVFluxReconstructionNS::setup()
{
  CFAUTOTRACE;
  
  // setup parent class
  LLAVFluxReconstruction::setup();
  
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
  
  

  
  // if needed, get the MS 
  
  m_gradsBackUp.resize(2);
  m_gradsBackUp[LEFT].resize(m_nbrSolPnts);
  m_gradsBackUp[RIGHT].resize(m_nbrSolPnts);

  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  { 
    for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
    {
      m_gradsBackUp[LEFT][iSol].push_back(new RealVector(m_dim));
      m_gradsBackUp[RIGHT][iSol].push_back(new RealVector(m_dim));
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void LLAVFluxReconstructionNS::unsetup()
{
  CFAUTOTRACE;
  
  for (CFuint iSol= 0; iSol < m_nbrSolPnts; ++iSol)
  {
    for (CFuint iGrad = 0; iGrad < m_nbrEqs; ++iGrad)
    {
      deletePtr(m_gradsBackUp[LEFT][iSol][iGrad]);  
      deletePtr(m_gradsBackUp[RIGHT][iSol][iGrad]);
    }
    m_gradsBackUp[LEFT][iSol].clear();
    m_gradsBackUp[RIGHT][iSol].clear();
  }
  m_gradsBackUp[LEFT].clear();
  m_gradsBackUp[RIGHT].clear();
  m_gradsBackUp.clear();
  
  // unsetup parent class
  LLAVFluxReconstruction::unsetup();
}


//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD

