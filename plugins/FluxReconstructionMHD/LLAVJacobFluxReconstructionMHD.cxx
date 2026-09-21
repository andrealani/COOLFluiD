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

#include "FluxReconstructionMHD/LLAVJacobFluxReconstructionMHD.hh"
#include "FluxReconstructionMethod/FluxReconstruction.hh"
#include "FluxReconstructionMethod/FluxReconstructionElementData.hh"

#include "MHD/MHD3DProjectionVarSet.hh"

#include "FluxReconstructionMHD/FluxReconstructionMHD.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::MHD;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace FluxReconstructionMethod {
    
//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider< LLAVJacobFluxReconstructionMHD,
		       FluxReconstructionSolverData,
		       FluxReconstructionModule >
LLAVJacobFluxReconstructionMHDFluxReconstructionProvider("LLAVJacobMHD");

//////////////////////////////////////////////////////////////////////////////
  
LLAVJacobFluxReconstructionMHD::LLAVJacobFluxReconstructionMHD(const std::string& name) :
  LLAVJacobFluxReconstruction(name),
  m_varSet(CFNULL),
  m_pData(),
  m_pData2()
  {
  }

//////////////////////////////////////////////////////////////////////////////

void LLAVJacobFluxReconstructionMHD::configure ( Config::ConfigArgs& args )
{
  LLAVJacobFluxReconstruction::configure(args);
}  

//////////////////////////////////////////////////////////////////////////////

void LLAVJacobFluxReconstructionMHD::setFaceData(CFuint faceID)
{
  LLAVJacobFluxReconstruction::setFaceData(faceID);
  
  if (getMethodData().getUpdateVarStr() != "Puvt" && getMethodData().hasDiffTerm())
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

CFreal LLAVJacobFluxReconstructionMHD::computePeclet()
{
//  const CFreal machInf = m_eulerVarSet->getModel()->getMachInf();
//  const CFreal velInf = m_eulerVarSet->getModel()->getVelInf();
//  const CFreal pressInf = m_eulerVarSet->getModel()->getPressInf();
//  const CFreal gamma = m_eulerVarSet->getModel()->getGamma();
//  
////   cf_assert(machInf > 1.0);
////   cf_assert(velInf > 0.0);
////   cf_assert(pressInf > 0.0);
//  
//  const CFreal rhoInf = gamma*pressInf*machInf*machInf/(velInf*velInf);
//  
//  const CFreal factor = 90.0*4./3.*(m_order+1.)/(m_order+2.);
//  
//  //CFreal result = factor/(m_peclet*(machInf-1.0)*rhoInf);
  CFreal result = m_peclet;
  
  return result;
}

//////////////////////////////////////////////////////////////////////////////

void LLAVJacobFluxReconstructionMHD::computeSmoothness()
{ 
  /*CFreal sNum = 0.0;
  
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
            
      m_varSet->computePhysicalData(*((*m_cellStates)[iSol]),m_pData);
      m_varSet->computePhysicalData(m_statesPMinOne[iSol],m_pData2);
      
      stateP = m_pData[m_monitoredPhysVar];
      diffStatesPPMinOne = stateP - m_pData2[m_monitoredPhysVar];
    
      monPhysVar[(((*m_cellStates)[iSol]))->getLocalID()] = stateP;
    }
    else
    {
      stateP = (*((*m_cellStates)[iSol]))[m_monitoredVar];
      diffStatesPPMinOne = stateP - m_statesPMinOne[iSol][m_monitoredVar];
    }
    
    sNum += (diffStatesPPMinOne+1.e-7)*(diffStatesPPMinOne+1.e-7);
    sDenom += (stateP+1.e-7)*(stateP+1.e-7);
  }
  if (sNum <= MathTools::MathConsts::CFrealEps() || sDenom <= MathTools::MathConsts::CFrealEps())
  {
    m_s = -100.0;
  }
  else
  {
    m_s = log10(sNum/sDenom);
  }*/

  std::vector<CFreal> modalCoeffs(m_nbrSolPnts, 0.0); // modal coefficients for a single variable
  CFreal max_indicator = 0.0; // will store the max(E_var) across all variables
  CFuint m_nbrSolPntsMinTwo = (m_order > 1) ? ((m_order-1)*(m_order-1)*(m_order)/2) : 1 ;
  CFuint m_nbrSolPntsMinOne = (m_order)*(m_order)*(m_order+1)/2;
  //for (CFuint iEq = 0; iEq < 1; ++iEq) // loop over variables (equations)
  //{
    // Step 1: Load nodal values for this variable
    for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
    {
      CFreal u = (*((*m_cellStates)[iSol]))[1];
      CFreal v = (*((*m_cellStates)[iSol]))[2];
      CFreal w = (*((*m_cellStates)[iSol]))[3];
      m_tempSolPntVec[iSol] = sqrt(u*u+v*v+w*w);
      /*CFreal rho = (*((*m_cellStates)[iSol]))[0];
      CFreal rhoU = (*((*m_cellStates)[iSol]))[1];
      CFreal rhoV = (*((*m_cellStates)[iSol]))[2];
      //CFreal rhoW = (*((*m_cellStates)[iSol]))[3];
      CFreal rhoE = (*((*m_cellStates)[iSol]))[3];
      CFreal p = 0.4*(rhoE-0.5*(rhoU*rhoU+rhoV*rhoV)/rho); //(*((*m_cellStates)[iSol]))[7]; //
      m_tempSolPntVec[iSol] = rho*p;*/
    }

    // Step 2: Transform to modal space: modalCoeffs = V^{-1} * nodalValues
    m_tempSolPntVec2 = m_vdmInv * m_tempSolPntVec;

    for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
    {
      modalCoeffs[iSol] = m_tempSolPntVec2[iSol];
    }

    // Step 3: Compute modal energy and high-mode content
    CFreal energy_total = 0.0;
    CFreal energy_pMin1 = 0.0;
    CFreal mN2 = 0.0;
    CFreal mNminOne2 = 0.0;

    for (CFuint j = 0; j < m_nbrSolPnts; ++j)
    {
      CFreal mj2 = modalCoeffs[j] * modalCoeffs[j];
      energy_total += mj2;
      if (j > m_nbrSolPntsMinOne - 1) mN2 += mj2;
    }

    for (CFuint j = m_nbrSolPntsMinTwo - 1; j < m_nbrSolPntsMinOne; ++j)
    {
      CFreal mj2 = modalCoeffs[j] * modalCoeffs[j];
      mNminOne2 += mj2;
    }

    for (CFuint j = 0; j < m_nbrSolPntsMinOne; ++j)
    {
      CFreal mj2 = modalCoeffs[j] * modalCoeffs[j];
      energy_pMin1 += mj2;
    }

    CFreal E1 = mN2 / std::max(energy_total, MathTools::MathConsts::CFrealEps());
    CFreal E2 = mNminOne2 / std::max(energy_pMin1, MathTools::MathConsts::CFrealEps());
    CFreal E_var = std::max(E1, E2); // Energy indicator for this variable
    E1 =  mN2 / mNminOne2;
    if (true) {E_var = E1 ;}

    max_indicator = std::max(max_indicator, E_var); // take max across variables
  //}

  m_s = max_indicator; // final smoothness value for this cell

  if (m_s>0)
  {
    m_s = log10(m_s);
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

void LLAVJacobFluxReconstructionMHD::computeSmoothness(const CFuint side)
{ 
  CFreal sNum = 0.0;
  
  CFreal sDenom = 0.0;
  
  // get datahandle
  DataHandle< CFreal > monPhysVar = socket_monPhysVar.getDataHandle();
  
  /*for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {
    CFreal stateP = 0.0;
    CFreal diffStatesPPMinOne = 0.0;
    
    if (m_monitoredPhysVar < m_pData.size())
    {
      RealVector statePMinOne = *((*m_states[side])[iSol]->getData()) - m_statesPMinOne[iSol];
            
      m_varSet->computePhysicalData(*((*m_states[side])[iSol]),m_pData);
      m_varSet->computePhysicalData(m_statesPMinOne[iSol],m_pData2);
      
      stateP = m_pData[m_monitoredPhysVar];
      diffStatesPPMinOne = stateP - m_pData2[m_monitoredPhysVar];
    
      monPhysVar[(((*m_states[side])[iSol]))->getLocalID()] = stateP;
    }
    else
    {
      stateP = (*((*m_states[side])[iSol]))[m_monitoredVar];
      diffStatesPPMinOne = stateP - m_statesPMinOne[iSol][m_monitoredVar];
    }
    
    sNum += (diffStatesPPMinOne+1.e-7)*(diffStatesPPMinOne+1.e-7);
    sDenom += (stateP+1.e-7)*(stateP+1.e-7);
  }
  if (sNum <= MathTools::MathConsts::CFrealEps() || sDenom <= MathTools::MathConsts::CFrealEps())
  {
    m_s = -100.0;
  }
  else
  {
    m_s = log10(sNum/sDenom);
  }*/

  std::vector<CFreal> modalCoeffs(m_nbrSolPnts, 0.0); // modal coefficients for a single variable
  CFreal max_indicator = 0.0; // will store the max(E_var) across all variables
  CFuint m_nbrSolPntsMinTwo = (m_order > 1) ? ((m_order-1)*(m_order-1)*(m_order)/2) : 1 ;
  CFuint m_nbrSolPntsMinOne = (m_order)*(m_order)*(m_order+1)/2;
  //for (CFuint iEq = 0; iEq < 1; ++iEq) // loop over variables (equations)
  //{
    // Step 1: Load nodal values for this variable
    for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
    {
      CFreal u = (*((*m_states[side])[iSol]))[1];
      CFreal v = (*((*m_states[side])[iSol]))[2];
      CFreal w = (*((*m_states[side])[iSol]))[3];
      m_tempSolPntVec[iSol] = sqrt(u*u+v*v+w*w);
      /*CFreal rho = (*((*m_cellStates)[iSol]))[0];
      CFreal rhoU = (*((*m_cellStates)[iSol]))[1];
      CFreal rhoV = (*((*m_cellStates)[iSol]))[2];
      //CFreal rhoW = (*((*m_cellStates)[iSol]))[3];
      CFreal rhoE = (*((*m_cellStates)[iSol]))[3];
      CFreal p = 0.4*(rhoE-0.5*(rhoU*rhoU+rhoV*rhoV)/rho); //(*((*m_cellStates)[iSol]))[7]; //
      m_tempSolPntVec[iSol] = rho*p;*/
    }

    // Step 2: Transform to modal space: modalCoeffs = V^{-1} * nodalValues
    m_tempSolPntVec2 = m_vdmInv * m_tempSolPntVec;

    for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
    {
      modalCoeffs[iSol] = m_tempSolPntVec2[iSol];
    }

    // Step 3: Compute modal energy and high-mode content
    CFreal energy_total = 0.0;
    CFreal energy_pMin1 = 0.0;
    CFreal mN2 = 0.0;
    CFreal mNminOne2 = 0.0;

    for (CFuint j = 0; j < m_nbrSolPnts; ++j)
    {
      CFreal mj2 = modalCoeffs[j] * modalCoeffs[j];
      energy_total += mj2;
      if (j > m_nbrSolPntsMinOne - 1) mN2 += mj2;
    }

    for (CFuint j = m_nbrSolPntsMinTwo - 1; j < m_nbrSolPntsMinOne; ++j)
    {
      CFreal mj2 = modalCoeffs[j] * modalCoeffs[j];
      mNminOne2 += mj2;
    }

    for (CFuint j = 0; j < m_nbrSolPntsMinOne; ++j)
    {
      CFreal mj2 = modalCoeffs[j] * modalCoeffs[j];
      energy_pMin1 += mj2;
    }

    CFreal E1 = mN2 / std::max(energy_total, MathTools::MathConsts::CFrealEps());
    CFreal E2 = mNminOne2 / std::max(energy_pMin1, MathTools::MathConsts::CFrealEps());
    CFreal E_var = std::max(E1, E2); // Energy indicator for this variable
    E1 =  mN2 / mNminOne2;
    if (true) {E_var = E1 ;}

    max_indicator = std::max(max_indicator, E_var); // take max across variables
  //}

  m_s = max_indicator; // final smoothness value for this cell

  if (m_s>0)
  {
    m_s = log10(m_s);
  }
}

//////////////////////////////////////////////////////////////////////////////

void LLAVJacobFluxReconstructionMHD::setup()
{
  CFAUTOTRACE;
  
  // setup parent class
  LLAVJacobFluxReconstruction::setup();
  
  // get the update varset
  m_updateVarSet = getMethodData().getUpdateVar();

  vector< FluxReconstructionElementData* >& frLocalData = getMethodData().getFRLocalData();

  m_vdmInv = *(frLocalData[0]->getVandermondeMatrixInv());
  
  // get 3D varset
  m_varSet = getMethodData().getUpdateVar().d_castTo< MHD3DProjectionVarSet >();

  if (m_varSet.isNull())
  {
    throw Common::ShouldNotBeHereException (FromHere(),"Update variable set is not MHD3DProjectionVarSet in LLAVJacobFRMHD!\n");
  }

  // resize the physical data for internal and ghost solution points
  m_varSet->getModel()->resizePhysicalData(m_pData);
  m_varSet->getModel()->resizePhysicalData(m_pData2);

}

//////////////////////////////////////////////////////////////////////////////

void LLAVJacobFluxReconstructionMHD::unsetup()
{
  CFAUTOTRACE;
  
  // unsetup parent class
  LLAVJacobFluxReconstruction::unsetup();
}


//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD

