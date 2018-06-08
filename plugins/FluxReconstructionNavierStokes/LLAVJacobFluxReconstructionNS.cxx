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

#include "FluxReconstructionNavierStokes/LLAVJacobFluxReconstructionNS.hh"
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

MethodCommandProvider< LLAVJacobFluxReconstructionNS,
		       FluxReconstructionSolverData,
		       FluxReconstructionModule >
LLAVJacobFluxReconstructionNSFluxReconstructionProvider("LLAVJacobNS");

//////////////////////////////////////////////////////////////////////////////
  
LLAVJacobFluxReconstructionNS::LLAVJacobFluxReconstructionNS(const std::string& name) :
  LLAVJacobFluxReconstruction(name),
  m_eulerVarSet(CFNULL),
  m_msEulerTerm(CFNULL),
  m_nbrSpecies(),
  m_pData(),
  m_pData2(),
  m_eulerVarSet2(CFNULL)
  {
  }

//////////////////////////////////////////////////////////////////////////////

void LLAVJacobFluxReconstructionNS::configure ( Config::ConfigArgs& args )
{
  LLAVJacobFluxReconstruction::configure(args);
}  

//////////////////////////////////////////////////////////////////////////////

void LLAVJacobFluxReconstructionNS::setFaceData(CFuint faceID)
{
  LLAVJacobFluxReconstruction::setFaceData(faceID);
  
  if (getMethodData().getUpdateVarStr() == "Cons" && getMethodData().hasDiffTerm())
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

void LLAVJacobFluxReconstructionNS::setCellData()
{
  LLAVJacobFluxReconstruction::setCellData();
  
  if (getMethodData().getUpdateVarStr() == "Cons" && getMethodData().hasDiffTerm())
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

void LLAVJacobFluxReconstructionNS::setFaceNeighbourGradients()
{
  if (getMethodData().getUpdateVarStr() == "Cons" && getMethodData().hasDiffTerm())
  {
    // get the gradients datahandle
    DataHandle< vector< RealVector > > gradientsAV = socket_gradientsAV.getDataHandle();

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
            m_faceNghbrGrads[iSide][iFace][iSide2][iState] = &gradientsAV[stateID];
          }
        }
      }
    }
  }
  else
  {
    LLAVJacobFluxReconstruction::setFaceNeighbourGradients();
  }
}

//////////////////////////////////////////////////////////////////////////////

CFreal LLAVJacobFluxReconstructionNS::computePeclet()
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

void LLAVJacobFluxReconstructionNS::computeSmoothness()
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

void LLAVJacobFluxReconstructionNS::setup()
{
  CFAUTOTRACE;
  
  // setup parent class
  LLAVJacobFluxReconstruction::setup();
  
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

void LLAVJacobFluxReconstructionNS::unsetup()
{
  CFAUTOTRACE;
  
  // unsetup parent class
  LLAVJacobFluxReconstruction::unsetup();
}


//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD

