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

#include "FluxReconstructionMultiFluidMHD/LLAVJacobFluxReconstructionMFMHD.hh"
#include "FluxReconstructionMethod/FluxReconstruction.hh"
#include "FluxReconstructionMethod/FluxReconstructionElementData.hh"

#include "MultiFluidMHD/DiffMFMHD2DVarSet.hh"
#include "MultiFluidMHD/MultiFluidMHDVarSet.hh"
#include "MultiFluidMHD/EulerMFMHDTerm.hh"
#include "Maxwell/Maxwell2DProjectionVarSet.hh"

#include "FluxReconstructionMultiFluidMHD/FluxReconstructionMultiFluidMHD.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::NavierStokes;
using namespace COOLFluiD::Physics::MultiFluidMHD;
using namespace COOLFluiD::Physics::Maxwell;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace FluxReconstructionMethod {
    
//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider< LLAVJacobFluxReconstructionMFMHD,
		       FluxReconstructionSolverData,
		       FluxReconstructionModule >
LLAVJacobFluxReconstructionMFMHDFluxReconstructionProvider("LLAVJacobMFMHD");

//////////////////////////////////////////////////////////////////////////////
  
LLAVJacobFluxReconstructionMFMHD::LLAVJacobFluxReconstructionMFMHD(const std::string& name) :
  LLAVJacobFluxReconstructionNS(name),
  m_varSet(CFNULL)
  {
  }

//////////////////////////////////////////////////////////////////////////////

void LLAVJacobFluxReconstructionMFMHD::configure ( Config::ConfigArgs& args )
{
  LLAVJacobFluxReconstructionNS::configure(args);
}  

//////////////////////////////////////////////////////////////////////////////

void LLAVJacobFluxReconstructionMFMHD::computeSmoothness()
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

void LLAVJacobFluxReconstructionMFMHD::computeSmoothness(const CFuint side)
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

void LLAVJacobFluxReconstructionMFMHD::setup()
{
  CFAUTOTRACE;
  
  // setup parent class
  LLAVJacobFluxReconstructionNS::setup();

  // get Euler 2D varset
  m_varSet = getMethodData().getUpdateVar().d_castTo< MultiFluidMHDVarSet<Maxwell2DProjectionVarSet> >();

  if (m_varSet.isNull())
  {
    throw Common::ShouldNotBeHereException (FromHere(),"Update variable set is not MultiFluid2DMHDVarSet in LLAVJacobFRMFMHD!\n");
  }

  // resize the physical data for internal and ghost solution points
  m_varSet->getModel()->resizePhysicalData(m_pData);
  m_varSet->getModel()->resizePhysicalData(m_pData2);
  
}

//////////////////////////////////////////////////////////////////////////////

void LLAVJacobFluxReconstructionMFMHD::unsetup()
{
  CFAUTOTRACE;
  
  // unsetup parent class
  LLAVJacobFluxReconstructionNS::unsetup();
}


//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD

