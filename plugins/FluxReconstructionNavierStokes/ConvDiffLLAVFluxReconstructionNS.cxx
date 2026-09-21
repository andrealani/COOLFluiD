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

#include "FluxReconstructionNavierStokes/ConvDiffLLAVFluxReconstructionNS.hh"
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

MethodCommandProvider< ConvDiffLLAVFluxReconstructionNS,
		       FluxReconstructionSolverData,
		       FluxReconstructionNavierStokesModule >
convDiffLLAVRHSNSFluxReconstructionProvider("ConvDiffLLAVRHSNS");
  
//////////////////////////////////////////////////////////////////////////////
  
ConvDiffLLAVFluxReconstructionNS::ConvDiffLLAVFluxReconstructionNS(const std::string& name) :
  ConvDiffLLAVFluxReconstruction(name),
  m_tempGradTerm(),
  m_tempGradTermL(),
  m_tempGradTermR(),
  m_diffusiveVarSetNS(CFNULL),
  m_tempStatesL(),
  m_tempStatesR(),
  m_tempStatesL2(),
  m_tempStatesR2(),
  m_tempStatesCell(),
  m_dampCoeffDiff(),
  m_eulerVarSet(CFNULL),
  m_msEulerTerm(CFNULL),
  m_nbrSpecies(),
  m_pData(),
  m_pData2(),
  m_eulerVarSet2(CFNULL),
  m_tempSolVarState(),
  m_tempSolVarState2()
{
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffLLAVFluxReconstructionNS::configure ( Config::ConfigArgs& args )
{
  ConvDiffLLAVFluxReconstruction::configure(args);
} 
  
//////////////////////////////////////////////////////////////////////////////

void ConvDiffLLAVFluxReconstructionNS::computeWaveSpeedUpdates(vector< CFreal >& waveSpeedUpd)
{
  // compute the wave speed updates for the neighbouring cells
  cf_assert(waveSpeedUpd.size() == 2);
  
  // here convective and artificial parts are added!
  ConvDiffLLAVFluxReconstruction::computeWaveSpeedUpdates(waveSpeedUpd);
          
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

void ConvDiffLLAVFluxReconstructionNS::computeInterfaceFlxCorrection()
{
  ConvDiffLLAVFluxReconstruction::computeInterfaceFlxCorrection();
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffLLAVFluxReconstructionNS::computeCellGradTerm(RealMatrix& gradTerm)
{   
  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {
    m_tempStatesCell[iSol] = (*m_cellStates)[iSol]->getData();
  }
  
  m_diffusiveVarSetNS->setGradientVars(m_tempStatesCell,gradTerm,m_nbrSolPnts);
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffLLAVFluxReconstructionNS::computeFaceGradTerms(RealMatrix& gradTermL, RealMatrix& gradTermR)
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

void ConvDiffLLAVFluxReconstructionNS::prepareFluxComputation()
{
  const bool isPerturb = this->getMethodData().isPerturb();
  const CFuint iPerturbVar = this->getMethodData().iPerturbVar();

  m_diffusiveVarSetNS->setComposition(m_avgSol, isPerturb, iPerturbVar);
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffLLAVFluxReconstructionNS::computeGradients()
{
  ConvDiffLLAVFluxReconstruction::computeGradients();
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffLLAVFluxReconstructionNS::computeGradientFaceCorrections()
{
  ConvDiffLLAVFluxReconstruction::computeGradientFaceCorrections();
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffLLAVFluxReconstructionNS::computeSmoothness()
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

void ConvDiffLLAVFluxReconstructionNS::computeSmoothness(const CFuint side)
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

      m_eulerVarSet->computePhysicalData(*((*m_states[side])[iSol]),m_pData);
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
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffLLAVFluxReconstructionNS::computeEpsilon0()
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
  
  m_epsilon0 = max(h*wavespeed*(2.0/peclet - m_subcellRes/peclet),0.0);
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffLLAVFluxReconstructionNS::computeEpsilon0(const CFuint side)
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

void ConvDiffLLAVFluxReconstructionNS::setup()
{
  CFAUTOTRACE;
  
  // setup parent class
  ConvDiffLLAVFluxReconstruction::setup();
  
  // get damping coeff
  m_dampCoeffDiff = getMethodData().getDiffDampCoefficient();
  
  // get the diffusive varset
  m_diffusiveVarSetNS = (getMethodData().getDiffusiveVar()).d_castTo< NavierStokesVarSet >();
  
  m_tempSolVarState.resize(m_nbrEqs);
  m_tempSolVarState2.resize(m_nbrEqs);
  
  m_tempGradTerm.resize(m_nbrEqs,m_nbrSolPnts);
  
  m_tempGradTermL.resize(m_nbrEqs,m_nbrFaceFlxPnts);
  m_tempGradTermR.resize(m_nbrEqs,m_nbrFaceFlxPnts);
  
  m_tempStatesL.resize(m_nbrFaceFlxPnts);
  m_tempStatesR.resize(m_nbrFaceFlxPnts);
  m_tempStatesL2.resize(m_nbrFaceFlxPnts);
  m_tempStatesR2.resize(m_nbrFaceFlxPnts);
  
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

