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

#include "NavierStokes/EulerTerm.hh"

#include "NavierStokes/NavierStokesVarSet.hh"
#include "NavierStokes/EulerVarSet.hh"

#include "FluxReconstructionNavierStokes/LLAVDiffFluxReconstructionNS.hh"
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

MethodCommandProvider< LLAVDiffFluxReconstructionNS,
		       FluxReconstructionSolverData,
		       FluxReconstructionModule >
LLAVDiffFluxReconstructionNSFluxReconstructionProvider("LLAVDiffNS");

//////////////////////////////////////////////////////////////////////////////
  
LLAVDiffFluxReconstructionNS::LLAVDiffFluxReconstructionNS(const std::string& name) :
  LLAVDiffFluxReconstruction(name),
  m_gradsBackUp(),
  m_eulerVarSet(CFNULL),
  m_msEulerTerm(CFNULL),
  m_nbrSpecies(),
  m_pData(),
  m_pData2(),
  m_eulerVarSet2(CFNULL),
  m_tempGradTermL(),
  m_tempGradTermR(),
  m_diffusiveVarSet(CFNULL),
  m_tempStatesL(),
  m_tempStatesR(),
  m_cellGradsAV(),
  m_avgGradAV(),
  m_cellGradFlxPntAV(),
  m_dampCoeff()
  {
  }

//////////////////////////////////////////////////////////////////////////////

void LLAVDiffFluxReconstructionNS::configure ( Config::ConfigArgs& args )
{
  LLAVDiffFluxReconstruction::configure(args);
}  

//////////////////////////////////////////////////////////////////////////////

void LLAVDiffFluxReconstructionNS::prepareFluxComputation()
{
  const bool isPerturb = this->getMethodData().isPerturb();
  const CFuint iPerturbVar = this->getMethodData().iPerturbVar();

  m_diffusiveVarSet->setComposition(m_avgSol, isPerturb, iPerturbVar);
}

//////////////////////////////////////////////////////////////////////////////

void LLAVDiffFluxReconstructionNS::computeWaveSpeedUpdates(vector< CFreal >& waveSpeedUpd)
{
  // compute the wave speed updates for the neighbouring cells
  cf_assert(waveSpeedUpd.size() == 2);
  CFreal visc = 1.0;

  const CFreal dynVisc = m_diffusiveVarSet->getCurrDynViscosity();
  const CFreal factorPr = min(m_diffusiveVarSet->getModel().getPrandtl(),1.0); 
  cf_assert(factorPr>0.0);
  
  for (CFuint iSide = 0; iSide < 2; ++iSide)
  {
    waveSpeedUpd[iSide] = 0.0;
    //for (CFuint iFlx = 0; iFlx < m_cellFlx[iSide].size(); ++iFlx)
    for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
    {
      const CFreal jacobXJacobXIntCoef = m_faceJacobVecAbsSizeFlxPnts[iFlx]*
                                 m_faceJacobVecAbsSizeFlxPnts[iFlx]*
                                   (*m_faceIntegrationCoefs)[iFlx]*
                                   m_cflConvDiffRatio;
      const CFreal rho = m_diffusiveVarSet->getDensity(*(m_cellStatesFlxPnt[iSide][iFlx]));
      visc = dynVisc/rho/factorPr;
      
      if (m_addUpdCoeff)
      {
        const CFreal epsilon = 0.5*(m_epsilonLR[LEFT][iFlx]+m_epsilonLR[RIGHT][iFlx]);
        const CFreal viscCoef = computeViscCoef(m_cellStatesFlxPnt[iSide][iFlx]);
        visc += epsilon*viscCoef;
      }
      
      // transform update states to physical data to calculate eigenvalues
      waveSpeedUpd[iSide] += visc*jacobXJacobXIntCoef/m_cellVolume[iSide];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void LLAVDiffFluxReconstructionNS::setFaceData(CFuint faceID)
{
  LLAVDiffFluxReconstruction::setFaceData(faceID);

  // get the gradients datahandle
  DataHandle< vector< RealVector > > gradientsAV = socket_gradientsAV.getDataHandle();

  for (CFuint iSide = 0; iSide < 2; ++iSide)
  {
    for (CFuint iState = 0; iState < m_nbrSolPnts; ++iState)
    {
      const CFuint stateID = (*(m_states[iSide]))[iState]->getLocalID();
      m_cellGradsAV[iSide][iState] = &gradientsAV[stateID];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void LLAVDiffFluxReconstructionNS::setCellData()
{
  LLAVDiffFluxReconstruction::setCellData();

  // get the gradients datahandle
  DataHandle< vector< RealVector > > gradientsAV = socket_gradientsAV.getDataHandle();

  // get the grads in the current cell
  for (CFuint iState = 0; iState < m_nbrSolPnts; ++iState)
  {
    const CFuint stateID = (*(m_cellStates))[iState]->getLocalID();
    m_cellGradsAV[0][iState] = &gradientsAV[stateID];
  }
}

//////////////////////////////////////////////////////////////////////////////

void LLAVDiffFluxReconstructionNS::computeInterfaceFlxCorrection()
{
  LLAVDiffFluxReconstruction::computeInterfaceFlxCorrection();
}

//////////////////////////////////////////////////////////////////////////////

CFreal LLAVDiffFluxReconstructionNS::computeViscCoef(RealVector* state)
{
  CFreal result = 1.0;
  
  return result;
}

//////////////////////////////////////////////////////////////////////////////

void LLAVDiffFluxReconstructionNS::computeSmoothness()
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
  const RealVector& coords = *((*m_cellNodes)[0]->getData());
  const CFreal d = sqrt(coords[0]*coords[0]+coords[1]*coords[1]);

  if (sNum <= MathTools::MathConsts::CFrealEps() || sDenom <= MathTools::MathConsts::CFrealEps() ) //  || d < 1.2 || coords[0] < 0.1
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

CFreal LLAVDiffFluxReconstructionNS::computePeclet()
{ 
  CFreal result = m_peclet;
  
  return result;
}

//////////////////////////////////////////////////////////////////////////////

void LLAVDiffFluxReconstructionNS::computeFlxPntStatesAndGrads()
{
  LLAVDiffFluxReconstruction::computeFlxPntStatesAndGrads();
}

//////////////////////////////////////////////////////////////////////////////

void LLAVDiffFluxReconstructionNS::computeDivDiscontFlx(vector< RealVector >& residuals)
{
  LLAVDiffFluxReconstruction::computeDivDiscontFlx(residuals);
}

//////////////////////////////////////////////////////////////////////////////

void LLAVDiffFluxReconstructionNS::setup()
{
  CFAUTOTRACE;
  
  // setup parent class
  LLAVDiffFluxReconstruction::setup();
  
  const bool RhoivtTv = getMethodData().getUpdateVarStr() == "RhoivtTv";
  
  // get damping coeff
  m_dampCoeff = getMethodData().getDiffDampCoefficient();
  
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
  
  // get the diffusive varset
  m_diffusiveVarSet = (getMethodData().getDiffusiveVar()).d_castTo< NavierStokesVarSet >();
  
  m_tempGradTermL.resize(m_nbrEqs,m_nbrFaceFlxPnts);
  m_tempGradTermR.resize(m_nbrEqs,m_nbrFaceFlxPnts);
  
  m_tempStatesL.resize(m_nbrFaceFlxPnts);
  m_tempStatesR.resize(m_nbrFaceFlxPnts);
  
  m_cellGradsAV.resize(2);
  m_cellGradsAV[LEFT].resize(m_nbrSolPnts);
  m_cellGradsAV[RIGHT].resize(m_nbrSolPnts);
  
  m_cellGradFlxPntAV.resize(2);
  m_avgGradAV.reserve(m_nbrEqs);
  m_cellGradFlxPntAV[LEFT].resize(m_nbrFaceFlxPnts);
  m_cellGradFlxPntAV[RIGHT].resize(m_nbrFaceFlxPnts);
  
  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  { 
    for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
    {
      m_cellGradFlxPntAV[LEFT][iFlx].push_back(new RealVector(m_dim));
      m_cellGradFlxPntAV[RIGHT][iFlx].push_back(new RealVector(m_dim));
    }
  }
  
  for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
  {
    m_avgGradAV[iVar] = new RealVector(m_dim);
  }
}

//////////////////////////////////////////////////////////////////////////////

void LLAVDiffFluxReconstructionNS::unsetup()
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
  for (CFuint iFlx= 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    for (CFuint iGrad = 0; iGrad < m_nbrEqs; ++iGrad)
    {
      deletePtr(m_cellGradFlxPntAV[LEFT][iFlx][iGrad]);  
      deletePtr(m_cellGradFlxPntAV[RIGHT][iFlx][iGrad]);
    }
    m_cellGradFlxPntAV[LEFT][iFlx].clear();
    m_cellGradFlxPntAV[RIGHT][iFlx].clear();
  }
  
  m_gradsBackUp[LEFT].clear();
  m_gradsBackUp[RIGHT].clear();
  m_gradsBackUp.clear();
  
  for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
  {
    deletePtr(m_avgGradAV[iVar]); 
  }
  m_avgGradAV.clear();
  m_cellGradFlxPntAV[LEFT].clear();
  m_cellGradFlxPntAV[RIGHT].clear();
  m_cellGradFlxPntAV.clear();
  
  // unsetup parent class
  LLAVDiffFluxReconstruction::unsetup();
}


//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD

