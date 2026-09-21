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

#include "FluxReconstructionNavierStokes/ConvDiffJacobFluxReconstructionNS.hh"
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

MethodCommandProvider< ConvDiffJacobFluxReconstructionNS,
		       FluxReconstructionSolverData,
		       FluxReconstructionNavierStokesModule >
convDiffRHSJacobNSFluxReconstructionProvider("ConvDiffRHSJacobNS");
  
//////////////////////////////////////////////////////////////////////////////
  
ConvDiffJacobFluxReconstructionNS::ConvDiffJacobFluxReconstructionNS(const std::string& name) :
  ConvDiffJacobFluxReconstruction(name),
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
  m_tempGradTermJacob(),
  m_tempStatesJacob(),
  m_tempGradTermJacob2(),
  m_tempStatesJacob2(),
  m_unpertGradVars(),
  m_pertGradVars()
{
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffJacobFluxReconstructionNS::configure ( Config::ConfigArgs& args )
{
  ConvDiffJacobFluxReconstruction::configure(args);
} 
  
//////////////////////////////////////////////////////////////////////////////

void ConvDiffJacobFluxReconstructionNS::computeWaveSpeedUpdates(vector< CFreal >& waveSpeedUpd)
{
  // compute the wave speed updates for the neighbouring cells
  cf_assert(waveSpeedUpd.size() == 2);
  
  // here convective and artificial parts are added!
  ConvDiffJacobFluxReconstruction::computeWaveSpeedUpdates(waveSpeedUpd);
          
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

void ConvDiffJacobFluxReconstructionNS::computeInterfaceFlxCorrection()
{
  ConvDiffJacobFluxReconstruction::computeInterfaceFlxCorrection();
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffJacobFluxReconstructionNS::computeRiemannFluxJacobianNum(const CFreal resFactor)
{
  ConvDiffJacobFluxReconstruction::computeRiemannFluxJacobianNum(resFactor);
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffJacobFluxReconstructionNS::computeRiemannFluxToGradJacobianNum(const CFreal resFactor)
{
  ConvDiffJacobFluxReconstruction::computeRiemannFluxToGradJacobianNum(resFactor);
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffJacobFluxReconstructionNS::computeGradVarsToStateJacobianNum()
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
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffJacobFluxReconstructionNS::prepareFluxComputation()
{
  const bool isPerturb = this->getMethodData().isPerturb();
  const CFuint iPerturbVar = this->getMethodData().iPerturbVar();

  m_diffusiveVarSetNS->setComposition(m_avgSol, isPerturb, iPerturbVar);
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffJacobFluxReconstructionNS::computeGradients()
{
  ConvDiffJacobFluxReconstruction::computeGradients();
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffJacobFluxReconstructionNS::computeGradientFaceCorrections()
{
  ConvDiffJacobFluxReconstruction::computeGradientFaceCorrections();
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffJacobFluxReconstructionNS::setup()
{
  CFAUTOTRACE;
  
  // setup parent class
  ConvDiffJacobFluxReconstruction::setup();
  
  // get damping coeff
  m_dampCoeffDiff = getMethodData().getDiffDampCoefficient();
  
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
