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

#include "FluxReconstructionNavierStokes/DiffRHSFluxReconstructionNS.hh"
#include "FluxReconstructionNavierStokes/FluxReconstructionNavierStokes.hh"
#include "FluxReconstructionMethod/FluxReconstructionElementData.hh"
#include "NavierStokes/NavierStokesVarSet.hh"

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

MethodCommandProvider< DiffRHSFluxReconstructionNS,
		       FluxReconstructionSolverData,
		       FluxReconstructionNavierStokesModule >
diffRHSNSFluxReconstructionProvider("DiffRHSNS");
  
//////////////////////////////////////////////////////////////////////////////
  
DiffRHSFluxReconstructionNS::DiffRHSFluxReconstructionNS(const std::string& name) :
  DiffRHSFluxReconstruction(name),
  m_tempGradTermL(),
  m_tempGradTermR(),
  m_diffusiveVarSet(CFNULL),
  m_tempStatesL(),
  m_tempStatesR(),
  m_dampCoeff()
{
  //addConfigOptionsTo(this);
}
  
//////////////////////////////////////////////////////////////////////////////

void DiffRHSFluxReconstructionNS::computeWaveSpeedUpdates(vector< CFreal >& waveSpeedUpd)
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
    for (CFuint iFlx = 0; iFlx < m_cellFlx[iSide].size(); ++iFlx)
    {
      const CFreal jacobXJacobXIntCoef = m_faceJacobVecAbsSizeFlxPnts[iFlx]*
                                 m_faceJacobVecAbsSizeFlxPnts[iFlx]*
                                   (*m_faceIntegrationCoefs)[iFlx]*
                                   m_cflConvDiffRatio;
      const CFreal rho = m_diffusiveVarSet->getDensity(*(m_cellStatesFlxPnt[iSide][iFlx]));
      visc = dynVisc/rho/factorPr;
      
      // transform update states to physical data to calculate eigenvalues
      waveSpeedUpd[iSide] += visc*jacobXJacobXIntCoef/m_cellVolume[iSide];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void DiffRHSFluxReconstructionNS::prepareFluxComputation()
{
  const bool isPerturb = this->getMethodData().isPerturb();
  const CFuint iPerturbVar = this->getMethodData().iPerturbVar();

  m_diffusiveVarSet->setComposition(m_avgSol, isPerturb, iPerturbVar);
}

//////////////////////////////////////////////////////////////////////////////

void DiffRHSFluxReconstructionNS::computeInterfaceFlxCorrection()
{
    
  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    m_tempStatesL[iFlx] = m_cellStatesFlxPnt[LEFT][iFlx]->getData();
    m_tempStatesR[iFlx] = m_cellStatesFlxPnt[RIGHT][iFlx]->getData();
  }
  
  m_diffusiveVarSet->setGradientVars(m_tempStatesL,m_tempGradTermL,m_nbrFaceFlxPnts);
  m_diffusiveVarSet->setGradientVars(m_tempStatesR,m_tempGradTermR,m_nbrFaceFlxPnts);
    
  // Loop over the flux points to calculate FI
  for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
  { 
    // compute the average sol and grad to use the BR2 scheme
    for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
    {
      *(m_avgGrad[iVar]) = (*(m_cellGradFlxPnt[LEFT][iFlxPnt][iVar]) + *(m_cellGradFlxPnt[RIGHT][iFlxPnt][iVar]))/2.0;
       
      m_avgSol[iVar] = ((*(m_cellStatesFlxPnt[LEFT][iFlxPnt]))[iVar] + (*(m_cellStatesFlxPnt[RIGHT][iFlxPnt]))[iVar])/2.0; 
    }

    // damping factor
    const CFreal dampFactor = m_dampCoeff*m_faceInvCharLengths[iFlxPnt];

    // compute averaged (damped) gradients
    for (CFuint iGrad = 0; iGrad < m_nbrEqs; ++iGrad)
    {
      // compute damping term
      const RealVector dGradVarXNormal = (m_tempGradTermL(iGrad,iFlxPnt) - m_tempGradTermR(iGrad,iFlxPnt))*m_unitNormalFlxPnts[iFlxPnt];
      *m_avgGrad[iGrad] -= dampFactor*dGradVarXNormal;
    }
    
    prepareFluxComputation();
     
    // compute the Riemann flux
//     m_flxPntRiemannFlux[iFlxPnt] = m_diffusiveVarSet->getFlux(m_avgSol,m_avgGrad,m_unitNormalFlxPnts[iFlxPnt],0);
    computeFlux(m_avgSol,m_avgGrad,m_unitNormalFlxPnts[iFlxPnt],0,m_flxPntRiemannFlux[iFlxPnt]);
     
    // compute FI in the mapped coord frame
    m_cellFlx[LEFT][iFlxPnt] = (m_flxPntRiemannFlux[iFlxPnt])*m_faceJacobVecSizeFlxPnts[iFlxPnt][LEFT];
    m_cellFlx[RIGHT][iFlxPnt] = (m_flxPntRiemannFlux[iFlxPnt])*m_faceJacobVecSizeFlxPnts[iFlxPnt][RIGHT];
  }
}

//////////////////////////////////////////////////////////////////////////////

void DiffRHSFluxReconstructionNS::setup()
{
  CFAUTOTRACE;
  
  // setup parent class
  DiffRHSFluxReconstruction::setup();
  
  // get damping coeff
  m_dampCoeff = getMethodData().getDiffDampCoefficient();
  
  // get the diffusive varset
  m_diffusiveVarSet = (getMethodData().getDiffusiveVar()).d_castTo< NavierStokesVarSet >();
  
  m_tempGradTermL.resize(m_nbrEqs,m_nbrFaceFlxPnts);
  m_tempGradTermR.resize(m_nbrEqs,m_nbrFaceFlxPnts);
  
  m_tempStatesL.resize(m_nbrFaceFlxPnts);
  m_tempStatesR.resize(m_nbrFaceFlxPnts);
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD

