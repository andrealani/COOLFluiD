#include "Framework/MethodCommandProvider.hh"

#include "Framework/MeshData.hh"
#include "Framework/BaseTerm.hh"

#include "MathTools/MathFunctions.hh"

#include "FluxReconstructionNavierStokes/DiffBndCorrectionsRHSFluxReconstructionNS.hh"
#include "FluxReconstructionNavierStokes/FluxReconstructionNavierStokes.hh"
#include "FluxReconstructionMethod/FluxReconstructionElementData.hh"
#include "NavierStokes/NavierStokesVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider< DiffBndCorrectionsRHSFluxReconstructionNS, 
		       FluxReconstructionSolverData, 
		       FluxReconstructionNavierStokesModule >
DiffBndCorrectionsRHSNSFluxReconstructionProvider("DiffBndCorrectionsRHSNS");

//////////////////////////////////////////////////////////////////////////////

DiffBndCorrectionsRHSFluxReconstructionNS::DiffBndCorrectionsRHSFluxReconstructionNS(const std::string& name) :
  DiffBndCorrectionsRHSFluxReconstruction(name)
{
}

//////////////////////////////////////////////////////////////////////////////

DiffBndCorrectionsRHSFluxReconstructionNS::~DiffBndCorrectionsRHSFluxReconstructionNS()
{
}

//////////////////////////////////////////////////////////////////////////////

void DiffBndCorrectionsRHSFluxReconstructionNS::computeWaveSpeedUpdates(CFreal& waveSpeedUpd)
{
  CFreal visc = 1.0;
  SafePtr< NavierStokesVarSet > navierStokesVarSet = m_diffusiveVarSet.d_castTo< NavierStokesVarSet >();
  const CFreal dynVisc = navierStokesVarSet->getCurrDynViscosity();
  
  waveSpeedUpd = 0.0;
  for (CFuint iFlx = 0; iFlx < m_cellFlx.size(); ++iFlx)
  {
    const CFreal jacobXJacobXIntCoef = m_faceJacobVecAbsSizeFlxPnts[iFlx]*
                                 m_faceJacobVecAbsSizeFlxPnts[iFlx]*
                                   (*m_faceIntegrationCoefs)[iFlx]*
                                   m_cflConvDiffRatio;
    const CFreal rho = navierStokesVarSet->getDensity(*m_cellStatesFlxPnt[iFlx]);
    visc = dynVisc/rho;
				   
    // transform update states to physical data to calculate eigenvalues
    waveSpeedUpd += visc*jacobXJacobXIntCoef/m_cellVolume;
  }

}

//////////////////////////////////////////////////////////////////////////////

void DiffBndCorrectionsRHSFluxReconstructionNS::prepareFluxComputation()
{
  const bool isPerturb = this->getMethodData().isPerturb();
  const CFuint iPerturbVar = this->getMethodData().iPerturbVar();
  SafePtr< NavierStokesVarSet > navierStokesVarSet = m_diffusiveVarSet.d_castTo< NavierStokesVarSet >();
  navierStokesVarSet->setComposition(m_avgSol, isPerturb, iPerturbVar);
}

//////////////////////////////////////////////////////////////////////////////

//void DiffBndCorrectionsRHSFluxReconstructionNS::computeInterfaceFlxCorrection()
//{
//    
//  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
//  {
//    m_tempStatesL[iFlx] = m_cellStatesFlxPnt[LEFT][iFlx]->getData();
//    m_tempStatesR[iFlx] = m_cellStatesFlxPnt[RIGHT][iFlx]->getData();
//  }
//  
//  m_diffusiveVarSet->setGradientVars(m_tempStatesL,m_tempGradTermL,m_nbrFaceFlxPnts);
//  m_diffusiveVarSet->setGradientVars(m_tempStatesR,m_tempGradTermR,m_nbrFaceFlxPnts);
//    
//  // Loop over the flux points to calculate FI
//  for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
//  { 
//    // compute the average sol and grad to use the BR2 scheme
//    for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
//    {
//      *(m_avgGrad[iVar]) = (*(m_cellGradFlxPnt[LEFT][iFlxPnt][iVar]) + *(m_cellGradFlxPnt[RIGHT][iFlxPnt][iVar]))/2.0;
//       
//      m_avgSol[iVar] = ((*(m_cellStatesFlxPnt[LEFT][iFlxPnt]))[iVar] + (*(m_cellStatesFlxPnt[RIGHT][iFlxPnt]))[iVar])/2.0; 
//    }
//
//    // damping factor
//    const CFreal dampFactor = 1.0*m_faceInvCharLengths[iFlxPnt];
//
//    // compute averaged (damped) gradients
//    for (CFuint iGrad = 0; iGrad < m_nbrEqs; ++iGrad)
//    {
//      // compute damping term
//      const RealVector dGradVarXNormal = (m_tempGradTermL(iGrad,iFlxPnt) - m_tempGradTermR(iGrad,iFlxPnt))*m_unitNormalFlxPnts[iFlxPnt];
//      *m_avgGrad[iGrad] -= dampFactor*dGradVarXNormal;
//    }
//    
//    prepareFluxComputation();
//     
//    // compute the Riemann flux
////     m_flxPntRiemannFlux[iFlxPnt] = m_diffusiveVarSet->getFlux(m_avgSol,m_avgGrad,m_unitNormalFlxPnts[iFlxPnt],0);
//    computeFlux(m_avgSol,m_avgGrad,m_unitNormalFlxPnts[iFlxPnt],0,m_flxPntRiemannFlux[iFlxPnt]);
//     
//    // compute FI in the mapped coord frame
//    m_cellFlx[LEFT][iFlxPnt] = (m_flxPntRiemannFlux[iFlxPnt])*m_faceJacobVecSizeFlxPnts[iFlxPnt][LEFT];
//    m_cellFlx[RIGHT][iFlxPnt] = (m_flxPntRiemannFlux[iFlxPnt])*m_faceJacobVecSizeFlxPnts[iFlxPnt][RIGHT];
//  }
//}

//////////////////////////////////////////////////////////////////////////////

void DiffBndCorrectionsRHSFluxReconstructionNS::setup()
{
  CFAUTOTRACE;
  
  // setup parent class
  DiffBndCorrectionsRHSFluxReconstruction::setup();
  
  // get the diffusive varset
  m_diffusiveVarSet = (getMethodData().getDiffusiveVar()).d_castTo< NavierStokesVarSet >();
  
  m_tempGradTerm.resize(m_nbrEqs,m_nbrFaceFlxPnts);
  
  m_tempStates.resize(m_nbrFaceFlxPnts);
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluxReconstructionMethod

} // namespace COOLFluiD
