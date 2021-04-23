#include "Framework/MethodCommandProvider.hh"

#include "Framework/MeshData.hh"
#include "Framework/BaseTerm.hh"

#include "MathTools/MathFunctions.hh"

#include "FluxReconstructionTurb/DiffBndCorrectionsRHSJacobFluxReconstructionGammaAlpha.hh"
#include "FluxReconstructionNavierStokes/FluxReconstructionNavierStokes.hh"
#include "FluxReconstructionMethod/FluxReconstructionElementData.hh"
#include "NavierStokes/NavierStokesVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::NavierStokes;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider< DiffBndCorrectionsRHSJacobFluxReconstructionGammaAlpha, 
		       FluxReconstructionSolverData, 
		       FluxReconstructionNavierStokesModule >
DiffBndCorrectionsRHSJacobGammaAlphaFluxReconstructionProvider("DiffBndCorrectionsRHSJacobGammaAlpha");

//////////////////////////////////////////////////////////////////////////////

DiffBndCorrectionsRHSJacobFluxReconstructionGammaAlpha::DiffBndCorrectionsRHSJacobFluxReconstructionGammaAlpha(const std::string& name) :
  DiffBndCorrectionsRHSJacobFluxReconstructionTurb(name)
{
}

//////////////////////////////////////////////////////////////////////////////

DiffBndCorrectionsRHSJacobFluxReconstructionGammaAlpha::~DiffBndCorrectionsRHSJacobFluxReconstructionGammaAlpha()
{
}

//////////////////////////////////////////////////////////////////////////////

void DiffBndCorrectionsRHSJacobFluxReconstructionGammaAlpha::setup()
{
  DiffBndCorrectionsRHSJacobFluxReconstructionTurb::setup();
}

//////////////////////////////////////////////////////////////////////////////

void DiffBndCorrectionsRHSJacobFluxReconstructionGammaAlpha::unsetup()
{
  DiffBndCorrectionsRHSJacobFluxReconstructionTurb::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

void DiffBndCorrectionsRHSJacobFluxReconstructionGammaAlpha::computeFlxPntStates()
{
  // Loop over flux points to extrapolate the states and gradients to the flux points
  for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
  {
    // reset the states in flx pnts
    *(m_cellStatesFlxPnt[iFlxPnt]) = 0.0;
    
    // reset the grads in flx pnts
    for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
    {
      *(m_cellGradFlxPnt[iFlxPnt][iVar]) = 0.0;
    }
  
    // index of current flx pnt
    const CFuint currFlxIdx = (*m_faceFlxPntConn)[m_orient][iFlxPnt];
    
    // Loop over sol points to add the contributions to each sol pnt
    for (CFuint iSol = 0; iSol < m_nbrSolDep; ++iSol)
    {
      const CFuint solIdx = (*m_flxSolDep)[currFlxIdx][iSol];

      *(m_cellStatesFlxPnt[iFlxPnt]) += (*m_solPolyValsAtFlxPnts)[currFlxIdx][solIdx]*(*((*m_cellStates)[solIdx]));
      
      for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
      {
        *(m_cellGradFlxPnt[iFlxPnt][iVar]) += (*m_solPolyValsAtFlxPnts)[currFlxIdx][solIdx]*((*(m_cellGrads[solIdx]))[iVar]);
      }
    }
  }
  
  SafePtr< NavierStokesVarSet > navierStokesVarSet = m_diffusiveVarSet.d_castTo< NavierStokesVarSet >();
  
  // Loop over flux points to check the transition criterion
  for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
  {
    const CFreal rho = navierStokesVarSet->getDensity(*(m_cellStatesFlxPnt[iFlxPnt]));
    
    const CFreal muTot = navierStokesVarSet->getDynViscosity(*(m_cellStatesFlxPnt[iFlxPnt]), m_cellGradFlxPnt[iFlxPnt]);
    
//    const CFreal tau = muTot*(m_unitNormalFlxPnts[iFlxPnt][YY]*MathFunctions::innerProd(*(m_cellGradFlxPnt[iFlxPnt][1]),m_unitNormalFlxPnts[iFlxPnt]) -
//                       m_unitNormalFlxPnts[iFlxPnt][XX]*MathFunctions::innerProd(*(m_cellGradFlxPnt[iFlxPnt][2]),m_unitNormalFlxPnts[iFlxPnt]));
    
    CFreal tau = 0.0;
    
    if (m_dim == 2)
    {
      tau = muTot*(m_unitNormalFlxPnts[iFlxPnt][YY]*MathFunctions::innerProd(*(m_cellGradFlxPnt[iFlxPnt][1]),m_unitNormalFlxPnts[iFlxPnt]) -
            m_unitNormalFlxPnts[iFlxPnt][XX]*MathFunctions::innerProd(*(m_cellGradFlxPnt[iFlxPnt][2]),m_unitNormalFlxPnts[iFlxPnt]));
    }
    else if (fabs(m_unitNormalFlxPnts[iFlxPnt][ZZ]) <= fabs(m_unitNormalFlxPnts[iFlxPnt][XX]))
    {
      const CFreal tauT1 = muTot*(m_unitNormalFlxPnts[iFlxPnt][YY]*MathFunctions::innerProd(*(m_cellGradFlxPnt[iFlxPnt][1]),m_unitNormalFlxPnts[iFlxPnt]) -
                           m_unitNormalFlxPnts[iFlxPnt][XX]*MathFunctions::innerProd(*(m_cellGradFlxPnt[iFlxPnt][2]),m_unitNormalFlxPnts[iFlxPnt]));
      
      const CFreal tauT2 = muTot*(m_unitNormalFlxPnts[iFlxPnt][XX]*m_unitNormalFlxPnts[iFlxPnt][ZZ]*MathFunctions::innerProd(*(m_cellGradFlxPnt[iFlxPnt][1]),m_unitNormalFlxPnts[iFlxPnt]) +
                           m_unitNormalFlxPnts[iFlxPnt][YY]*m_unitNormalFlxPnts[iFlxPnt][ZZ]*MathFunctions::innerProd(*(m_cellGradFlxPnt[iFlxPnt][2]),m_unitNormalFlxPnts[iFlxPnt]) - 
                           (m_unitNormalFlxPnts[iFlxPnt][XX]*m_unitNormalFlxPnts[iFlxPnt][XX]+m_unitNormalFlxPnts[iFlxPnt][YY]*m_unitNormalFlxPnts[iFlxPnt][YY])*MathFunctions::innerProd(*(m_cellGradFlxPnt[iFlxPnt][3]),m_unitNormalFlxPnts[iFlxPnt]));
      tau = sqrt(tauT1*tauT1+tauT2*tauT2);
    }
    else
    {
      const CFreal tauT1 = muTot*(-m_unitNormalFlxPnts[iFlxPnt][ZZ]*MathFunctions::innerProd(*(m_cellGradFlxPnt[iFlxPnt][2]),m_unitNormalFlxPnts[iFlxPnt]) +
                           m_unitNormalFlxPnts[iFlxPnt][YY]*MathFunctions::innerProd(*(m_cellGradFlxPnt[iFlxPnt][3]),m_unitNormalFlxPnts[iFlxPnt]));
      
      const CFreal tauT2 = muTot*(-m_unitNormalFlxPnts[iFlxPnt][XX]*m_unitNormalFlxPnts[iFlxPnt][YY]*MathFunctions::innerProd(*(m_cellGradFlxPnt[iFlxPnt][2]),m_unitNormalFlxPnts[iFlxPnt]) -
                           m_unitNormalFlxPnts[iFlxPnt][XX]*m_unitNormalFlxPnts[iFlxPnt][ZZ]*MathFunctions::innerProd(*(m_cellGradFlxPnt[iFlxPnt][3]),m_unitNormalFlxPnts[iFlxPnt]) + 
                           (m_unitNormalFlxPnts[iFlxPnt][YY]*m_unitNormalFlxPnts[iFlxPnt][YY]+m_unitNormalFlxPnts[iFlxPnt][ZZ]*m_unitNormalFlxPnts[iFlxPnt][ZZ])*MathFunctions::innerProd(*(m_cellGradFlxPnt[iFlxPnt][1]),m_unitNormalFlxPnts[iFlxPnt]));
      tau = sqrt(tauT1*tauT1+tauT2*tauT2);
    }
    
    const CFreal tauCrit = tau/sqrt(rho*muTot);
    
    if (tauCrit <= (*(m_cellStatesFlxPnt[iFlxPnt]))[7])
    {
      m_bcStateComputer->setTransitionCriterion(iFlxPnt,true);
    }
    else
    {
      m_bcStateComputer->setTransitionCriterion(iFlxPnt,false); 
    }
  }
  
  // compute ghost states
  m_bcStateComputer->computeGhostStates(m_cellStatesFlxPnt,m_flxPntGhostSol,m_unitNormalFlxPnts,m_flxPntCoords);
   
  // compute ghost gradients
  m_bcStateComputer->computeGhostGradients(m_cellGradFlxPnt,m_flxPntGhostGrads,m_unitNormalFlxPnts,m_flxPntCoords);
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluxReconstructionMethod

} // namespace COOLFluiD
