#include "Framework/MethodCommandProvider.hh"

#include "Framework/MeshData.hh"
#include "Framework/BaseTerm.hh"

#include "MathTools/MathFunctions.hh"

#include "FluxReconstructionNavierStokes/DiffBndCorrectionsRHSJacobFluxReconstructionNS.hh"
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

MethodCommandProvider< DiffBndCorrectionsRHSJacobFluxReconstructionNS, 
		       FluxReconstructionSolverData, 
		       FluxReconstructionNavierStokesModule >
DiffBndCorrectionsRHSJacobNSFluxReconstructionProvider("DiffBndCorrectionsRHSJacobNS");

//////////////////////////////////////////////////////////////////////////////

DiffBndCorrectionsRHSJacobFluxReconstructionNS::DiffBndCorrectionsRHSJacobFluxReconstructionNS(const std::string& name) :
  DiffBndCorrectionsRHSJacobFluxReconstruction(name),
  m_tempStates(),
  m_tempStatesSol()
{
}

//////////////////////////////////////////////////////////////////////////////

DiffBndCorrectionsRHSJacobFluxReconstructionNS::~DiffBndCorrectionsRHSJacobFluxReconstructionNS()
{
}

//////////////////////////////////////////////////////////////////////////////

void DiffBndCorrectionsRHSJacobFluxReconstructionNS::setup()
{
  DiffBndCorrectionsRHSJacobFluxReconstruction::setup();
    
  m_tempStates.resize(2);

  m_tempStates[LEFT].resize(m_nbrFaceFlxPnts);
  m_tempStates[RIGHT].resize(m_nbrFaceFlxPnts);
  m_tempStatesSol.resize(m_nbrSolPnts);
}

//////////////////////////////////////////////////////////////////////////////

void DiffBndCorrectionsRHSJacobFluxReconstructionNS::unsetup()
{
  DiffBndCorrectionsRHSJacobFluxReconstruction::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

void DiffBndCorrectionsRHSJacobFluxReconstructionNS::computeWaveSpeedUpdates(CFreal& waveSpeedUpd)
{
  CFreal visc = 1.0;
  /// @todo needs to be changed for non-NS
  SafePtr< NavierStokesVarSet > navierStokesVarSet = m_diffusiveVarSet.d_castTo< NavierStokesVarSet >();
  const CFreal dynVisc = navierStokesVarSet->getCurrDynViscosity();
  
  const CFreal factorPr = min(navierStokesVarSet->getModel().getPrandtl(),1.0);
  cf_assert(factorPr>0.0);
  
  waveSpeedUpd = 0.0;
  //for (CFuint iFlx = 0; iFlx < m_cellFlx.size(); ++iFlx)
  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    const CFreal jacobXJacobXIntCoef = m_faceJacobVecAbsSizeFlxPnts[iFlx]*
                                 m_faceJacobVecAbsSizeFlxPnts[iFlx]*
                                   (*m_faceIntegrationCoefs)[iFlx]*
                                   m_cflConvDiffRatio;
    const CFreal rho = navierStokesVarSet->getDensity(*m_cellStatesFlxPnt[iFlx]);
    visc = dynVisc/rho/factorPr;
				   
    // transform update states to physical data to calculate eigenvalues
    waveSpeedUpd += visc*jacobXJacobXIntCoef/m_cellVolume;
  }

}

//////////////////////////////////////////////////////////////////////////////

void DiffBndCorrectionsRHSJacobFluxReconstructionNS::computeBndGradTerms(RealMatrix& gradTerm, RealMatrix& ghostGradTerm)
{ 
  SafePtr< NavierStokesVarSet > navierStokesVarSet = m_diffusiveVarSet.d_castTo< NavierStokesVarSet >();

  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    m_tempStates[LEFT][iFlx] = (m_cellStatesFlxPnt[iFlx]->getData());
    m_tempStates[RIGHT][iFlx] = (m_flxPntGhostSol[iFlx]->getData());
  }
  
  navierStokesVarSet->setGradientVars(m_tempStates[LEFT],gradTerm,m_nbrFaceFlxPnts);
  navierStokesVarSet->setGradientVars(m_tempStates[RIGHT],ghostGradTerm,m_nbrFaceFlxPnts);
}

//////////////////////////////////////////////////////////////////////////////

void DiffBndCorrectionsRHSJacobFluxReconstructionNS::computeBndGradTerms2(RealMatrix& gradTerm, RealMatrix& ghostGradTerm)
{ 
  SafePtr< NavierStokesVarSet > navierStokesVarSet = m_diffusiveVarSet.d_castTo< NavierStokesVarSet >();

  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    m_tempStates[LEFT][iFlx] = (m_cellStatesFlxPnt2[iFlx]->getData());
    m_tempStates[RIGHT][iFlx] = (m_flxPntGhostSol[iFlx]->getData());
  }
  
  navierStokesVarSet->setGradientVars(m_tempStates[LEFT],gradTerm,m_nbrFaceFlxPnts);
  navierStokesVarSet->setGradientVars(m_tempStates[RIGHT],ghostGradTerm,m_nbrFaceFlxPnts);
}

//////////////////////////////////////////////////////////////////////////////

void DiffBndCorrectionsRHSJacobFluxReconstructionNS::computeCellGradTerm(RealMatrix& gradTerm)
{
  SafePtr< NavierStokesVarSet > navierStokesVarSet = m_diffusiveVarSet.d_castTo< NavierStokesVarSet >();
  
  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {
    m_tempStatesSol[iSol] = ((*m_cellStates)[iSol]->getData());
  }
  
  navierStokesVarSet->setGradientVars(m_tempStatesSol,gradTerm,m_nbrSolPnts);
}

//////////////////////////////////////////////////////////////////////////////

void DiffBndCorrectionsRHSJacobFluxReconstructionNS::computeFaceGradTerms(RealMatrix& gradTermL, RealMatrix& gradTermR)
{
  SafePtr< NavierStokesVarSet > navierStokesVarSet = m_diffusiveVarSet.d_castTo< NavierStokesVarSet >();

  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    m_tempStates[LEFT][iFlx] = (m_pertCellStatesFlxPnt[LEFT][iFlx]->getData());
    m_tempStates[RIGHT][iFlx] = (m_pertCellStatesFlxPnt[RIGHT][iFlx]->getData());
  }
  
  navierStokesVarSet->setGradientVars(m_tempStates[LEFT],gradTermL,m_nbrFaceFlxPnts);
  navierStokesVarSet->setGradientVars(m_tempStates[RIGHT],gradTermR,m_nbrFaceFlxPnts);
}

//////////////////////////////////////////////////////////////////////////////

void DiffBndCorrectionsRHSJacobFluxReconstructionNS::prepareFluxComputation()
{
  const bool isPerturb = this->getMethodData().isPerturb();
  const CFuint iPerturbVar = this->getMethodData().iPerturbVar();
  SafePtr< NavierStokesVarSet > navierStokesVarSet = m_diffusiveVarSet.d_castTo< NavierStokesVarSet >();
  navierStokesVarSet->setComposition(m_avgSol, isPerturb, iPerturbVar);
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluxReconstructionMethod

} // namespace COOLFluiD
