#include "Framework/MethodCommandProvider.hh"

#include "Framework/MeshData.hh"
#include "Framework/BaseTerm.hh"

#include "MathTools/MathFunctions.hh"

#include "FluxReconstructionMHD/DiffBndCorrectionsRHSJacobFluxReconstructionMHD.hh"
#include "FluxReconstructionMHD/FluxReconstructionMHD.hh"
#include "FluxReconstructionMethod/FluxReconstructionElementData.hh"
#include "MHD/MHDProjectionDiffVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::MHD;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider< DiffBndCorrectionsRHSJacobFluxReconstructionMHD, 
		       FluxReconstructionSolverData, 
		       FluxReconstructionMHDModule >
DiffBndCorrectionsRHSJacobMHDFluxReconstructionProvider("DiffBndCorrectionsRHSJacobMHD");

//////////////////////////////////////////////////////////////////////////////

DiffBndCorrectionsRHSJacobFluxReconstructionMHD::DiffBndCorrectionsRHSJacobFluxReconstructionMHD(const std::string& name) :
  DiffBndCorrectionsRHSJacobFluxReconstruction(name),
  m_tempStates(),
  m_tempStatesSol()
{
}

//////////////////////////////////////////////////////////////////////////////

DiffBndCorrectionsRHSJacobFluxReconstructionMHD::~DiffBndCorrectionsRHSJacobFluxReconstructionMHD()
{
}

//////////////////////////////////////////////////////////////////////////////

void DiffBndCorrectionsRHSJacobFluxReconstructionMHD::setup()
{
  DiffBndCorrectionsRHSJacobFluxReconstruction::setup();
  
  // get the diffusive varset
  m_diffusiveVarSet = (getMethodData().getDiffusiveVar()).d_castTo< MHDProjectionDiffVarSet >();
    
  m_tempStates.resize(2);

  m_tempStates[LEFT].resize(m_nbrFaceFlxPnts);
  m_tempStates[RIGHT].resize(m_nbrFaceFlxPnts);
  m_tempStatesSol.resize(m_nbrSolPnts);
}

//////////////////////////////////////////////////////////////////////////////

void DiffBndCorrectionsRHSJacobFluxReconstructionMHD::unsetup()
{
  DiffBndCorrectionsRHSJacobFluxReconstruction::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

void DiffBndCorrectionsRHSJacobFluxReconstructionMHD::computeWaveSpeedUpdates(CFreal& waveSpeedUpd)
{
  CFreal visc = 1.0;
  /// @todo needs to be changed for non-MHD
  //SafePtr< NavierStokesVarSet > navierStokesVarSet = m_diffusiveVarSet.d_castTo< NavierStokesVarSet >();
  const CFreal dynVisc = m_diffusiveVarSet->getCurrDynViscosity();
  
  const CFreal factorPr = 1.0;//min(m_diffusiveVarSet->getModel().getPrandtl(),1.0);
  cf_assert(factorPr>0.0);
  
  waveSpeedUpd = 0.0;
  for (CFuint iFlx = 0; iFlx < m_cellFlx.size(); ++iFlx)
  {
    const CFreal jacobXJacobXIntCoef = m_faceJacobVecAbsSizeFlxPnts[iFlx]*
                                 m_faceJacobVecAbsSizeFlxPnts[iFlx]*
                                   (*m_faceIntegrationCoefs)[iFlx]*
                                   m_cflConvDiffRatio;
    const CFreal rho = m_diffusiveVarSet->getDensity(*m_cellStatesFlxPnt[iFlx]);
    visc = dynVisc/rho/factorPr;
				   
    // transform update states to physical data to calculate eigenvalues
    waveSpeedUpd += visc*jacobXJacobXIntCoef/m_cellVolume;
  }

}

//////////////////////////////////////////////////////////////////////////////

void DiffBndCorrectionsRHSJacobFluxReconstructionMHD::computeBndGradTerms(RealMatrix& gradTerm, RealMatrix& ghostGradTerm)
{ 
  //SafePtr< NavierStokesVarSet > navierStokesVarSet = m_diffusiveVarSet.d_castTo< NavierStokesVarSet >();

  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    m_tempStates[LEFT][iFlx] = (m_cellStatesFlxPnt[iFlx]->getData());
    m_tempStates[RIGHT][iFlx] = (m_flxPntGhostSol[iFlx]->getData());
  }
  
  m_diffusiveVarSet->setGradientVars(m_tempStates[LEFT],gradTerm,m_nbrFaceFlxPnts);
  m_diffusiveVarSet->setGradientVars(m_tempStates[RIGHT],ghostGradTerm,m_nbrFaceFlxPnts);
}

//////////////////////////////////////////////////////////////////////////////

void DiffBndCorrectionsRHSJacobFluxReconstructionMHD::computeBndGradTerms2(RealMatrix& gradTerm, RealMatrix& ghostGradTerm)
{ 
  //SafePtr< NavierStokesVarSet > navierStokesVarSet = m_diffusiveVarSet.d_castTo< NavierStokesVarSet >();

  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    m_tempStates[LEFT][iFlx] = (m_cellStatesFlxPnt2[iFlx]->getData());
    m_tempStates[RIGHT][iFlx] = (m_flxPntGhostSol[iFlx]->getData());
  }
  
  m_diffusiveVarSet->setGradientVars(m_tempStates[LEFT],gradTerm,m_nbrFaceFlxPnts);
  m_diffusiveVarSet->setGradientVars(m_tempStates[RIGHT],ghostGradTerm,m_nbrFaceFlxPnts);
}

//////////////////////////////////////////////////////////////////////////////

void DiffBndCorrectionsRHSJacobFluxReconstructionMHD::computeCellGradTerm(RealMatrix& gradTerm)
{
  //SafePtr< NavierStokesVarSet > navierStokesVarSet = m_diffusiveVarSet.d_castTo< NavierStokesVarSet >();
  
  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {
    m_tempStatesSol[iSol] = ((*m_cellStates)[iSol]->getData());
  }
  
  m_diffusiveVarSet->setGradientVars(m_tempStatesSol,gradTerm,m_nbrSolPnts);
}

//////////////////////////////////////////////////////////////////////////////

void DiffBndCorrectionsRHSJacobFluxReconstructionMHD::computeFaceGradTerms(RealMatrix& gradTermL, RealMatrix& gradTermR)
{
  //SafePtr< NavierStokesVarSet > navierStokesVarSet = m_diffusiveVarSet.d_castTo< NavierStokesVarSet >();

  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    m_tempStates[LEFT][iFlx] = (m_pertCellStatesFlxPnt[LEFT][iFlx]->getData());
    m_tempStates[RIGHT][iFlx] = (m_pertCellStatesFlxPnt[RIGHT][iFlx]->getData());
  }
  
  m_diffusiveVarSet->setGradientVars(m_tempStates[LEFT],gradTermL,m_nbrFaceFlxPnts);
  m_diffusiveVarSet->setGradientVars(m_tempStates[RIGHT],gradTermR,m_nbrFaceFlxPnts);
}

//////////////////////////////////////////////////////////////////////////////

void DiffBndCorrectionsRHSJacobFluxReconstructionMHD::prepareFluxComputation()
{
  const bool isPerturb = this->getMethodData().isPerturb();
  const CFuint iPerturbVar = this->getMethodData().iPerturbVar();
  //SafePtr< NavierStokesVarSet > navierStokesVarSet = m_diffusiveVarSet.d_castTo< NavierStokesVarSet >();
  m_diffusiveVarSet->setComposition(m_avgSol, isPerturb, iPerturbVar);
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluxReconstructionMethod

} // namespace COOLFluiD
