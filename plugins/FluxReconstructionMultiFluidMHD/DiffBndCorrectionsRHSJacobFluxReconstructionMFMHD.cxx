#include "Framework/MethodCommandProvider.hh"

#include "Framework/MeshData.hh"
#include "Framework/BaseTerm.hh"

#include "MathTools/MathFunctions.hh"

#include "FluxReconstructionMultiFluidMHD/DiffBndCorrectionsRHSJacobFluxReconstructionMFMHD.hh"
#include "FluxReconstructionMultiFluidMHD/FluxReconstructionMultiFluidMHD.hh"
#include "FluxReconstructionMethod/FluxReconstructionElementData.hh"
#include "MultiFluidMHD/MultiFluidMHDVarSet.hh"
#include "MultiFluidMHD/DiffMFMHDVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::MultiFluidMHD;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider< DiffBndCorrectionsRHSJacobFluxReconstructionMFMHD, 
		       FluxReconstructionSolverData, 
		       FluxReconstructionMultiFluidMHDModule >
DiffBndCorrectionsRHSJacobMFMHDFluxReconstructionProvider("DiffBndCorrectionsRHSJacobMFMHD");

//////////////////////////////////////////////////////////////////////////////

DiffBndCorrectionsRHSJacobFluxReconstructionMFMHD::DiffBndCorrectionsRHSJacobFluxReconstructionMFMHD(const std::string& name) :
  DiffBndCorrectionsRHSJacobFluxReconstruction(name),
  m_tempStates(),
  m_tempStatesSol()
{
}

//////////////////////////////////////////////////////////////////////////////

DiffBndCorrectionsRHSJacobFluxReconstructionMFMHD::~DiffBndCorrectionsRHSJacobFluxReconstructionMFMHD()
{
}

//////////////////////////////////////////////////////////////////////////////

void DiffBndCorrectionsRHSJacobFluxReconstructionMFMHD::setup()
{
  DiffBndCorrectionsRHSJacobFluxReconstruction::setup();
    
  m_tempStates.resize(2);

  m_tempStates[LEFT].resize(m_nbrFaceFlxPnts);
  m_tempStates[RIGHT].resize(m_nbrFaceFlxPnts);
  m_tempStatesSol.resize(m_nbrSolPnts);
}

//////////////////////////////////////////////////////////////////////////////

void DiffBndCorrectionsRHSJacobFluxReconstructionMFMHD::unsetup()
{
  DiffBndCorrectionsRHSJacobFluxReconstruction::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

void DiffBndCorrectionsRHSJacobFluxReconstructionMFMHD::computeWaveSpeedUpdates(CFreal& waveSpeedUpd)
{
  CFreal visc = 1.0;
  /// @todo needs to be changed for non-MFMHD
  SafePtr< DiffMFMHDVarSet > diffMFMHDVarSet = m_diffusiveVarSet.d_castTo< DiffMFMHDVarSet >();
  const CFreal dynVisc = diffMFMHDVarSet->getCurrDynViscosity();
  
  waveSpeedUpd = 0.0;
  for (CFuint iFlx = 0; iFlx < m_cellFlx.size(); ++iFlx)
  {
    const CFreal jacobXJacobXIntCoef = m_faceJacobVecAbsSizeFlxPnts[iFlx]*
                                 m_faceJacobVecAbsSizeFlxPnts[iFlx]*
                                   (*m_faceIntegrationCoefs)[iFlx]*
                                   m_cflConvDiffRatio;
    const CFreal rho = diffMFMHDVarSet->getDensity(*m_cellStatesFlxPnt[iFlx]);
    visc = dynVisc/rho;
				   
    // transform update states to physical data to calculate eigenvalues
    waveSpeedUpd += visc*jacobXJacobXIntCoef/m_cellVolume;
  }

}

//////////////////////////////////////////////////////////////////////////////

void DiffBndCorrectionsRHSJacobFluxReconstructionMFMHD::computeBndGradTerms(RealMatrix& gradTerm, RealMatrix& ghostGradTerm)
{ 
  SafePtr< DiffMFMHDVarSet > diffMFMHDVarSet = m_diffusiveVarSet.d_castTo< DiffMFMHDVarSet >();

  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    m_tempStates[LEFT][iFlx] = (m_cellStatesFlxPnt[iFlx]->getData());
    m_tempStates[RIGHT][iFlx] = (m_flxPntGhostSol[iFlx]->getData());
  }
  
  diffMFMHDVarSet->setGradientVars(m_tempStates[LEFT],gradTerm,m_nbrFaceFlxPnts);
  diffMFMHDVarSet->setGradientVars(m_tempStates[RIGHT],ghostGradTerm,m_nbrFaceFlxPnts);
}

//////////////////////////////////////////////////////////////////////////////

void DiffBndCorrectionsRHSJacobFluxReconstructionMFMHD::computeBndGradTerms2(RealMatrix& gradTerm, RealMatrix& ghostGradTerm)
{ 
  SafePtr< DiffMFMHDVarSet > diffMFMHDVarSet = m_diffusiveVarSet.d_castTo< DiffMFMHDVarSet >();

  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    m_tempStates[LEFT][iFlx] = (m_cellStatesFlxPnt2[iFlx]->getData());
    m_tempStates[RIGHT][iFlx] = (m_flxPntGhostSol[iFlx]->getData());
  }
  
  diffMFMHDVarSet->setGradientVars(m_tempStates[LEFT],gradTerm,m_nbrFaceFlxPnts);
  diffMFMHDVarSet->setGradientVars(m_tempStates[RIGHT],ghostGradTerm,m_nbrFaceFlxPnts);
}

//////////////////////////////////////////////////////////////////////////////

void DiffBndCorrectionsRHSJacobFluxReconstructionMFMHD::computeCellGradTerm(RealMatrix& gradTerm)
{
  SafePtr< DiffMFMHDVarSet > diffMFMHDVarSet = m_diffusiveVarSet.d_castTo< DiffMFMHDVarSet >();
  
  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {
    m_tempStatesSol[iSol] = ((*m_cellStates)[iSol]->getData());
  }
  
  diffMFMHDVarSet->setGradientVars(m_tempStatesSol,gradTerm,m_nbrSolPnts);
}

//////////////////////////////////////////////////////////////////////////////

void DiffBndCorrectionsRHSJacobFluxReconstructionMFMHD::computeFaceGradTerms(RealMatrix& gradTermL, RealMatrix& gradTermR)
{
  SafePtr< DiffMFMHDVarSet > diffMFMHDVarSet = m_diffusiveVarSet.d_castTo< DiffMFMHDVarSet >();

  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    m_tempStates[LEFT][iFlx] = (m_pertCellStatesFlxPnt[LEFT][iFlx]->getData());
    m_tempStates[RIGHT][iFlx] = (m_pertCellStatesFlxPnt[RIGHT][iFlx]->getData());
  }
  
  diffMFMHDVarSet->setGradientVars(m_tempStates[LEFT],gradTermL,m_nbrFaceFlxPnts);
  diffMFMHDVarSet->setGradientVars(m_tempStates[RIGHT],gradTermR,m_nbrFaceFlxPnts);
}

//////////////////////////////////////////////////////////////////////////////

void DiffBndCorrectionsRHSJacobFluxReconstructionMFMHD::prepareFluxComputation()
{
  const bool isPerturb = this->getMethodData().isPerturb();
  const CFuint iPerturbVar = this->getMethodData().iPerturbVar();
  SafePtr< DiffMFMHDVarSet > diffMFMHDVarSet = m_diffusiveVarSet.d_castTo< DiffMFMHDVarSet >();
  diffMFMHDVarSet->setComposition(m_avgSol, isPerturb, iPerturbVar);
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluxReconstructionMethod

} // namespace COOLFluiD
