#include "Framework/MethodCommandProvider.hh"

#include "Framework/MeshData.hh"
#include "Framework/BaseTerm.hh"

#include "MathTools/MathFunctions.hh"

#include "FluxReconstructionPoisson/DiffBndCorrectionsRHSJacobFluxReconstructionPoisson.hh"
#include "FluxReconstructionNavierStokes/FluxReconstructionNavierStokes.hh"
#include "FluxReconstructionMethod/FluxReconstructionElementData.hh"
#include "NavierStokes/NavierStokesVarSet.hh"

#include "Poisson/PoissonDiffVarSet.hh"
#include "Poisson/PoissonConvVarSet.hh"



//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::NavierStokes;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Physics::Poisson;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider< DiffBndCorrectionsRHSJacobFluxReconstructionPoisson, 
		       FluxReconstructionSolverData, 
		       FluxReconstructionNavierStokesModule >
DiffBndCorrectionsRHSJacobPoissonFluxReconstructionProvider("DiffBndCorrectionsRHSJacobPoisson");

//////////////////////////////////////////////////////////////////////////////

DiffBndCorrectionsRHSJacobFluxReconstructionPoisson::DiffBndCorrectionsRHSJacobFluxReconstructionPoisson(const std::string& name) :
  DiffBndCorrectionsRHSJacobFluxReconstructionNS(name),
  m_tempStates(),
  m_tempStatesSol()
{
}

//////////////////////////////////////////////////////////////////////////////

DiffBndCorrectionsRHSJacobFluxReconstructionPoisson::~DiffBndCorrectionsRHSJacobFluxReconstructionPoisson()
{
}

//////////////////////////////////////////////////////////////////////////////

void DiffBndCorrectionsRHSJacobFluxReconstructionPoisson::setup()
{
  DiffBndCorrectionsRHSJacobFluxReconstructionNS::setup();
  
  m_tempStates.resize(2);

  m_tempStates[LEFT].resize(m_nbrFaceFlxPnts);
  m_tempStates[RIGHT].resize(m_nbrFaceFlxPnts);
  m_tempStatesSol.resize(m_nbrSolPnts);
  
  // get the diffusive varset
  m_diffVarSetPoisson = m_diffusiveVarSet.d_castTo< Physics::Poisson::PoissonDiffVarSet >();
  cf_assert(m_diffVarSetPoisson.isNotNull());
  m_convVarSetPoisson = getMethodData().getUpdateVar().d_castTo<Physics::Poisson::PoissonConvVarSet>();
}

//////////////////////////////////////////////////////////////////////////////

void DiffBndCorrectionsRHSJacobFluxReconstructionPoisson::unsetup()
{
  DiffBndCorrectionsRHSJacobFluxReconstructionNS::unsetup();
}
//////////////////////////////////////////////////////////////////////////////

std::vector< Common::SafePtr< BaseDataSocketSink > >
    DiffBndCorrectionsRHSJacobFluxReconstructionPoisson::needsSockets()
{
  std::vector< Common::SafePtr< BaseDataSocketSink > > result = DiffBndCorrectionsRHSJacobFluxReconstructionNS::needsSockets();

  //result.push_back(&socket_wallDistance);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void DiffBndCorrectionsRHSJacobFluxReconstructionPoisson::computeWaveSpeedUpdates(CFreal& waveSpeedUpd)
{
  CFreal visc = 33.0;
  /// @todo needs to be changed for non-MFMHD
//  const CFreal dynVisc = diffMFMHDVarSet->getCurrDynViscosity();
  
  waveSpeedUpd = 0.0;
  //for (CFuint iFlx = 0; iFlx < m_cellFlx.size(); ++iFlx)
  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    const CFreal jacobXJacobXIntCoef = m_faceJacobVecAbsSizeFlxPnts[iFlx]*
                                 m_faceJacobVecAbsSizeFlxPnts[iFlx]*
                                   (*m_faceIntegrationCoefs)[iFlx]*
                                   m_cflConvDiffRatio;
    //const CFreal rho = diffMFMHDVarSet->getDensity(*m_cellStatesFlxPnt[iFlx]);
    //visc = dynVisc/rho;
				   
    // transform update states to physical data to calculate eigenvalues
    waveSpeedUpd += visc*jacobXJacobXIntCoef/m_cellVolume;
  }

}

//////////////////////////////////////////////////////////////////////////////

void DiffBndCorrectionsRHSJacobFluxReconstructionPoisson::computeBndGradTerms(RealMatrix& gradTerm, RealMatrix& ghostGradTerm)
{ 
  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    m_tempStates[LEFT][iFlx] = (m_cellStatesFlxPnt[iFlx]->getData());
    m_tempStates[RIGHT][iFlx] = (m_flxPntGhostSol[iFlx]->getData());
  }
  
  m_diffVarSetPoisson->setGradientVars(m_tempStates[LEFT],gradTerm,m_nbrFaceFlxPnts);
  m_diffVarSetPoisson->setGradientVars(m_tempStates[RIGHT],ghostGradTerm,m_nbrFaceFlxPnts);
}

//////////////////////////////////////////////////////////////////////////////

void DiffBndCorrectionsRHSJacobFluxReconstructionPoisson::computeBndGradTerms2(RealMatrix& gradTerm, RealMatrix& ghostGradTerm)
{ 
  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    m_tempStates[LEFT][iFlx] = (m_cellStatesFlxPnt2[iFlx]->getData());
    m_tempStates[RIGHT][iFlx] = (m_flxPntGhostSol[iFlx]->getData());
  }
  
  m_diffVarSetPoisson->setGradientVars(m_tempStates[LEFT],gradTerm,m_nbrFaceFlxPnts);
  m_diffVarSetPoisson->setGradientVars(m_tempStates[RIGHT],ghostGradTerm,m_nbrFaceFlxPnts);
}

//////////////////////////////////////////////////////////////////////////////

void DiffBndCorrectionsRHSJacobFluxReconstructionPoisson::computeCellGradTerm(RealMatrix& gradTerm)
{  
  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {
    m_tempStatesSol[iSol] = ((*m_cellStates)[iSol]->getData());
  }
  
  m_diffVarSetPoisson->setGradientVars(m_tempStatesSol,gradTerm,m_nbrSolPnts);
}

//////////////////////////////////////////////////////////////////////////////

void DiffBndCorrectionsRHSJacobFluxReconstructionPoisson::computeFaceGradTerms(RealMatrix& gradTermL, RealMatrix& gradTermR)
{
  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    m_tempStates[LEFT][iFlx] = (m_pertCellStatesFlxPnt[LEFT][iFlx]->getData());
    m_tempStates[RIGHT][iFlx] = (m_pertCellStatesFlxPnt[RIGHT][iFlx]->getData());
  }
  
  m_diffVarSetPoisson->setGradientVars(m_tempStates[LEFT],gradTermL,m_nbrFaceFlxPnts);
  m_diffVarSetPoisson->setGradientVars(m_tempStates[RIGHT],gradTermR,m_nbrFaceFlxPnts);
}

//////////////////////////////////////////////////////////////////////////////

void DiffBndCorrectionsRHSJacobFluxReconstructionPoisson::prepareFluxComputation()
{
  //const bool isPerturb = this->getMethodData().isPerturb();
  //const CFuint iPerturbVar = this->getMethodData().iPerturbVar();

  //m_diffVarSetPoisson->setComposition(m_avgSol, isPerturb, iPerturbVar);
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluxReconstructionMethod

} // namespace COOLFluiD
