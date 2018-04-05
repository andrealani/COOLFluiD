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

    } // namespace FluxReconstructionMethod

} // namespace COOLFluiD
