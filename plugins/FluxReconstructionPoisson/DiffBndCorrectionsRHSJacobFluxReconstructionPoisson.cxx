#include "Framework/MethodCommandProvider.hh"

#include "Framework/MeshData.hh"
#include "Framework/BaseTerm.hh"

#include "MathTools/MathFunctions.hh"

#include "FluxReconstructionPoisson/DiffBndCorrectionsRHSJacobFluxReconstructionPoisson.hh"
#include "FluxReconstructionPoisson/FluxReconstructionPoisson.hh"
#include "FluxReconstructionMethod/FluxReconstructionElementData.hh"

#include "Poisson/PoissonDiffVarSet.hh"
#include "Poisson/PoissonConvVarSet.hh"



//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Physics::Poisson;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider< DiffBndCorrectionsRHSJacobFluxReconstructionPoisson, 
		       FluxReconstructionSolverData, 
		       FluxReconstructionPoissonModule >
DiffBndCorrectionsRHSJacobPoissonFluxReconstructionProvider("DiffBndCorrectionsRHSJacobPoisson");

//////////////////////////////////////////////////////////////////////////////

DiffBndCorrectionsRHSJacobFluxReconstructionPoisson::DiffBndCorrectionsRHSJacobFluxReconstructionPoisson(const std::string& name) :
  DiffBndCorrectionsRHSJacobFluxReconstruction(name)
{
}

//////////////////////////////////////////////////////////////////////////////

DiffBndCorrectionsRHSJacobFluxReconstructionPoisson::~DiffBndCorrectionsRHSJacobFluxReconstructionPoisson()
{
}

//////////////////////////////////////////////////////////////////////////////

void DiffBndCorrectionsRHSJacobFluxReconstructionPoisson::setup()
{
  DiffBndCorrectionsRHSJacobFluxReconstruction::setup();
  
  // get the diffusive varset
  m_diffVarSetPoisson = m_diffusiveVarSet.d_castTo< Physics::Poisson::PoissonDiffVarSet >();
  cf_assert(m_diffVarSetPoisson.isNotNull());
  m_convVarSetPoisson = getMethodData().getUpdateVar().d_castTo<Physics::Poisson::PoissonConvVarSet>();
}

//////////////////////////////////////////////////////////////////////////////

void DiffBndCorrectionsRHSJacobFluxReconstructionPoisson::unsetup()
{
  DiffBndCorrectionsRHSJacobFluxReconstruction::unsetup();
}
//////////////////////////////////////////////////////////////////////////////

std::vector< Common::SafePtr< BaseDataSocketSink > >
    DiffBndCorrectionsRHSJacobFluxReconstructionPoisson::needsSockets()
{
  std::vector< Common::SafePtr< BaseDataSocketSink > > result = DiffBndCorrectionsRHSJacobFluxReconstruction::needsSockets();

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

void DiffBndCorrectionsRHSJacobFluxReconstructionPoisson::prepareFluxComputation()
{
  //const bool isPerturb = this->getMethodData().isPerturb();
  //const CFuint iPerturbVar = this->getMethodData().iPerturbVar();

  //m_diffVarSetPoisson->setComposition(m_avgSol, isPerturb, iPerturbVar);
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluxReconstructionMethod

} // namespace COOLFluiD
