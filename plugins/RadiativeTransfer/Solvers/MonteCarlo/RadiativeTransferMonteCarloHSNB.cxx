#include "RadiativeTransfer/Solvers/MonteCarlo/RadiativeTransferMonteCarloHSNB.hh"
#include "LagrangianSolver/ParticleTracking/ParticleTrackingAxi.hh"
#include "LagrangianSolver/ParticleTracking/ParticleTracking3D.hh"
#include "LagrangianSolver/ParticleTracking/ParticleTracking2D.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

namespace RadiativeTransfer {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<RadiativeTransferMonteCarloHSNB<LagrangianSolver::ParticleTrackingAxi>,
                      DataProcessingData,
              RadiativeTransferModule>
RadiativeTransferMonteCarloHSNBAxiFVMCCProvider("RadiativeTransferMonteCarloHSNBAxiFVMCC");


MethodCommandProvider<RadiativeTransferMonteCarloHSNB<LagrangianSolver::ParticleTracking3D>,
              DataProcessingData,
              RadiativeTransferModule>
RadiativeTransferMonteCarloHSNB3DFVMCCProvider("RadiativeTransferMonteCarloHSNB3DFVMCC");


MethodCommandProvider<RadiativeTransferMonteCarloHSNB<LagrangianSolver::ParticleTracking2D>,
              DataProcessingData,
RadiativeTransferModule>
RadiativeTransferMonteCarloHSNB2DFVMCCProvider("RadiativeTransferMonteCarloHSNB2DFVMCC");







//////////////////////////////////////////////////////////////////////////////

} // namespace RadiativeTransfer

} // namespace COOLFluiD

