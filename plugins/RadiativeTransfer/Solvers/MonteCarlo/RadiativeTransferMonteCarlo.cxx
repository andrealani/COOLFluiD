#include "RadiativeTransfer/Solvers/MonteCarlo/RadiativeTransferMonteCarlo.hh"
#include "LagrangianSolver/ParticleTracking/ParticleTrackingAxi.hh"
#include "LagrangianSolver/ParticleTracking/ParticleTracking3D.hh"
#include "LagrangianSolver/ParticleTracking/ParticleTracking2D.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

namespace RadiativeTransfer {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<RadiativeTransferMonteCarlo<LagrangianSolver::ParticleTrackingAxi>,
                      DataProcessingData, 
		      RadiativeTransferModule>
RadiativeTransferMonteCarloAxiFVMCCProvider("RadiativeTransferMonteCarloAxiFVMCC");


MethodCommandProvider<RadiativeTransferMonteCarlo<LagrangianSolver::ParticleTracking3D>,
		      DataProcessingData,    
		      RadiativeTransferModule>
RadiativeTransferMonteCarlo3DFVMCCProvider("RadiativeTransferMonteCarlo3DFVMCC");

  
MethodCommandProvider<RadiativeTransferMonteCarlo<LagrangianSolver::ParticleTracking2D>,
		      DataProcessingData,    
RadiativeTransferModule>
RadiativeTransferMonteCarlo2DFVMCCProvider("RadiativeTransferMonteCarlo2DFVMCC");



//////////////////////////////////////////////////////////////////////////////

} // namespace RadiativeTransfer

} // namespace COOLFluiD

