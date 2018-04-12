#include "Framework/MethodStrategyProvider.hh"

#include "FluxReconstructionNavierStokes/AUSMPlusFlux.hh"
#include "FluxReconstructionNavierStokes/AUSMPlusUpFlux.hh"
#include "FluxReconstructionNavierStokes/FluxReconstructionNavierStokes.hh"

#include "NavierStokes/Euler2DVarSet.hh"
#include "NavierStokes/Euler3DVarSet.hh"

#include "NavierStokes/MultiScalarVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

// AUSM+
Framework::MethodStrategyProvider< AUSMPlusFlux<Euler2DVarSet>,
                                   FluxReconstructionSolverData,
				   RiemannFlux,
				   FluxReconstructionNavierStokesModule > 
				   AUSMPlusFlux2DProvider("AUSMPlusFlux2D");
				   
Framework::MethodStrategyProvider< AUSMPlusFlux<Euler3DVarSet>,
                                   FluxReconstructionSolverData,
				   RiemannFlux,
				   FluxReconstructionNavierStokesModule > 
				   AUSMPlusFlux3DProvider("AUSMPlusFlux3D");
				   
Framework::MethodStrategyProvider<AUSMPlusFlux <MultiScalarVarSet<Euler2DVarSet> >,
		       FluxReconstructionSolverData,
                       RiemannFlux,
		       FluxReconstructionNavierStokesModule>
AUSMPlusFluxMS2DProvider("AUSMPlusFluxMS2D");
    
Framework::MethodStrategyProvider<AUSMPlusFlux <MultiScalarVarSet<Euler3DVarSet> >,
		       FluxReconstructionSolverData,
                       RiemannFlux,
		       FluxReconstructionNavierStokesModule>
AUSMPlusFluxMS3DProvider("AUSMPlusFluxMS3D");
				   
// AUSM+Up
Framework::MethodStrategyProvider< AUSMPlusUpFlux<Euler2DVarSet>,
                                   FluxReconstructionSolverData,
				   RiemannFlux,
				   FluxReconstructionNavierStokesModule > 
				   AUSMPlusUpFlux2DProvider("AUSMPlusUpFlux2D");
				   
Framework::MethodStrategyProvider< AUSMPlusUpFlux<Euler3DVarSet>,
                                   FluxReconstructionSolverData,
				   RiemannFlux,
				   FluxReconstructionNavierStokesModule > 
				   AUSMPlusUpFlux3DProvider("AUSMPlusUpFlux3D");
    


//////////////////////////////////////////////////////////////////////////////

  } // namespace FluxReconstructionMethod

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
