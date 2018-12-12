#ifndef COOLFluiD_Numerics_FluxReconstructionMethod_BCNoSlipWallHeatFluxNS2D_hh
#define COOLFluiD_Numerics_FluxReconstructionMethod_BCNoSlipWallHeatFluxNS2D_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/BaseMethodStrategyProvider.hh"
#include "FluxReconstructionMethod/BCStateComputer.hh"
#include "FluxReconstructionMethod/FluxReconstructionSolverData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {
    namespace NavierStokes {
      class Euler2DVarSet;
    }
  }

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a no-slip-wall boundary condition with imposed heat flux or wall temperature
 * for the 2D Euler/Navier-Stokes equations.
 *
 * @author Ray Vandenhoeck
 */
class BCNoSlipWallHeatFluxNS2D : public BCStateComputer {

public:  // methods

  /// Constructor
  BCNoSlipWallHeatFluxNS2D(const std::string& name);

  /// Destructor
  ~BCNoSlipWallHeatFluxNS2D();

  /// Gets the Class name
  static std::string getClassName()
  {
    return "BCNoSlipWallHeatFluxNS2D";
  }
  
  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /// Set up private data and data
  void setup();

  /**
   * Sets the ghost states in all the boundary points (depends on the boundary condition type)
   */
  void computeGhostStates(const std::vector< Framework::State* >& intStates,
                                  std::vector< Framework::State* >& ghostStates,
                                  const std::vector< RealVector >& normals,
                                  const std::vector< RealVector >& coords);

  /**
   * Sets the ghost gradients in all the boundary points (depends on the boundary condition type)
   */
  void computeGhostGradients(const std::vector< std::vector< RealVector* > >& intGrads,
                             std::vector< std::vector< RealVector* > >& ghostGrads,
                             const std::vector< RealVector >& normals,
                             const std::vector< RealVector >& coords);

protected: // data

  /// physical model (in conservative variables)
  Common::SafePtr<Physics::NavierStokes::Euler2DVarSet> m_eulerVarSet;

  /// variable for physical data of ghostSol
  RealVector m_ghostSolPhysData;

  /// variable for physical data of intSol
  RealVector m_intSolPhysData;
  
  /// boolean telling if the wall has constant heat flux
  bool m_heatFlux;
  
  /// boolean telling whether the strong ghost T should be used
  bool m_strongT;
  
  /// wall static temperature
  CFreal m_wallT;
  
  /// wall heat flux
  CFreal m_wallQ;
  
  /// iteration after which is changed to an isothermal wall BC
  CFuint m_changeToIsoT;
  
  /// wall x-velocity
  CFreal m_wallU;
  
  /// wall y-velocity
  CFreal m_wallV;

}; // class BCNoSlipWallHeatFluxNS2D

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod

}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif  // COOLFluiD_Numerics_FluxReconstructionMethod_BCNoSlipWallHeatFluxNS2D_hh

