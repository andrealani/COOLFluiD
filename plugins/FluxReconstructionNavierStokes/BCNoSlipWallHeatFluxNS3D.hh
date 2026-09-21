#ifndef COOLFluiD_Numerics_FluxReconstructionMethod_BCNoSlipWallHeatFluxNS3D_hh
#define COOLFluiD_Numerics_FluxReconstructionMethod_BCNoSlipWallHeatFluxNS3D_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/BaseMethodStrategyProvider.hh"
#include "FluxReconstructionMethod/BCStateComputer.hh"
#include "FluxReconstructionMethod/FluxReconstructionSolverData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {
    namespace NavierStokes {
      class Euler3DVarSet;
      class NavierStokesVarSet;
    }
  }

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a no-slip-wall boundary condition with imposed heat flux or wall temperature
 * for the 3D Euler/Navier-Stokes equations.
 *
 * Ghost state, with u_i the interior state at the wall flux point (the trace of the cell polynomial):
 * interior density, velocity reversed (v_g = -v_i, fixed wall), temperature reflected about the wall
 * temperature on an isothermal wall (T_g = 2 T_wall - T_i, at least 0.01 T_wall) or copied on a
 * heat-flux wall, pressure from density and temperature. The average of u_i and the ghost carries the
 * wall velocity and temperature, and no mass crosses the wall.
 * LegacyGhost = true restores the previous ghost (the wall state on an isothermal wall).
 *
 * @author Ray Vandenhoeck
 * @author Rayan Dhib
 */
class BCNoSlipWallHeatFluxNS3D : public BCStateComputer {

public:  // methods

  /// Constructor
  BCNoSlipWallHeatFluxNS3D(const std::string& name);

  /// Destructor
  ~BCNoSlipWallHeatFluxNS3D();

  /// Gets the Class name
  static std::string getClassName()
  {
    return "BCNoSlipWallHeatFluxNS3D";
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

  /**
   * Sets the boundary values of the gradient variables: the gradient variables
   * extrapolated to the flux points with zero velocity and, for an isothermal wall,
   * the wall temperature.
   */
  void computeBndGradVars(const std::vector< RealVector* >& gradVarsFlxPnt,
                          const std::vector< Framework::State* >& intStates,
                          const std::vector< Framework::State* >& ghostStates,
                          const std::vector< RealVector >& unitNormals,
                          const std::vector< RealVector >& flxPntCoords,
                          std::vector< RealVector* >& bndGradVars);

  /**
   * Sets the boundary states from (p, u, v, w, T) with the interior pressure, zero velocity
   * and the wall temperature, or the interior temperature when the heat flux is
   * prescribed.
   */
  void computeBndStates(const std::vector< Framework::State* >& intStates,
                        const std::vector< Framework::State* >& ghostStates,
                        const std::vector< RealVector >& unitNormals,
                        const std::vector< RealVector >& flxPntCoords,
                        std::vector< RealVector* >& bndStates);

  /**
   * Sets the boundary gradients: the compact face gradients, q_b = q.
   */
  void computeBndGrads(const std::vector< std::vector< RealVector* > >& intGrads,
                       std::vector< std::vector< RealVector* > >& bndGrads,
                       const std::vector< RealVector* >& bndStates,
                       const std::vector< RealVector >& unitNormals,
                       const std::vector< RealVector >& flxPntCoords);

  /**
   * Sets the normal temperature gradient from the prescribed wall heat flux, see
   * prescribeNSWallHeatFlux. Only when the heat flux is prescribed.
   */
  void constrainBndGrads(const RealVector& bndState,
                         std::vector< RealVector* >& bndGrads,
                         const RealVector& unitNormal,
                         const RealVector& flxPntCoord);

protected: // data

  /// physical model (in conservative variables)
  Common::SafePtr<Physics::NavierStokes::Euler3DVarSet> m_eulerVarSet;

  /// variable for physical data of ghostSol
  RealVector m_ghostSolPhysData;

  /// variable for physical data of intSol
  RealVector m_intSolPhysData;

  /// boundary values of the primitive variables (p, u, v, w, T)
  RealVector m_bndPrimState;

  /// diffusive variable set
  Common::SafePtr< Physics::NavierStokes::NavierStokesVarSet > m_diffusiveVarSet;
  
  /// boolean telling if the wall has constant heat flux
  bool m_heatFlux;
  
  /// wall static temperature
  CFreal m_wallT;
  
  /// wall heat flux
  CFreal m_wallQ;

  /// iteration after which is changed to an isothermal wall BC
  CFuint m_changeToIsoT;

  /// use the previous ghost state instead of the reflected one
  bool m_legacyGhost;

}; // class BCNoSlipWallHeatFluxNS3D

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod

}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif  // COOLFluiD_Numerics_FluxReconstructionMethod_BCNoSlipWallHeatFluxNS3D_hh

