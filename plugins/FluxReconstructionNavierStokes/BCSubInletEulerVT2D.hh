#ifndef COOLFluiD_FluxReconstructionMethod_BCSubInletEulerVT2D_hh
#define COOLFluiD_FluxReconstructionMethod_BCSubInletEulerVT2D_hh

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
 * This class represents a subsonic inlet boundary condition for the 2D Euler/Navier-Stokes equations
 * with imposed velocity vector and temperature
 *
 * @author Ray Vandenhoeck
 * @author Rayan Dhib
 */
class BCSubInletEulerVT2D : public BCStateComputer {

public:  // methods

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /// Constructor
  BCSubInletEulerVT2D(const std::string& name);

  /// Destructor
  ~BCSubInletEulerVT2D();

  /// Gets the Class name
  static std::string getClassName()
  {
    return "BCSubInletEulerVT2D";
  }

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
   * extrapolated to the flux points with the prescribed velocity and temperature.
   */
  void computeBndGradVars(const std::vector< RealVector* >& gradVarsFlxPnt,
                          const std::vector< Framework::State* >& intStates,
                          const std::vector< Framework::State* >& ghostStates,
                          const std::vector< RealVector >& unitNormals,
                          const std::vector< RealVector >& flxPntCoords,
                          std::vector< RealVector* >& bndGradVars);

  /**
   * Sets the boundary states from the (p, u, v, T) of computeInletPrimState.
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

protected: // functions

  /**
   * Computes the boundary values of (p, u, v, T): the interior pressure and the
   * prescribed velocity and temperature.
   */
  void computeInletPrimState(const Framework::State& intState,
                             RealVector& primState);

protected: // data

  /// physical model (in conservative variables)
  Common::SafePtr<Physics::NavierStokes::Euler2DVarSet> m_eulerVarSet;

  /// variable for physical data of ghostSol
  RealVector m_ghostSolPhysData;

  /// variable for physical data of intSol
  RealVector m_intSolPhysData;

  /// boundary values of the primitive variables (p, u, v, T)
  RealVector m_bndPrimState;

  /// x-velocity
  CFreal     m_u;

  /// y-velocity
  CFreal     m_v;

  /// temperature
  CFreal     m_T;

}; // class BCSubInletEulerVT2D

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod

}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif  // COOLFluiD_Numerics_FluxReconstructionMethod_BCSubInletEulerVT2D_hh

