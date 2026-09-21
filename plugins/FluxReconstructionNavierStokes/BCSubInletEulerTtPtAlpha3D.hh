#ifndef COOLFluiD_Numerics_FluxReconstructionMethod_BCSubInletEulerTtPtAlpha3D_hh
#define COOLFluiD_Numerics_FluxReconstructionMethod_BCSubInletEulerTtPtAlpha3D_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/BaseMethodStrategyProvider.hh"
#include "FluxReconstructionMethod/BCStateComputer.hh"
#include "FluxReconstructionMethod/FluxReconstructionSolverData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {
    namespace NavierStokes {
      class Euler3DVarSet;
    }
  }

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a subsonic inlet boundary condition for the 3D Euler/Navier-Stokes equations
 *
 * @author Kris Van den Abeele
 * @author Rayan Dhib
 */
class BCSubInletEulerTtPtAlpha3D : public BCStateComputer {

public:  // methods

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /// Constructor
  BCSubInletEulerTtPtAlpha3D(const std::string& name);

  /// Destructor
  ~BCSubInletEulerTtPtAlpha3D();

  /// Gets the Class name
  static std::string getClassName()
  {
    return "BCSubInletEulerTtPtAlpha3D";
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
   * Sets the boundary values of the gradient variables: the (p, u, v, w, T) of
   * computeInletPrimState.
   */
  void computeBndGradVars(const std::vector< RealVector* >& gradVarsFlxPnt,
                          const std::vector< Framework::State* >& intStates,
                          const std::vector< Framework::State* >& ghostStates,
                          const std::vector< RealVector >& unitNormals,
                          const std::vector< RealVector >& flxPntCoords,
                          std::vector< RealVector* >& bndGradVars);

  /**
   * Sets the boundary states from the (p, u, v, w, T) of computeInletPrimState.
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
   * Computes the boundary values of (p, u, v, w, T) from the total temperature and
   * pressure, the flow angles and the Mach number M of the interior state:
   *
   *   c = 1 + 0.5 (gamma-1) M^2,   T = Tt/c,   p = pt/c^(gamma/(gamma-1)),
   *   |u| = M sqrt(gamma R T),   v = tan(alphaXY) u, w = tan(alphaXZ) u
   */
  void computeInletPrimState(const Framework::State& intState,
                             RealVector& primState);

protected: // data

  /// physical model (in conservative variables)
  Common::SafePtr<Physics::NavierStokes::Euler3DVarSet> m_eulerVarSet;

  /// variable for physical data of ghostSol
  RealVector m_ghostSolPhysData;

  /// variable for physical data of intSol
  RealVector m_intSolPhysData;

  /// boundary values of the primitive variables (p, u, v, w, T)
  RealVector m_bndPrimState;

  /// total temperature
  CFreal     m_tTotal;

  /// total pressure
  CFreal     m_pTotal;

  /// alpha (for v/u (Y/X))
  CFreal     m_alphaXY;

  /// alpha (for w/u (Z/X))
  CFreal     m_alphaXZ;

}; // class BCSubInletEulerTtPtAlpha3D

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod

}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif  // COOLFluiD_Numerics_FluxReconstructionMethod_BCSubInletEulerTtPtAlpha3D_hh

