#ifndef COOLFluiD_Numerics_FluxReconstructionMethod_BCMirrorVelocity_hh
#define COOLFluiD_Numerics_FluxReconstructionMethod_BCMirrorVelocity_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/BaseMethodStrategyProvider.hh"
#include "FluxReconstructionMethod/BCStateComputer.hh"
#include "FluxReconstructionMethod/FluxReconstructionSolverData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a slip wall (mirror) boundary condition: the ghost
 * state copies the interior state except for its normal velocity v_n.
 *
 * Default ghost: the wall state, u_ghost = u_int - v_n n/|n|^2. With
 * ReflectedGhost = true: the mirror image, u_ghost = u_int - 2 v_n n/|n|^2.
 * The boundary state of the diffusive terms is the wall state in both cases.
 *
 * @author Ray Vandenhoeck
 * @author Firas Ben Ameur
 * @author Rayan Dhib
 */
class BCMirrorVelocity : public BCStateComputer {

public:  // methods

  /// Constructor
  BCMirrorVelocity(const std::string& name);

  /// Destructor
  ~BCMirrorVelocity();

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Gets the Class name
  static std::string getClassName()
  {
    return "BCMirrorVelocity";
  }

  /// Set up private data
  void setup();
  
  /// Unset up private data
  void unsetup();

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
   * extrapolated to the flux points without the normal velocity.
   */
  void computeBndGradVars(const std::vector< RealVector* >& gradVarsFlxPnt,
                          const std::vector< Framework::State* >& intStates,
                          const std::vector< Framework::State* >& ghostStates,
                          const std::vector< RealVector >& unitNormals,
                          const std::vector< RealVector >& flxPntCoords,
                          std::vector< RealVector* >& bndGradVars);

  /**
   * Sets the boundary states: the ghost states, U_b = U_ghost. The ghost state
   * already carries the projected wall velocity.
   */
  void computeBndStates(const std::vector< Framework::State* >& intStates,
                        const std::vector< Framework::State* >& ghostStates,
                        const std::vector< RealVector >& unitNormals,
                        const std::vector< RealVector >& flxPntCoords,
                        std::vector< RealVector* >& bndStates);

  /**
   * Sets the boundary gradients with the slip wall rule of
   * BCStateComputer::setSlipWallBndGrads.
   */
  void computeBndGrads(const std::vector< std::vector< RealVector* > >& intGrads,
                       std::vector< std::vector< RealVector* > >& bndGrads,
                       const std::vector< RealVector* >& bndStates,
                       const std::vector< RealVector >& unitNormals,
                       const std::vector< RealVector >& flxPntCoords);

protected:
  /// array of flags telling if a variable is a velocity component
  std::valarray<bool> m_isVelocityComp;

  /// velocity IDs
  std::vector< CFuint > m_velocityIDs;
  
  /// tangent vector
  RealVector m_tangent;

  /// temporary velocity gradient
  RealVector m_velocityNGradI;

  /// temporary velocity gradient
  RealVector m_velocityTGradI;

  /// temporary velocity gradient
  RealVector m_velocityNGradG;
  
  /// temporary velocity gradient
  RealVector m_velocityTGradG;

  /// ghost velocity is the mirror image of the interior velocity (default false: the wall state)
  bool m_reflectedGhost;
  
}; // class BCMirrorVelocity
    
    //////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod

}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif  // COOLFluiD_Numerics_FluxReconstructionMethod_BCMirrorVelocity_hh

