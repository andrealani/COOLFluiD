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
 * This class represents a supersonic outlet boundary condition
 *
 * @author Ray Vandenhoeck
 * @author Firas Ben Ameur
 */
class BCMirrorVelocity : public BCStateComputer {

public:  // methods

  /// Constructor
  BCMirrorVelocity(const std::string& name);

  /// Destructor
  ~BCMirrorVelocity();

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
  
}; // class BCMirrorVelocity
    
    //////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod

}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif  // COOLFluiD_Numerics_FluxReconstructionMethod_BCMirrorVelocity_hh

