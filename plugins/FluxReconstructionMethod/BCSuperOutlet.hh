#ifndef COOLFluiD_Numerics_FluxReconstructionMethod_BCSuperOutlet_hh
#define COOLFluiD_Numerics_FluxReconstructionMethod_BCSuperOutlet_hh

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
 * @author Kris Van den Abeele
 * @author Ray Vandenhoeck
 * @author Rayan Dhib
 */
class BCSuperOutlet : public BCStateComputer {

public:  // methods

  /// Constructor
  BCSuperOutlet(const std::string& name);

  /// Destructor
  ~BCSuperOutlet();

  /// Gets the Class name
  static std::string getClassName()
  {
    return "BCSuperOutlet";
  }

  /// Set up private data
  void setup();
  
  /// Unset up private data
  void unsetup();
  
  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);
  
  /**
   * Configures the command.
   */
  void configure ( Config::ConfigArgs& args );

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
   * extrapolated to the flux points.
   */
  void computeBndGradVars(const std::vector< RealVector* >& gradVarsFlxPnt,
                          const std::vector< Framework::State* >& intStates,
                          const std::vector< Framework::State* >& ghostStates,
                          const std::vector< RealVector >& unitNormals,
                          const std::vector< RealVector >& flxPntCoords,
                          std::vector< RealVector* >& bndGradVars);

  /**
   * Sets the boundary states: the interior states, U_b = U.
   */
  void computeBndStates(const std::vector< Framework::State* >& intStates,
                        const std::vector< Framework::State* >& ghostStates,
                        const std::vector< RealVector >& unitNormals,
                        const std::vector< RealVector >& flxPntCoords,
                        std::vector< RealVector* >& bndStates);

  /**
   * Sets the boundary gradients: the compact face gradients, q_b = q, with the
   * normal component removed when ZeroGrad is true.
   */
  void computeBndGrads(const std::vector< std::vector< RealVector* > >& intGrads,
                       std::vector< std::vector< RealVector* > >& bndGrads,
                       const std::vector< RealVector* >& bndStates,
                       const std::vector< RealVector >& unitNormals,
                       const std::vector< RealVector >& flxPntCoords);
  
  protected: // data
      
  // boolean telling whether to put the normal gradients to zero
  bool m_zeroGrad;

}; // class BCSuperOutlet

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod

}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif  // COOLFluiD_Numerics_FluxReconstructionMethod_BCSuperOutlet_hh

