#ifndef COOLFluiD_Numerics_FluxReconstructionMethod_BCDirichlet_hh
#define COOLFluiD_Numerics_FluxReconstructionMethod_BCDirichlet_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/BaseMethodStrategyProvider.hh"
#include "Framework/VarSetTransformer.hh"
#include "Framework/VectorialFunction.hh"
#include "FluxReconstructionMethod/BCStateComputer.hh"
#include "FluxReconstructionMethod/FluxReconstructionSolverData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a Dirichlet boundary condition
 *
 * @author Kris Van den Abeele
 * @author Ray Vandenhoeck
 */
class BCDirichlet : public BCStateComputer {

public:  // methods

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /// Constructor
  BCDirichlet(const std::string& name);

  /// Destructor
  ~BCDirichlet();

  /// Gets the Class name
  static std::string getClassName()
  {
    return "BCDirichlet";
  }

  /// Setup private data
  void setup();

  /// Unsetup private data
  void unsetup();

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

protected: // data

  /// physical model var set
  Common::SafePtr<Framework::ConvectiveVarSet> m_varSet;

  /// Transformer from Update to Linear Variables
  Common::SelfRegistPtr<Framework::VarSetTransformer> m_inputToUpdateVar;

  /// a vector of string to hold the functions
  std::vector<std::string> m_functions;

  /// a vector of for the variable names
  std::vector<std::string> m_vars;

  /// the VectorialFunction to use to parse the user conditions
  Framework::VectorialFunction m_vFunction;

  /// a string to hold the name of the input variables
  std::string m_inputVarStr;

  /// input state
  Framework::State* m_inputState;
  
  /// vector for space coords and time
  RealVector m_spaceTime;
  
  /// dimensional state
  Framework::State* m_dimState;

}; // class BCDirichlet

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod

}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif  // COOLFluiD_Numerics_FluxReconstructionMethod_BCDirichlet_hh

