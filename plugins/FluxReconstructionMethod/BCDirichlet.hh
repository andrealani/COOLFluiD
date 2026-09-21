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
 * @author Rayan Dhib
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

  /**
   * Sets the boundary values of the gradient variables: the gradient variables
   * of the prescribed state at the flux point coordinates and the current time.
   */
  void computeBndGradVars(const std::vector< RealVector* >& gradVarsFlxPnt,
                          const std::vector< Framework::State* >& intStates,
                          const std::vector< Framework::State* >& ghostStates,
                          const std::vector< RealVector >& unitNormals,
                          const std::vector< RealVector >& flxPntCoords,
                          std::vector< RealVector* >& bndGradVars);

  /**
   * Sets the boundary states: the average of the interior and the mirrored
   * ghost state, that is the prescribed state.
   */
  void computeBndStates(const std::vector< Framework::State* >& intStates,
                        const std::vector< Framework::State* >& ghostStates,
                        const std::vector< RealVector >& unitNormals,
                        const std::vector< RealVector >& flxPntCoords,
                        std::vector< RealVector* >& bndStates);

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

  /// prescribed states in update variables at the flux points
  std::vector< RealVector > m_prescStates;

  /// pointers to m_prescStates, passed to the diffusive variable set
  std::vector< RealVector* > m_prescStatePtrs;

  /// gradient variables of the prescribed states
  RealMatrix m_prescGradVars;

}; // class BCDirichlet

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod

}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif  // COOLFluiD_Numerics_FluxReconstructionMethod_BCDirichlet_hh

