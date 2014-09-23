#ifndef COOLFluiD_Numerics_FluctSplit_SuperInlet_hh
#define COOLFluiD_Numerics_FluctSplit_SuperInlet_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/VectorialFunction.hh"
#include "Framework/DataSocketSink.hh"
#include "MathTools/FunctionParser.hh"

#include "FluctSplit/FluctuationSplitData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/// This class is a Numerical command that implements a supersonic inlet
/// boundary condition for the Fluctuation Splitting method.
/// It has the possibility of only being applied if a certain
/// condition specified by a user defined expression is bigger than zero.
/// @author Tiago Quintino
/// @author Thomas Wuilbaut
class FluctSplit_API SuperInlet : public FluctuationSplitCom {
public:

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Constructor
  SuperInlet(const std::string& name);

  /// Default destructor
  ~SuperInlet();

  /// Configures the command.
  virtual void configure ( Config::ConfigArgs& args );

  /// Set up the member data
  virtual void setup();

  /// UnSet up private data and data of the aggregated classes
  /// in this command after processing phase
  virtual void unsetup();

  /// Returns the DataSocket's that this command needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

protected: // methods

  /// Execute on the current TRS
  void executeOnTrs();

  /// Setup the FunctionParser that will parse the expression for the
  /// condition to apply this boundary condition.
  /// @throw Common::ParserException if the condition string or variables are badly set
  void setCondition();

  /// Helper function to throw the exception with the error message.
  /// @param add a string to add to the error output
  /// @throw Common::ParserException if the condition string or variables are badly set
  void throwConditionException(const Common::CodeLocation& where, const std::string& add = std::string());

protected: // data

  /// the socket to the data handle of the rhs
  Framework::DataSocketSink<CFreal> socket_rhs;

  /// the socket to the data handle of the isUpdated flag
  Framework::DataSocketSink<bool> socket_isUpdated;

  /// the socket to the data handle of the state's
  Framework::DataSocketSink < Framework::State* , Framework::GLOBAL > socket_states;

  /// user defined string for the condition
  bool m_checkCondition;

  /// user defined string for the condition
  std::string m_conditionStr;

  /// condition to specify wether to apply the BC or not
  MathTools::FunctionParser m_condition;

  /// a vector of string to hold the functions
  std::vector<std::string> m_functions;

  /// a vector of string to hold the functions
  std::vector<std::string> m_vars;

  // the VectorialFunction to use
  Framework::VectorialFunction m_vFunction;

  /// Used for skipping the zeroing of the rhs for the lasts equations
  CFuint _nbEquationsToSkip;

  /// physical model var set
  Common::SafePtr<Framework::ConvectiveVarSet> _varSet;

  /// Transformer from Update to Linear Variables
  Common::SelfRegistPtr<Framework::VarSetTransformer> _inputToUpdateVar;

  /// RealVector holding the input variables
  Framework::State* _input;

  /// Flag to input directly adimensional values
  bool _inputAdimensionalValues;

  /// a string to hold the name of the input variables
  std::string _inputVarStr;

  /// a string to hold the name of the update variables
  std::string _updateVarStr;

  /// IDs of the variables that will be changed interactively
  std::vector<CFuint> _interVarIDs;

  /// Factor to multiply the selected InteractiveVarIDs (should be < 1)
  CFreal _interFactor;

}; // end of class SuperInlet

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_SuperInlet_hh
