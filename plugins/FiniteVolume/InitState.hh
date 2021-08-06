#ifndef COOLFluiD_Numerics_FiniteVolume_InitState_hh
#define COOLFluiD_Numerics_FiniteVolume_InitState_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/VectorialFunction.hh"
#include "Framework/DataSocketSink.hh"
#include "Framework/VarSetTransformer.hh"
#include "FiniteVolume/CellCenterFVMData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a initalizing solution command
   *
   * @author Tiago Quintino
   * @author Andrea Lani
   *
   */
class InitState : public CellCenterFVMCom {
public:

  typedef Framework::VarSetTransformer VectorTransformer;

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor.
   */
  explicit InitState(const std::string& name);

  /**
   * Destructor.
   */
  virtual ~InitState();

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> >
  needsSockets();

  /**
   * Set up private data
   */
  virtual void setup();

  /**
   * Unsetup private data
   */
  virtual void unsetup();

protected:

  /**
   * Execute Processing actions
   */
  virtual void executeOnTrs();

  /**
   * Configures the command.
   */
  virtual void configure ( Config::ConfigArgs& args );

  /**
   * Read the file defining the vectorial functions
   */
  virtual void readFunctionsFile();
  
protected: // data
  
  /// handle to the face normals
  Framework::DataSocketSink< CFreal> socket_normals;
  
  /// storage of the states
  Framework::DataSocketSink < Framework::State* , Framework::GLOBAL > socket_states;
  
  /// storage of the ghost states
  Framework::DataSocketSink<Framework::State*> socket_gstates;
  
  /// storage of the nodes
  Framework::DataSocketSink < Framework::Node* , Framework::GLOBAL > socket_nodes;
  
  /// physical model var set
  Common::SafePtr<Framework::ConvectiveVarSet> _varSet;
  
  /// Transformer from input variables to update variables
  Common::SelfRegistPtr<VectorTransformer> _inputToUpdateVar;
  
  /// RealVector holding the input variables
  Framework::State* _input;

  /// Flag to input directly adimensional values
  bool _inputAdimensionalValues;

  /// a vector of string to hold the functions
  std::vector<std::string> _functions;

  /// a vector of string to hold the functions
  std::vector<std::string> _vars;

  /// name of the file where strings holding the functions are defined
  std::string _functionsFileName;
  
  // the VectorialFunction to use
  Framework::VectorialFunction _vFunction;

  /// a string to hold the name of the input variables
  std::string _inputVarStr;

}; // class InitState

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_InitState_hh

