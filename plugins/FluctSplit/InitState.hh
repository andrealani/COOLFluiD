#ifndef COOLFluiD_Numerics_FluctSplit_InitState_hh
#define COOLFluiD_Numerics_FluctSplit_InitState_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/VectorialFunction.hh"
#include "Framework/VarSetTransformer.hh"
#include "FluctuationSplitData.hh"
#include "Framework/DataSocketSink.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a initalizing solution command
/// It is possible to introduce functions of x,y,z (one per state variable) 
/// or to input a first set of "named" functions as additional variables 
/// (e.g. radius, area, ... whatever you like) that will then be used later 
/// in the final expressions. The latter become functions of x,y,z, radius, area, 
/// and all other user-defined variables ... 
/// 
/// @author Tiago Quintino
/// @author Andrea Lani
class FluctSplit_API InitState : public FluctuationSplitCom {
public:

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  typedef Framework::VarSetTransformer VectorTransformer;

  /// Constructor.
  explicit InitState(const std::string& name);

  /// Destructor.
  ~InitState();

  /// Returns the DataSocket's that this command needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();


  /// Set up the member data
  virtual void setup();

  /// UnSet up private data and data of the aggregated classes
  /// in this command after processing phase
  virtual void unsetup();

protected:

  /// Execute Processing actions
  void executeOnTrs();

  /// Configures the command.
  virtual void configure ( Config::ConfigArgs& args );

protected: // data

  /// the socket to the data handle of the state's
  Framework::DataSocketSink < Framework::State* , Framework::GLOBAL > socket_states;

  /// physical model var set
  Common::SafePtr<Framework::ConvectiveVarSet> _varSet;

  /// Transformer from Update to Linear Variables
  Common::SelfRegistPtr<VectorTransformer> _inputToUpdateVar;

  /// RealVector holding the input variables
  Framework::State* _input;
  
  /// array for storing temporary function results 
  RealVector _tmpFun;
  
  /// array for storing temporary variables to pass to the 
  /// function parser 
  RealVector _tmpVars;
  
  /// the VectorialFunction to use
  Framework::VectorialFunction _vFunction;

  // the VectorialFunction to use
  Framework::VectorialFunction _vInitFunction;
  
  /// Flag to input directly adimensional values
  bool _inputAdimensionalValues;

  /// a vector of string to hold the functions
  std::vector<std::string> _functions;

  /// a vector of string to hold the functions
  std::vector<std::string> _vars;
 
  /// a vector of string to hold the functions
  std::vector<std::string> _initFunctions;
  
  /// a vector of string to hold the functions
  std::vector<std::string> _initVars;
   
  /// a string to hold the name of the input variables
  std::string _inputVarStr;
  
  /// IDs of the variables that will be changed interactively
  std::vector<CFuint> _interVarIDs;
    
  /// Factor to multiply the selected InteractiveVarIDs (should be < 1)
  CFreal _interFactor;
  
  /// Variable for initial Blasius profile
  /// Switch to decide whether Blasius inflow profile should be used 
  bool m_useBlasius;

  /// Reference length used to make all lengths non-dimensional
  CFreal m_ReferenceLength;
  
  ///(dimensional) position upstream of the leading edge
  CFreal m_UpstreamPosition;
  
  
  /// Reynolds number based on free stream quantities
  CFreal m_ReynoldsNumber;
  
   /// Mach number based on free stream quantities
  CFreal m_MachNumber;
  
   /// adiabatic coefficient gamma=c_p/c_v 
  CFreal m_gamma;
  
  ///Blasius inflow velocity in u-direction
  CFreal blasius_inflow;
  
  ///Blasius inflow energy flux
  CFreal blasius_energy_inflow;
  
}; // class InitState

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_InitState_hh

