#ifndef COOLFluiD_DiscontGalerkin_InitState_hh
#define COOLFluiD_DiscontGalerkin_InitState_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/VectorialFunction.hh"
#include "Framework/VarSetTransformer.hh"
#include "Framework/DataSocketSink.hh"
#include "DiscontGalerkin/DiscontGalerkinSolverData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace DiscontGalerkin {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a initalizing solution command
 *
 * @author Tiago Quintino
 * @author Martin Holik
 */
class InitState : public DiscontGalerkinSolverCom {

public: // functions

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  typedef Framework::VarSetTransformer VectorTransformer;

  /**
   * Constructor.
   */
  explicit InitState(const std::string& name);

  /**
   * Destructor.
   */
  ~InitState();

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();


  /**
   * Set up the member data
   */
  virtual void setup();

  /**
   * UnSet up private data and data of the aggregated classes
   * in this command after processing phase
   */
  virtual void unsetup();

protected: // functions

  /**
   * Execute Processing actions
   */
  void executeOnTrs();

  /**
   * Configures the command.
   */
  virtual void configure ( Config::ConfigArgs& args );

protected: // data

  /// the socket to the data handle of the state's
  Framework::DataSocketSink < Framework::State* , Framework::GLOBAL > socket_states;

  /// a vector of string to hold the functions
  std::vector<std::string> _functions;

  /// a vector of string to hold the functions
  std::vector<std::string> _vars;

  /// the VectorialFunction to use
  Framework::VectorialFunction _vFunction;

}; // class InitState

//////////////////////////////////////////////////////////////////////////////

    } // namespace DiscontGalerkin

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_DiscontGalerkin_InitState_hh

