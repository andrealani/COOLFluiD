#ifndef COOLFluiD_Physics_LinEuler_CreateMeanFlowAnalytic_hh
#define COOLFluiD_Physics_LinEuler_CreateMeanFlowAnalytic_hh

//////////////////////////////////////////////////////////////////////////////

#include "MathTools/FunctionParser.hh"
#include "Framework/DataProcessingData.hh"
#include "Framework/DataSocketSource.hh"
#include "Framework/DataSocketSink.hh"
#include "Framework/VectorialFunction.hh"

#include "Framework/PhysicalModel.hh"
//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

   namespace LinearizedEuler {

//////////////////////////////////////////////////////////////////////////////

///
/// This class create a socket for the mean flow variables
/// @author Lilla Koloszar


class CreateMeanFlowAnalytic : public Framework::DataProcessingCom {
public: // functions

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor.
   */
  CreateMeanFlowAnalytic(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~CreateMeanFlowAnalytic();

  /**
   * Configure the command
   */
  virtual void configure ( Config::ConfigArgs& args );

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();


// This has to be deleted *****************************************************
  /**
   * Returns the DataSocket's that this command provides as sources
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSource> > providesSockets();
// This has to be deleted *****************************************************


  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  virtual void setup();

  /**
   * Unset up private data and data of the aggregated classes
   * in this command
   */
  virtual void unsetup();

protected: // functions

  /**
   * Execute on a set of dofs
   */
  void executeOnTrs();

private: // data

  /// the socket to the data handle of the state's
  Framework::DataSocketSink < Framework::State* , Framework::GLOBAL > socket_states;

// This has to be changed *****************************************************

  /// the socket stores the data of the mean flow
  Framework::DataSocketSource<RealVector> socket_meanflow;

// This has to be changed *****************************************************

  /// the socket to the data handle of the meanflow
//   Framework::DataSocketSink < Framework::State* , Framework::GLOBAL > socket_meanflow;

// This has to be changed *****************************************************

  /// temporry for the variables
  RealVector m_var_values;

  /// a vector of string to hold the functions
  std::vector<std::string> m_function_meanflow;

  /// a vector of string to hold the variables
  std::vector<std::string> m_vars_meanflow;

  /// function for the mean flow
  Framework::VectorialFunction m_function_parser_meanflow;

}; // end of class CreateMeanFlowAnalytic

//////////////////////////////////////////////////////////////////////////////

    } // namespace LinEuler

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_LinEuler_CreateMeanFlowAnalytic_hh
