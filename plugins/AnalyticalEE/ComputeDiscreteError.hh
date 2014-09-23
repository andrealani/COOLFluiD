#ifndef COOLFluiD_Numerics_AnalyticalEE_ComputeDiscreteError_hh
#define COOLFluiD_Numerics_AnalyticalEE_ComputeDiscreteError_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/SingleBehaviorFactory.hh"
#include "Environment/FileHandlerOutput.hh"
#include "Framework/DataSocketSink.hh"
#include "Framework/ProxyDofIterator.hh"

#include "AnalyticalEE/AnalyticEEData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace AnalyticalEE {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a MethodCommand that computes the solution error
/// based on a user-defined analytic function.
/// @author Tiago Quintino
class ComputeDiscreteError : public AnalyticEECom {
public: // methods

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Constructor.
  explicit ComputeDiscreteError(std::string name);

  /// Destructor.
  ~ComputeDiscreteError();

  /// Set up private data and data of the aggregated classes
  /// in this command before processing phase
  void setup();

  /// Unset up private data and data of the aggregated classes
  /// in this command
  void unsetup();

  /// Configure the command
  virtual void configure ( Config::ConfigArgs& args );

  /// Returns the DataSocket's that this command provides as sources
  /// @return a vector of SafePtr with the DataSockets
  std::vector<Common::SafePtr<Framework::BaseDataSocketSource> > providesSockets();

  /// Returns the DataSocket's that this command needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

  /// Execute Processing actions
  void execute();

protected: // methods

  /// Open the output file, write the header and close the file handle
  void prepareOutputFile();

  /// Open the output file
  /// @param bool indicates if reopen should be done with append
  void openFile(bool append);

protected: // data

  /// the socket to the data handle of the node's
  Framework::DataSocketSink<Framework::Node*, Framework::GLOBAL> socket_nodes;

  /// the socket to the data handle of the nStatesProxy
  Framework::DataSocketSink<Framework::State*, Framework::GLOBAL> socket_states;

  /// a vector of string to hold the functions
  std::vector<std::string> m_functions;

  /// a vector of string to hold the variables
  std::vector<std::string> m_vars;

  /// the VectorialFunction to use
  Framework::VectorialFunction m_vFunction;

  /// vector to store the temporary value of the analytical solution
  RealVector m_exact;

  /// vector to store the temporary value of the error
  RealVector m_error;

  /// vector to store the temporary value of the L2_norm
  RealVector m_L2;

/// vector to store the temporary value of the logarithm of the L2_norm
  RealVector m_result;

  /// Output file for recording the evolution of the Error norm
  Common::SelfRegistPtr<Environment::FileHandlerOutput> m_file;

  /// Name of output file
  std::string m_nameOutputFile;

}; // class Setup

//////////////////////////////////////////////////////////////////////////////

    } // namespace AnalyticalEE

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_AnalyticalEE_ComputeDiscreteError_hh

