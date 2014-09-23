#ifndef COOLFluiD_Numerics_AnalyticalEE_ComputeDiscreteError_hh
#define COOLFluiD_Numerics_AnalyticalEE_ComputeDiscreteError_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/SingleBehaviorFactory.hh"
#include "Environment/FileHandlerOutput.hh"
#include "Framework/DataSocketSink.hh"
#include "Framework/ProxyDofIterator.hh"

// #include "AnalyticalEE/AnalyticEEData.hh"
#include "AnalyticalEE/ComputeIntegralError.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace AnalyticalEE {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a MethodCommand that computes the solution error
/// based on a user-defined analytic function.
/// @author Tiago Quintino
class IntegralErrorFunction : public ComputeIntegralError {
public: // methods

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Constructor.
  explicit IntegralErrorFunction(std::string name);

  /// Destructor.
  ~IntegralErrorFunction();

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

  /// Compute the integral L2 norm in element
  /// @param cell is the celll in which we integrate
  /// @param error stores the computed values of error
  void integralErrorInCell(Framework::GeometricEntity *const cell, RealVector& error);

}; // class Setup

//////////////////////////////////////////////////////////////////////////////

    } // namespace AnalyticalEE

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_AnalyticalEE_ComputeDiscreteError_hh

