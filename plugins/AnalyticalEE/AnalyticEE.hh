#ifndef COOLFluiD_Numerics_AnalyticalEE_AnalyticEE_hh
#define COOLFluiD_Numerics_AnalyticalEE_AnalyticEE_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/ErrorEstimatorMethod.hh"
#include "AnalyticEEData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework { class NumericalCommand; }

    namespace AnalyticalEE {

//////////////////////////////////////////////////////////////////////////////

/// This is a method which estimates the error of the computed solution
/// by compring it to a user defined analytical solution.
/// @author Tiago Quintino
class AnalyticEE : public Framework::ErrorEstimatorMethod
{
public: // methods

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Default constructor without arguments
  /// @param name missing documentation
  explicit AnalyticEE(const std::string& name);

  /// Default destructor
  ~AnalyticEE();

  /// Configures the method, by allocating the it's dynamic members.
  /// @param args arguments from where to read the configuration
  virtual void configure ( Config::ConfigArgs& args );

  /// Gets the Data aggregator of this method
  /// @return SafePtr to the MethodData
  virtual Common::SafePtr< Framework::MethodData > getMethodData () const;

protected: // abstract interface implementations

  /// Estimate the Error
  /// @see ErrorEstimatorMethod::estimate()
  virtual void estimateImpl();

  /// Sets up the data for the method commands to be applied.
  /// @see Method::unsetMethod()
  virtual void unsetMethodImpl();

  /// UnSets the data of the method.
  /// @see Method::setMethod()
  virtual void setMethodImpl();

private: // member data

  /// The Setup command to use
  Common::SelfRegistPtr<AnalyticEECom> m_setup;

  /// The UnSetup command to use
  Common::SelfRegistPtr<AnalyticEECom> m_unSetup;

  /// The command to compute the error
  Common::SelfRegistPtr<AnalyticEECom> m_compute;

  /// The string for configuration of Setup command
  std::string m_setupStr;

  /// The string for configuration of Unsetup command
  std::string m_unSetupStr;

  /// The string for configuration of Compute command
  std::string m_computeStr;

  /// The data to share between AnalyticalEEMethod commands
  Common::SharedPtr<AnalyticEEData> m_data;

}; // class AnalyticEE

//////////////////////////////////////////////////////////////////////////////

    } // namespace AnalyticalEE

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_AnalyticalEE_hh
