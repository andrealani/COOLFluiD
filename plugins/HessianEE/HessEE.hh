#ifndef COOLFluiD_Numerics_HessianEE_HessEE_hh
#define COOLFluiD_Numerics_HessianEE_HessEE_hh

//////////////////////////////////////////////////////////////////////////////



#include "Framework/ErrorEstimatorMethod.hh"
#include "HessEEData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
    class NumericalCommand;
  }

  namespace Numerics {

    namespace HessianEE {

//////////////////////////////////////////////////////////////////////////////

/**
 *
 *
 * @author Tiago Quintino
 * @author Jurek Majevski
 *
 */
class HessEE : public Framework::ErrorEstimatorMethod
{
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Default constructor without arguments
   *
   * @param name missing documentation
   */
  explicit HessEE(const std::string& name);

  /**
   * Default destructor
   */
  ~HessEE();

  /**
   * Configures the method, by allocating the it's dynamic members.
   *
   * @param args arguments from where to read the configuration
   */
  virtual void configure ( Config::ConfigArgs& args );

protected: // abstract interface implementations

  /**
   * Gets the Data aggregator of this method
   * @return SafePtr to the MethodData
   */
  virtual Common::SafePtr< Framework::MethodData > getMethodData () const;

  /**
   * Estimate the Error
   * @see ErrorEstimatorMethod::estimate()
   */
  virtual void estimateImpl();

  /**
   * Sets up the data for the method commands to be applied.
   * @see Method::unsetMethod()
   */
  virtual void unsetMethodImpl();

  /**
   * UnSets the data of the method.
   * @see Method::setMethod()
   */
  virtual void setMethodImpl();

private: // member data

  ///The Setup command to use
  Common::SelfRegistPtr<HessEECom> _setup;

  ///The UnSetup command to use
  Common::SelfRegistPtr<HessEECom> _unSetup;

  ///The computeFunction command to use
  Common::SelfRegistPtr<HessEECom> _computeFunction;

  ///The computeHessian command to use
  Common::SelfRegistPtr<HessEECom> _computeHessian;

  ///The computeMetric command to use
  Common::SelfRegistPtr<HessEECom> _computeMetric;

  ///The smoother command to use
  Common::SelfRegistPtr<HessEECom> _smoother;

  ///The command used for updating global metric field
  Common::SelfRegistPtr<HessEECom> _updater;


  ///The Setup string for configuration
  std::string _setupStr;

  ///The UnSetup string for configuration
  std::string _unSetupStr;

  ///The updateSolution for configuration
  std::string _computeFunctionStr;

  ///The updateSolution for configuration
  std::string _computeHessianStr;

  ///The updateSolution for configuration
  std::string _computeMetricStr;

  ///The updateSolution for configuration
  std::string _smootherStr;

  ///The updateSolution for configuration
  std::string _updaterStr;

  ///The data to share between HessianEEMethod commands
  Common::SharedPtr<HessEEData> _data;

}; // class HessEE

//////////////////////////////////////////////////////////////////////////////

    } // namespace HessianEE

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_HessianEE_hh
