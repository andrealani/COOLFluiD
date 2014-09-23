#ifndef COOLFluiD_Numerics_ParMetisBalancer_ParMetisBalancer_hh
#define COOLFluiD_Numerics_ParMetisBalancer_ParMetisBalancer_hh

//////////////////////////////////////////////////////////////////////////////



#include "Framework/DynamicBalancerMethod.hh"
#include "ParMetisBalancerData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace ParMetisBalancer {

//////////////////////////////////////////////////////////////////////////////

/**
 * This clas is suposed to pervorm load balancing with use of ParMetis functions
 *
 * @author
 *
 */

class ParMetisBalancer : public Framework::DynamicBalancerMethod{
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor.
   *
   * @param name missing documentation
   */
  explicit ParMetisBalancer(const std::string& name);


  /**
   * Default destructor.
   */
  ~ParMetisBalancer();

  /**
   * Configures the method, by allocating its dynamic members.
   *
   * @param args missing documentation
   */
  virtual void configure (Config::ConfigArgs& args);

protected: // interface implementation functions

  /**
   * Gets the Data aggregator of this method
   * @return SafePtr to the MethodData
   */
   virtual Common::SafePtr< Framework::MethodData > getMethodData () const;

  /**
   * Sets up the data, commands and strategies of this Method
   *
   * @see Method::setMethod()
   */
  virtual void setMethodImpl();

  /**
   * Unsets the data, commands and strategies of this Method
   *
   * @see Method::unsetMethod()
   */
  virtual void unsetMethodImpl();

   /**
   * Do Repartitioning
   */
  virtual void doDynamicBalanceImpl();

protected: // helper functions

private: // member data

  ///The Setup string for configuration
  std::string m_setupStr;

  ///The UnSetup string for configuration
  std::string m_unSetupStr;

  ///The Setup command to use
  Common::SelfRegistPtr<ParMetisBalancerCom> m_setup;

  ///The UnSetup command to use
  Common::SelfRegistPtr<ParMetisBalancerCom> m_unSetup;

  /// string for configuring the  command
  std::string m_repartStr;

  /// The command to use
  Common::SelfRegistPtr<ParMetisBalancerCom> m_repart;

 /// The data to share between CFmeshReader commands
  Common::SharedPtr<ParMetisBalancerData> m_data;

}; // class FluctuationSplit

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_FluctuationSplit_hh
