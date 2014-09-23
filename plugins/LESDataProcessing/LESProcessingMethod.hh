#ifndef COOLFluiD_Numerics_LESProcessing_LESProcessing_hh
#define COOLFluiD_Numerics_LESProcessing_LESProcessing_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/DataProcessingMethod.hh"
#include "LESProcessingData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
    class NumericalCommand;
  }

  namespace Numerics {

    namespace LESDataProcessing {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class defines a LESProcessingMethod that implements
 * the LESProcessinging of a solution
 *
 * @author Willem Deconinck
 *
 */
class Framework_API LESProcessing : public Framework::DataProcessingMethod {
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
  explicit LESProcessing(const std::string& name);

  /**
   * Default destructor
   */
  ~LESProcessing();

  /** 
   * Gets the Class name
   */
  static std::string getClassName()
  {
    return "LESProcessing";
  }

  /**
   * Configures the method, by allocating its dynamic members.
   *
   * @param args arguments from where to read the configuration
   */
  virtual void configure ( Config::ConfigArgs& args );

  /**
   * Gets the Data aggregator of this method
   * @return SafePtr to the MethodData
   */
  virtual Common::SafePtr<Framework::MethodData> getMethodData() const;

protected: // abstract interface implementations

  virtual void processDataImpl();

  /**
   * Sets up the data for the method commands to be applied.
   * @see Method::unsetMethod()
   */
  virtual void setMethodImpl();

  /**
   * UnSets the data of the method.
   * @see Method::setMethod()
   */
  virtual void unsetMethodImpl();
   
private: // helper functions

  // void clearPrepareComs();
  
private: // member data

  /// The Prepare command to use
  Common::SelfRegistPtr<LESProcessingComBase> m_prepare;
  
  /// The Prepare Type for the configuration
  std::string m_prepareTypeStr;
  
  /// The process commands
  std::vector<Common::SelfRegistPtr<LESProcessingComBase> > m_processes;

  ///The process Types for configuration
  std::vector<std::string> m_processTypeStr;

  ///The data to share between LESProcessingMethod commands
  Common::SharedPtr<LESProcessingData> m_data;

}; // class ExpLESProcessings

//////////////////////////////////////////////////////////////////////////////

    } // namespace LESProcessing

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_LESProcessing_LESProcessing_hh
