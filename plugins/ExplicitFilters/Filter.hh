#ifndef COOLFluiD_Numerics_ExplicitFilters_Filter_hh
#define COOLFluiD_Numerics_ExplicitFilters_Filter_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/DataProcessingMethod.hh"
#include "FilterData.hh"
// #include "FilterCom.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
    class NumericalCommand;
  }

  namespace Numerics {

    namespace ExplicitFilters {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class defines a FilterMethod that implements
 * the filtering of a solution
 *
 * @author Willem Deconinck
 *
 */
class Framework_API Filter : public Framework::DataProcessingMethod {
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
  explicit Filter(const std::string& name);

  /**
   * Default destructor
   */
  ~Filter();

  /** 
   * Gets the Class name
   */
  static std::string getClassName()
  {
    return "Filter";
  }

  /**
   * Configures the method, by allocating its dynamic members.
   *
   * @param args arguments from where to read the configuration
   */
  virtual void configure ( Config::ConfigArgs& args );

  // /**
  //  * Gets a vector with all the NumericalStrategy's this method will use.
  //  * @return vector with the strategy pointers.
  //  */
  //  virtual std::vector<Common::SafePtr<Framework::NumericalStrategy> > getStrategyList () const;

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
  Common::SelfRegistPtr<FilterCom> m_prepare;
  
  /// The Prepare Type for the configuration
  std::string m_prepareTypeStr;
  
  /// The process commands
  std::vector<Common::SelfRegistPtr<FilterCom> > m_processes;

  ///The process Types for configuration
  std::vector<std::string> m_processTypeStr;

  ///The process Names for configuration
  std::vector<std::string> m_processNameStr;

  ///The data to share between ExplicitFiltersMethod commands
  Common::SharedPtr<FilterData> m_data;

}; // class ExpFilters

//////////////////////////////////////////////////////////////////////////////

    } // namespace ExplicitFilters

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_ExplicitFilters_Filter_hh
