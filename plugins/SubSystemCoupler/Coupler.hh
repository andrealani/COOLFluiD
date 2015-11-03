#ifndef COOLFluiD_Numerics_SubSystemCoupler_Coupler_hh
#define COOLFluiD_Numerics_SubSystemCoupler_Coupler_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/CouplerMethod.hh"
#include "SubSysCouplerData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
    class NumericalCommand;
  }

  namespace Numerics {

    namespace SubSystemCoupler {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a Coupler.
 * It couples TRSs from a Namespace of a subsystem with TRSs from another Namespace/Subsystem
 *
 * @author Thomas Wuilbaut
 *
 */
class Coupler : public Framework::CouplerMethod {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Default constructor without arguments.
   */
  Coupler(const std::string& name);

  /**
   * Default destructor.
   */
  ~Coupler();

  /**
   * Configures the method, by allocating the it's dynamic members.
   *
   * @param args arguments from where to read the configuration
   */
  virtual void configure ( Config::ConfigArgs& args );

  /**
   * Read the info sent by other subsystems to be able to communicate
   */
  void getInfoFromOtherSubSystem();

  /**
   * Transfer the info needed by other subsystems to be able to communicate
   */
  void setInfoToOtherSubSystem();

  /**
   * Finish configuring this Method.
   */
  virtual void postConfigure( Config::ConfigArgs& args );

  /**
   * Gets a vector with all the NumericalStrategy's this method will use.
   * @return vector with the strategy pointers.
   */
  virtual std::vector<Common::SafePtr<Framework::NumericalStrategy> > getStrategyList () const;

protected: // abstract interface implementations

  /**
   * Gets the Data aggregator of this method
   * @return SafePtr to the MethodData
   */
  virtual Common::SafePtr< Framework::MethodData > getMethodData () const;

  /**
   * Writes the data necessary for the preProcessing
   */
  virtual void preProcessWriteImpl();

  /**
   * Reads the data necessary for the preProcessing
   */
  virtual void preProcessReadImpl();

  /**
   * Writes the data necessary for the mesh matching
   */
  virtual void meshMatchingWriteImpl();

  /**
   * Reads the data necessary for the mesh matching
   */
  virtual void meshMatchingReadImpl();

  /**
   * Reads the data from the Coupled SubSystems
   */
  virtual void dataTransferReadImpl();

  /**
   * Writes the data for the Coupled SubSystems
   */
  virtual void dataTransferWriteImpl();

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
  
  /// Finalize the coupling
  virtual void finalizeImpl() {}
  
private: // helper methods

  /**
   * Read the info sent by other subsystems to be able to communicate
   */
  void readInfoFromOtherSubSystem();

  /**
   * Configures the CommandGroups as Interfaces
   */
  void configureInterfaces();

  /**
   * Clear Commands
   */
  void clearComds();

  /**
   * Helper function that gets a line from a file and puts
   * into a string incrementing a supplied counter
   * @param fin file stream
   * @param line string to put the line into
   * @param lineNb line counter
   */
  void getWordsFromLine(std::ifstream& fin,
          std::string& line,
          CFuint&  lineNb,
          std::vector<std::string>& words);


private: // member data

  ///The Setup command to use
  Common::SelfRegistPtr<CouplerCom> m_setup;

  ///The PreProcess command to use
  std::vector<Common::SelfRegistPtr<CouplerCom> > m_preProcessRead;
  std::vector<Common::SelfRegistPtr<CouplerCom> > m_preProcessWrite;

  ///The coupler command to use
  std::vector<Common::SelfRegistPtr<CouplerCom> > m_matchMeshesRead;
  std::vector<Common::SelfRegistPtr<CouplerCom> > m_matchMeshesWrite;

  ///The command to use for data transfer between the subsystems
  std::vector<Common::SelfRegistPtr<CouplerCom> > m_interfacesRead;
  std::vector<Common::SelfRegistPtr<CouplerCom> > m_interfacesWrite;

  ///The PostProcess command to use
  std::vector<Common::SelfRegistPtr<CouplerCom> > m_postProcess;

  ///The Setup string for configuration
  std::string m_setupStr;

  ///The MeshMatcher string for configuration
  std::vector<std::string> m_matchMeshesReadStr;
  std::vector<std::string> m_matchMeshesReadNameStr;
  std::vector<std::string> m_matchMeshesWriteStr;
  std::vector<std::string> m_matchMeshesWriteNameStr;

  ///The preProcessing commands for configuration
  std::vector<std::string> m_preProcessReadStr;
  std::vector<std::string> m_preProcessReadNameStr;
  std::vector<std::string> m_preProcessWriteStr;
  std::vector<std::string> m_preProcessWriteNameStr;


  ///The postProcessing commands for configuration
  std::vector<std::string> m_postProcessStr;
  std::vector<std::string> m_postProcessNameStr;

  ///The coupled TRS command types and names
  std::vector<std::string> m_interfacesReadStr;
  std::vector<std::string> m_interfacesReadNameStr;

  std::vector<std::string> m_interfacesWriteStr;
  std::vector<std::string> m_interfacesWriteNameStr;

  ///The names of the interfaces to be coupled
  std::vector<std::string> m_interfaceNameStr;

  ///The names  of the subsytems to be coupled with
  std::vector<std::string> m_coupledSubSystemsStr;

  ///The names  of the subsytems to be coupled with
  std::vector<std::string> m_coupledNamespacesStr;

  ///The data to share between SubSystemCoupler commands
  Common::SharedPtr<SubSysCouplerData> m_data;

  //Used to check that the method is fully configured (config + setInfo + getInfo + postConfig)
  CFuint m_fullConfigure;

  ///Rate at which the data should be read/written for each of the interfaces
  std::vector<CFuint> m_transferRates;

}; // end of class Coupler

//////////////////////////////////////////////////////////////////////////////

    } // namespace SubSystemCoupler

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_SubSystemCoupler_Coupler_hh
