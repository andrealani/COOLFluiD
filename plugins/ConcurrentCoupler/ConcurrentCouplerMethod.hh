#ifndef COOLFluiD_Numerics_ConcurrentCoupler_ConcurrentCouplerMethod_hh
#define COOLFluiD_Numerics_ConcurrentCoupler_ConcurrentCouplerMethod_hh

//////////////////////////////////////////////////////////////////////////////

#include <memory>

#include "Framework/CouplerMethod.hh"
#include "ConcurrentCoupler/ConcurrentCouplerData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
    class NumericalCommand;
  }

  namespace Numerics {

    namespace ConcurrentCoupler {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a Coupler able to handle multiple subsystems (or methods 
 * operating in different namespaces within each subsystem) running concurrently
 * It couples TRSs from a Namespace of a subsystem with TRSs from another Namespace/Subsystem
 *
 * @author Andrea Lani
 *
 */
class ConcurrentCouplerMethod : public Framework::CouplerMethod {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Default constructor without arguments.
   */
  ConcurrentCouplerMethod(const std::string& name);

  /**
   * Default destructor.
   */
  ~ConcurrentCouplerMethod();

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
  
  /// Definition of an enumerator for I/O
  enum TypeIO {READ=0, WRITE=1};
  
  /// Tell if a coupling step has to be accomplished
  bool isCouplingIter(std::pair<std::ifstream*, std::ofstream*>& file, 
		      TypeIO tio);
    
  /// Reset given status file to 0 
  /// @param fout  file handle
  /// @param tio   type of I/O operation 
  void resetStatusFile(std::ofstream* fout, TypeIO tio) 
  {
    const CFuint start = tio*m_coupledNamespacesStr.size();
    const CFuint end = (tio+1)*m_coupledNamespacesStr.size();
    CFLog(INFO, "start/end = " << start << "/" << end << "\n");
    fout->seekp(start);
    for (CFuint f = start; f < end; ++f) {(*fout) << (bool)0;}
  }
  
  /// @return the status filename for this coupling namespace
  std::string getStatusFilename() const 
  {
    return std::string("status." + getNamespace());
  }
  
private: // member data
  
  ///The Setup command to use
  Common::SelfRegistPtr<ConcurrentCouplerCom> m_setup;
  
  ///The command to use for data transfer between the subsystems
  std::vector<Common::SelfRegistPtr<ConcurrentCouplerCom> > m_interfacesRead;
  std::vector<Common::SelfRegistPtr<ConcurrentCouplerCom> > m_interfacesWrite;
  
  /// file handles for I/O (reading and writing)
  std::pair<std::ifstream*, std::ofstream*> m_fileRW;
  
  /// flag telling if this is a root for I/O operations
  bool m_ioRoot;
    
  ///The data to share between ConcurrentCoupler commands
  Common::SharedPtr<ConcurrentCouplerData> m_data;
      
  ///The Setup string for configuration
  std::string m_setupStr;
  
  ///The coupled TRS command types and names
  std::vector<std::string> m_interfacesReadStr;
  std::vector<std::string> m_interfacesReadNameStr;
  
  std::vector<std::string> m_interfacesWriteStr;
  std::vector<std::string> m_interfacesWriteNameStr;
  
  ///The names of the subsytems to be coupled with
  std::vector<std::string> m_coupledSubSystemsStr;
  
  ///The names  of the subsytems to be coupled with
  std::vector<std::string> m_coupledNamespacesStr;
  
  ///Rate at which the data should be read/written for each of the interfaces
  std::vector<CFuint> m_transferRates;
  
}; // end of class ConcurrentCouplerMethod

//////////////////////////////////////////////////////////////////////////////

    } // namespace ConcurrentCoupler

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_ConcurrentCoupler_ConcurrentCouplerMethod_hh
