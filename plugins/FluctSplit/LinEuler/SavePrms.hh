#ifndef COOLFluiD_FluctSplit_SavePrms_hh
#define COOLFluiD_FluctSplit_SavePrms_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/SubSystemStatus.hh"
#include "Framework/DataProcessingData.hh"
#include "Environment/FileHandlerOutput.hh"
#include "Framework/PhysicalModel.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace Physics {

      namespace LinearizedEuler {

//////////////////////////////////////////////////////////////////////////////

class SavePrms : public Framework::DataProcessingCom {
public: // functions

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor.
   */
  SavePrms(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~SavePrms();

  /**
   * Sockets
   */

  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

  /**
   * Configure the command
   */
  virtual void configure ( Config::ConfigArgs& args );

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

  /// Execute this command on the TRS
  void execute();

private: // data

  /// the socket to the data handle of the state's
  Framework::DataSocketSink < Framework::State* , Framework::GLOBAL > socket_states;

  RealVector Prms;
  
  CFuint m_saveRate;
  
  
  /// handle the file where sampling is dumped
  /// handle is opened in setup() closed in unsetup()
  Common::SelfRegistPtr<Environment::FileHandlerOutput> m_fhandle;
  
  /// File name
  std::string m_OutputFile;  
  
  
  
  
}; 

//////////////////////////////////////////////////////////////////////////////

    } // namespace LinEuler

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif 
