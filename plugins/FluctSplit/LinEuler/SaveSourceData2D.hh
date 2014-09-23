#ifndef COOLFluiD_FluctSplit_SaveSourceData2D_hh
#define COOLFluiD_FluctSplit_SaveSourceData2D_hh

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

class SaveSourceData2D : public Framework::DataProcessingCom {
public: // functions

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor.
   */
  SaveSourceData2D(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~SaveSourceData2D();

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

  /// the socket stores the data of the decomposed sources
  Framework::DataSocketSink<RealVector> socket_sourceSave;

}; 

//////////////////////////////////////////////////////////////////////////////

    } // namespace LinEuler

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif 
