#ifndef COOLFluiD_FluctSplit_SaveSourceData_hh
#define COOLFluiD_FluctSplit_SaveSourceData_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/SubSystemStatus.hh"
#include "Framework/DataProcessingData.hh"
#include "Environment/FileHandlerOutput.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/DataSocketSink.hh"
#include "FluctuationSplitData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  
    namespace Physics {
    namespace NavierStokes {
      class Euler2DVarSet;
      class Euler3DVarSet;     
    }
  }

    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

class SaveSourceData : public Framework::DataProcessingCom {
public: // functions

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor.
   */
  SaveSourceData(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~SaveSourceData();

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
  
    /**
   * Compute the values at the wall and write them to file
   */
  void computeReynoldsGrad();
  
    /**
   * Compute the mean value of the sources
   */  
  void writeMean();

    /**
   * Open the Output File and Write the header
   */
  void prepareOutputFile();

private: // data

  /// the socket to the data handle of the state's
  Framework::DataSocketSink < Framework::State* , Framework::GLOBAL > socket_states;
  
  /// the socket to the data handle of the state's
  Framework::DataSocketSink < Framework::Node* , Framework::GLOBAL > socket_nodes;
  
    /// The sockets for the update coefficients
  Framework::DataSocketSink<CFreal> socket_volumes;
  
    /// the socket to the data handle of the normals
  Framework::DataSocketSink<InwardNormalsData*> socket_normals;
  
  CFuint DIM;
  
  // ///flag for appending iteration
  bool m_appendIter;

  // ///flag for appending time
  bool m_appendTime;
  
   // ///flag for appending iteration
  CFuint m_saveratemean;

  // ///flag for appending time
  bool m_writecoordinates; 
  
   /// Name of Output File where to write the coeficients.
  std::string m_nameOutputFile;
  
  Common::SelfRegistPtr<Environment::FileHandlerOutput> m_file;

  std::vector<RealVector> sources;
  std::vector<RealVector> sources_mean;
  
  CFint meanCount;

}; 

//////////////////////////////////////////////////////////////////////////////

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif 
