/*
 * This program was created to receive user specified values (and coefficients = need implementation) which should 
 * be time averaged during an unsteady simulation using LES or DES.  
 */
#ifndef COOLFluiD_Physics_NavierStokes_RMS_hh
#define COOLFluiD_Physics_NavierStokes_RMS_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/FileHandlerOutput.hh"
#include "Common/Trio.hh"
#include "MathTools/FunctionParser.hh"
#include "Environment/SingleBehaviorFactory.hh"
#include "Framework/DataProcessingData.hh"
#include "Framework/DynamicDataSocketSet.hh"
#include "Framework/DataSocketSink.hh"
#include "Framework/DataSocketSource.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
    class ConvectiveVarSet;
  }
  
  namespace Physics {

    namespace NavierStokes {
      
      
//////////////////////////////////////////////////////////////////////////////

/// This class samples the solution at a point for every iteration
/// @author
      
class RMS : public Framework::DataProcessingCom {
public: // functions

  /**
   * Defines the Config Option's of this class
   * @param options an OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor.
   */
  RMS(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~RMS();

  /**
   * Configure the command
   */
  virtual void configure ( Config::ConfigArgs& args );

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
 virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > 
 needsSockets();

  /**
   * Returns the DataSocket's that this command provides as sources
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSource> > providesSockets();
  
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
  void computeRMS(bool issave );

    /**
   * Open the Output File and Write the header
   */
  void prepareOutputFileRMS();

private: // data

  /// handle the file where sampling is dumped
  /// handle is opened in setup() closed in unsetup()
  Common::SelfRegistPtr<Environment::FileHandlerOutput> m_fhandle;
  
  /// the socket to the data handle of the state's
  Framework::DataSocketSink < Framework::State* , Framework::GLOBAL > socket_states;
  
  /// socket for stencil
  Framework::DataSocketSource<CFreal> socket_rms;
  
  /// Output File for Wall Values
  Common::SelfRegistPtr<Environment::FileHandlerOutput> m_fileRMS;
  
  /// Storage for choosing when to save the wall values file
  CFuint m_saveRateRMS;
  CFuint m_compRateRMS;
  
  /// Name of Output File where to write the coeficients.
  std::string m_nameOutputFileRMS;
      
  /// physical data array 
  RealVector m_physicalData;
  
  /// user options
  std::vector<std::string> m_rmsOpt;
  
  // vector to save the rms results
  RealVector m_rmsresult;

  //RealVector m_AvgSkinFriction;
  //RealVector m_AvgHeatFlux;
  CFuint m_nbStep;
  
  // ///flag for appending iteration
  bool m_appendIter;

  // ///flag for appending time
  bool m_appendTime;

  bool m_restart;
  
  //We add these options to let the user define which values will be time averaged
  
  
  ///@todo for aerocoef
  //bool m_DoAvgHeatFlux; 
  //bool m_DoAvgSkinFric;


  CFuint m_InitSteps;

}; // end of class SamplingPoint

//////////////////////////////////////////////////////////////////////////////

    } // namespace NavierStokes

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_NavierStokes_RMS_hh
