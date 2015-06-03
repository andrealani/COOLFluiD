#ifndef COOLFluiD_Physics_NavierStokes_RMSTurb_hh
#define COOLFluiD_Physics_NavierStokes_RMSTurb_hh

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

      
/// This class samples the solution at a point and computesb turbulent statistical data like fluctuations.
/// It makes sense to use for unsteady simulations
/// @ATTENTION for the moment it can be used only for fixed time step (Conformal Meshes).
/// @ATTENTION the process rate and the stopIter (that is the iteration at which it will stop) for the moment
/// will be defined by the DataProcessingMethod as
      //Simulator.SubSystem.DataProcessing2.ProcessRate = 
      //Simulator.SubSystem.DataProcessing2.StopIter = 
/// @author Gkoudesnes Christos
      
class RMSTurb : public Framework::DataProcessingCom
{
public: // functions

  /**
   * Defines the Config Option's of this class
   * @param options an OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor.
   */
  RMSTurb(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~RMSTurb();

  /**
   * Configure the command
   */
  virtual void configure ( Config::ConfigArgs& args );

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
 virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();
 
   /**
      Returns the DataSocket's that this command provides as sources
      @return a vector of SafePtr with the DataSockets
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

  /// Execute RMSTurb
  void execute();


  /**
   * Compute the RMS and turbulent stastistics values
   */
  void computeRMS();

  /**
   * Open the Output File and Write the header
   */
  //void prepareOutputFileRMS();

private: // data
  
  /// the socket to the data handle of the state's
  Framework::DataSocketSink < Framework::State* , Framework::GLOBAL > socket_states;
  
  /// socket for stencil
  //Framework::DataSocketSink<CFreal> socket_rmsturb;
  Framework::DataSocketSource<CFreal> socket_rmsturb;
  
  /// the dynamic sockets where the initial solutions will be read from the CFmesh file
  Framework::DynamicDataSocketSet<> m_dynamicSockets;
      
  /// physical data array 
  RealVector m_physicalData;
  
  // vector to save the additive rms results
  RealVector m_rmsresult;

  // number of rms iterations
  CFuint m_nbStep;
  
  bool m_restart;

  //flag to say if turbulent intensity and kinetic energy will also be computed
  bool m_comp;
  
  //flag to say if the net values will be also be computed (like net velocity fluctuation)
  bool m_net;
  
  CFuint m_InitSteps;
  
  // it is an assistant parameter
  CFuint resetIter; 
  
  // number of the variables to be saved
  CFuint nbOpts;
  
  // number of the auxiliary variables
  CFuint nbrms;
  
  // name of the initial solution for the RMS (Used when we restart the simulation)
  std::string m_rmsInitSocketName;
  
  ///@Note By default it will save the results in the CFmesh files and in the .plt with the other states values
  /// If it needs to be stored in other file uncomment the below comments
  
  /// handle the file where sampling is dumped
  /// handle is opened in setup() closed in unsetup()
  //Common::SelfRegistPtr<Environment::FileHandlerOutput> m_fhandle;
  
  // ///flag for appending iteration
  //bool m_appendIter;

  // ///flag for appending time
  //bool m_appendTime;
  
  /// Output File for RMSTurb values
  //Common::SelfRegistPtr<Environment::FileHandlerOutput> m_fileRMS;
  
  /// Storage for choosing when to save the wall values file
  //CFuint m_saveRateRMS;
  
  /// Name of Output File where to write the coeficients.
  //std::string m_nameOutputFileRMS;

}; // end of class RMSTurb

//////////////////////////////////////////////////////////////////////////////

    } // namespace NavierStokes

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_NavierStokes_RMS_hh
