#ifndef COOLFluiD_Numerics_LESDataProcessing_TimeAveraging_hh
#define COOLFluiD_Numerics_LESDataProcessing_TimeAveraging_hh

//////////////////////////////////////////////////////////////////////////////

#include "LESProcessingData.hh"
#include "Framework/SubSystemStatus.hh"
#include "LESProcessingCom.hh"

//////////////////////////////////////////////////////////////////////////////


namespace COOLFluiD {
  
  namespace Numerics {

    namespace LESDataProcessing {

//////////////////////////////////////////////////////////////////////////////

/**
 * A Data post processing command that time-averages the solution, and 
 * velocity products. This allows calculation of the mean turbulent fluctuations
 * <u'u'> = <uu> - <u><u>
 * This is stored in sockets for restart or plotting purposes.
 * 
 * Write for Restart: 
 * Simulator.SubSystem.CFmesh.WriteSol = ParWriteSolution
 * Simulator.SubSystem.CFmesh.Data.ExtraStateVarNames = averageSolution averageVelocityProducts
 * Simulator.SubSystem.CFmesh.Data.ExtraStateVarStrides = 'nbEqs' 'nbVelProducts (1D:1, 2D:3, 3D:6)'
 * Simulator.SubSystem.CFmesh.Data.ExtraVarNames = turbulenceAveragingSteps
 * Simulator.SubSystem.CFmesh.Data.ExtraVarStrides = 1
 * 
 * Read restart --> make sure that FirstTimeCreation = false
 * Simulator.SubSystem.CFmeshFileReader.ReadCFmesh = ParReadCFmesh
 * Simulator.SubSystem.CFmeshFileReader.Data.ExtraStateVarNames = averageSolution averageVelocityProducts
 * Simulator.SubSystem.CFmeshFileReader.Data.ExtraStateVarStrides = 'nbEqs' 'nbVelProducts (1D:1, 2D:3, 3D:6)'
 * Simulator.SubSystem.CFmeshFileReader.Data.ExtraVarNames = turbulenceAveragingSteps
 * Simulator.SubSystem.CFmeshFileReader.Data.ExtraVarStrides = 1
 * 
 * If restart file doesn't contain the required data 
 * --> set FirstTimeCreation = true
 * 
 * @author Willem Deconinck
 *
 */
class TimeAveraging : public LESProcessingCom {
public:

  /**
   * Constructor.
   */
  TimeAveraging(const std::string& name);

  /**
   * Destructor.
   */
  virtual ~TimeAveraging()
  {
  }

  static void defineConfigOptions(Config::OptionList& options);

  virtual void setup();
  
  virtual void unsetup();

  virtual void configure ( Config::ConfigArgs& args );

  Framework::DataHandle<CFreal> getDataHandle(const std::string& socketName) ;
  
  void makeSocketAvailable(const std::string& socketName, SocketType socketType = UNDEFINED);

  /**
   * Execute Processing actions
   */
  virtual void execute();

  // /**
  //  * needsSockets()
  //  * @return a vector of SafePtr with the DataSockets needed as sinks
  //  */
  //   virtual std::vector< Common::SafePtr< Framework::BaseDataSocketSink > >
  //     needsSockets();
  //     
  //   virtual std::vector< Common::SafePtr< Framework::BaseDataSocketSource > >
  //     providesSockets();

  bool isSourceSocket(const std::string& socketName) { return !(m_socketMap.find(socketName)); }

protected:
  
  Common::CFMap<std::string,bool>                    m_socketMap;

  /// arrray of gradients
  std::vector<RealVector> m_gradients;

private:
  
  std::vector<Common::SafePtr <TurbulenceFunction> > m_turbulenceFunctions;
        
  /// flag that defines if averaging is performed
  bool m_averaging;
  
  // stride for velocityProducts
  CFuint m_stride;
  
  // Counter for number of steps that has been averaged
  CFuint m_avgStepsCounter;
  
  /// State to dimensionalize updateVar
  RealVector m_primState;
    
  /// Flag that defines if averaging starts from zero.
  bool m_resetFlag;
  
  /// Flag that defines if source sockets must be created instead of sink sockets 
  bool m_firstTimeCreation;
  
  /// Flag that tells if averaging happens from nodal values or state values
  bool m_nodal;
  
  /// number of states
  CFuint m_nbStates;
  
}; // class Setup

//////////////////////////////////////////////////////////////////////////////

    } // namespace LESDataProcessing

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_LESDataProcessing_TimeAveraging_hh

