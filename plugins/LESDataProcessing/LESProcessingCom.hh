#ifndef COOLFluiD_Numerics_LESDataProcessing_LESProcessingCom_hh
#define COOLFluiD_Numerics_LESDataProcessing_LESProcessingCom_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/MethodCommand.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/DynamicDataSocketSet.hh"
#include "Framework/BaseDataSocketSink.hh"
#include "Framework/BaseDataSocketSource.hh"

//////////////////////////////////////////////////////////////////////////////


namespace COOLFluiD {
  
  namespace Numerics {

    namespace LESDataProcessing {
      
      class LESProcessingData;
      
      typedef Framework::MethodCommand<LESProcessingData> LESProcessingComBase;

      typedef Framework::MethodCommand<LESProcessingData>::PROVIDER LESProcessingComProvider;

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
class LESProcessingCom : public LESProcessingComBase {
public:

  enum SocketType {SOURCE, SINK, UNDEFINED};

  /**
   * Constructor.
   */
  LESProcessingCom(const std::string& name) : LESProcessingComBase(name) {}

  /**
   * Destructor.
   */
  virtual ~LESProcessingCom()
  {
  }

  static void defineConfigOptions(Config::OptionList& options) 
  {
  }

  virtual void setup()
  {
    LESProcessingComBase::setup();
  }
  
  virtual void unsetup()
  {
    LESProcessingComBase::setup();
  }

  virtual void configure ( Config::ConfigArgs& args )
  {
    LESProcessingComBase::configure(args);
  }

  /**
   * Execute Processing actions
   */
  virtual void execute()
  {
  }

	/**
	 * needsSockets()
	 * @return	a vector of SafePtr with the DataSockets needed as sinks
	 */
  virtual std::vector< Common::SafePtr< Framework::BaseDataSocketSink > >
    needsSockets()
    {
      CFAUTOTRACE;
      std::vector< Common::SafePtr< Framework::BaseDataSocketSink > > result = m_sockets.getAllSinkSockets();
      std::vector< Common::SafePtr< Framework::BaseDataSocketSink > > globalSockets = m_globalSockets.getAllSinkSockets();
      for(CFuint i=0; i<globalSockets.size(); ++i) {
        result.push_back(globalSockets[i]);
      }
      return result;
    }
    
  virtual std::vector< Common::SafePtr< Framework::BaseDataSocketSource > >
    providesSockets()
    {
      std::vector< Common::SafePtr< Framework::BaseDataSocketSource > > result = m_sockets.getAllSourceSockets();
      std::vector< Common::SafePtr< Framework::BaseDataSocketSource > > globalSockets = m_globalSockets.getAllSourceSockets();
      for(CFuint i=0; i<globalSockets.size(); ++i) {
        result.push_back(globalSockets[i]);
      }
      return result;
    }

  
protected:
    
  /// Dynamic data sockets
  Framework::DynamicDataSocketSet<>                  m_sockets;
  Framework::DynamicDataSocketSet<Framework::GLOBAL> m_globalSockets;
  
}; // class Setup

//////////////////////////////////////////////////////////////////////////////

    } // namespace LESDataProcessing

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_LESDataProcessing_LESProcessingCom_hh

