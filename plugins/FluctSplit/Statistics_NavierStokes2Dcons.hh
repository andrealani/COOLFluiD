#ifndef COOLFluiD_FluctSplit_Statistics_NavierStokes2Dcons_hh
#define COOLFluiD_FluctSplit_Statistics_NavierStokes2Dcons_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/SubSystemStatus.hh"
#include "Framework/DataProcessingData.hh"

#include "NavierStokes/Euler2DCons.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {


    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/// This class samples the solution at a point for every iteration
/// @author Tiago Quintino
/// @author Erik Torres
class Statistics_NavierStokes2Dcons : public Framework::DataProcessingCom {
public: // functions

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor.
   */
  Statistics_NavierStokes2Dcons(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~Statistics_NavierStokes2Dcons();

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


private: // data

  /// handle the file where sampling is dumped
  /// handle is opened in setup() closed in unsetup()
  Common::SelfRegistPtr<Environment::FileHandlerOutput> m_fhandle;

  /// Index of the cell that contains the sampling point
  std::vector<CFuint> m_SamplingCellIndex;

  /// Pointer to TRS
  std::vector< Common::SafePtr< Framework::TopologicalRegionSet> > m_PointerToSampledTRS;

  /// the socket to the data handle of the state's
  Framework::DataSocketSink < Framework::Node* , Framework::GLOBAL > socket_nodes;

  /// the socket to the data handle of the state's
  Framework::DataSocketSink < Framework::State* , Framework::GLOBAL > socket_states;

  /// the socket stores the data of the mean flow
  Framework::DataSocketSource<RealVector> socket_meanflowNS;
  
  /// File name
  std::string m_OutputFile;

  /// Save Rate
  CFuint m_saveRate;


  /// My Path name
  boost::filesystem::path m_path_name;


  /// Average varuables that we want to compute

  std::vector< CFreal> m_rhobar ;
  std::vector< CFreal> m_Ubar;
  std::vector< CFreal> m_Vbar;
  std::vector< CFreal> m_UVbar;
  std::vector< CFreal> m_UUbar;
  std::vector< CFreal> m_VVbar;
  std::vector< CFreal> m_pbar ;
  std::vector< CFreal> m_rhoEbar ;
  std::vector< CFreal> m_rms ;
  std::vector< CFreal> m_devRho ;
  CFuint m_nbStep;
  
  std::vector< CFreal> mean_rhobar ;
  std::vector< CFreal> mean_Ubar;
  std::vector< CFreal> mean_Vbar;
  std::vector< CFreal> mean_UVbar;
  std::vector< CFreal> mean_UUbar;
  std::vector< CFreal> mean_VVbar;
  std::vector< CFreal> mean_pbar ;
  std::vector< CFreal> mean_rhoEbar ;
  std::vector< CFreal> mean_rms ;
  std::vector< CFreal> mean_devRho ;
  CFuint mean_nbStep;  
  
  /// physical model (in conservative variables)
  Common::SelfRegistPtr<Physics::NavierStokes::Euler2DCons> m_varSet;

  // restart from prev. statistics
  bool m_Restart;
  
  
}; // end of class SamplingPoint

//////////////////////////////////////////////////////////////////////////////

    }
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluctSplit_SamplingPoint_hh
