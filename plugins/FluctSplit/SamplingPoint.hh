#ifndef COOLFluiD_FluctSplit_SamplingPoint_hh
#define COOLFluiD_FluctSplit_SamplingPoint_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/SubSystemStatus.hh"
#include "Framework/DataProcessingData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {


    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/// This class samples the solution at a point for every iteration
/// @author Tiago Quintino
/// @author Erik Torres
class SamplingPoint : public Framework::DataProcessingCom {
public: // functions

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor.
   */
  SamplingPoint(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~SamplingPoint();

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

  /// TakeSample
  void TakeSample();

private: // data

  /// handle the file where sampling is dumped
  /// handle is opened in setup() closed in unsetup()
  Common::SelfRegistPtr<Environment::FileHandlerOutput> m_fhandle;

  /// Index of the cell that contains the sampling point
  CFuint m_SamplingCellIndex;

  /// Pointer to TRS
  Common::SafePtr< Framework::TopologicalRegionSet> m_PointerToSampledTRS;

  /// the socket to the data handle of the state's
  Framework::DataSocketSink < Framework::Node* , Framework::GLOBAL > socket_nodes;

  /// the socket to the data handle of the state's
  Framework::DataSocketSink < Framework::State* , Framework::GLOBAL > socket_states;

  /// File name
  std::string m_OutputFile;

  /// Probe location
  std::vector<CFreal> m_location;

  /// Save Rate
  CFuint m_SaveRate;

  /// Builder for standard TRS GeometricEntities
  Framework::GeometricEntityPool< Framework::StdTrsGeoBuilder >  m_cellbuilder;

  /// Stores whether the Element to probe has been found
  bool m_success;

  /// My Path name
  boost::filesystem::path m_path_name;


}; // end of class SamplingPoint

//////////////////////////////////////////////////////////////////////////////

    }
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluctSplit_SamplingPoint_hh
