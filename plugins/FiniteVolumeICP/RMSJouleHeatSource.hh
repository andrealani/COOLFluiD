#ifndef COOLFluiD_Numerics_FiniteVolumeICP_RMSJouleHeatSource_hh
#define COOLFluiD_Numerics_FiniteVolumeICP_RMSJouleHeatSource_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/SingleBehaviorFactory.hh"
#include "Environment/FileHandlerOutput.hh"
#include "Framework/DataProcessingData.hh"
#include "Framework/DataSocketSink.hh"
#include "NavierStokes/EulerVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  
  namespace Framework {
    class PhysicalChemicalLibrary;
  }

  namespace Numerics {

    namespace FiniteVolumeICP {

//////////////////////////////////////////////////////////////////////////////

/**
 *
 * This class computes the Joule heat source
 *
 * @author Radek Honzatko
 *
 */
class RMSJouleHeatSource : public Framework::DataProcessingCom {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor
   */
  RMSJouleHeatSource(const std::string& name);

  /**
   * Default destructor
   */
  ~RMSJouleHeatSource();

  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  void setup();

  /**
   * Unset up private data and data of the aggregated classes
   * in this command
   */
  void unsetup();

  /**
   * Execute on a set of dofs
   */
  void execute();

  /**
   * Configures this object with supplied arguments.
   */
  virtual void configure ( Config::ConfigArgs& args );

  /**
   * Returns the DataSocket's that this command needs as sources
   * @return a vector of SafePtr with the DataSockets
   */
  std::vector<Common::SafePtr<Framework::BaseDataSocketSource> > providesSockets();

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

private: //function

  /**
   * This function returns values of the complete elliptic integral of the first kind
   * to the modulus k
   */
  CFreal ellipticIntegralFirstKind(CFreal const& k);

  /**
   * This function returns values of the complete elliptic integral of the second kind
   * to the modulus k
   */
  CFreal ellipticIntegralSecondKind(CFreal const& k);

  /**
   * This function returns values of the combination of the complete elliptic
   * integrals of the first and second kind to the modulus k
   */
  CFreal ellipticIntegralCombined(CFreal const& k);

  /**
   * Prepare the output file for writing
   */
  void prepareOutputFile(std::ofstream& outputFile);

  /**
   * Write the electromagnetic field in the output file
   */
  void writeOutputFile();

  /**
   * Construct the file path to which to write the RMS source term
   */
  boost::filesystem::path constructFilename();

  /**
   * Outputs the quality of the cells to a file
   */
  void printElCurrentToFile();

private: //data

  /// socket for the time-averaged Joule heat source storage
  Framework::DataSocketSource<CFreal> socket_rmsJouleHeatSource;

  /// storage of volumes
  Framework::DataSocketSink<CFreal> socket_volumes;

  /// storage of states
  Framework::DataSocketSink < Framework::Node* , Framework::GLOBAL > socket_nodes;

  /// pointer to the physical-chemical library
  Common::SafePtr<Framework::PhysicalChemicalLibrary> m_library;
  
  /// builder of geometric entities
  Framework::GeometricEntityPool<Framework::StdTrsGeoBuilder> m_geoBuilder;
  
  /// Euler physical data
  RealVector m_physicalData;

  /// coefficient rescaling the electric current
  CFreal m_rescaleElCurrent;

  /// electric current amplitude - implicitly set to 1.0
  CFreal m_elCurrent;

  /// array storing the real component of electric field in the nodes
  RealVector m_elFieldReInNodes;

  /// array storing the imaginary component of electric field in the nodes
  RealVector m_elFieldImInNodes;

  /// array storing the Joule heat source in the nodes
  RealVector m_rmsJouleHeatSourceInNodes;

  /// frequency frequency of the torch [MHz]
  CFreal m_freqMHz;

  /// desired torch power [kW]
  CFreal m_desiredPowerkW;

  /// permeability of the free space
  CFreal m_permeability;

  /// number of coils
  CFuint m_nbCoils;

  /// radius of each circular coil
  std::vector<CFreal> m_radiusCoils;

  /// position on z-axis of each coil
  std::vector<CFreal> m_zPositionCoils;

  /// save rate
  CFuint  m_saveRate;

  /// Save each iteration to a different Tecplot file (with suffix _iter#).
  bool  m_appendIter;

  /// Save each timestep to a different Tecplot file (with suffix _iter#).
  bool  m_appendTime;

  /// name of the Output File where to write the electric current
  std::string m_nameOutputFileElCurrent;

  /// name of the Output File where to write the electric current
  std::string m_nameOutputFileEMField;

}; // end of class RMSJouleHeatSource

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolumeICP

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolumeICP_RMSJouleHeatSource_hh
