#ifndef COOLFluiD_Numerics_FiniteVolumeICP_RMSJouleHeatSourceCoupling_hh
#define COOLFluiD_Numerics_FiniteVolumeICP_RMSJouleHeatSourceCoupling_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/SingleBehaviorFactory.hh"
#include "Environment/FileHandlerOutput.hh"
#include "Framework/DataProcessingData.hh"
#include "Framework/DataSocketSink.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  
  namespace Framework {
    class PhysicalChemicalLibrary;
  }
    
  namespace Numerics {

    namespace FiniteVolumeICP 
    {
      class VectorPotential;
    }

    namespace FiniteVolumeICP {

//////////////////////////////////////////////////////////////////////////////

/**
 *
 * This class computes the Joule heat source to be used when coupling with the Induction Equation
 *
 * @author Radek Honzatko
 * @author Emanuele Sartori
 *
 */
template <typename ST>
class RMSJouleHeatSourceCoupling : public Framework::DataProcessingCom {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor
   */
  RMSJouleHeatSourceCoupling(const std::string& name);

  /**
   * Default destructor
   */
  ~RMSJouleHeatSourceCoupling();

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
   * Outputs the quality of the cells to a file
   */
  void printElCurrentToFile();

  /**
   *  Ev in cells and node is computed in these functions, and
   *  stored in a socket.
   *  Since current intensity in coils is not changing
   *  (we scale everithing at the end!) 
   *  these functions will be called only once (from ::setup method)
   *  
   *  Of course if you need it, you can compute Ev any time you like
   *  just run the correct method!
   */
  void prepareEvInCells();
  void prepareEvInNodes();

  /**
   * el conductivity is obtained from chemical library.
   * we may save some time if we call the library
   * just one time per iteration, saving results in a socket.
   */
  void prepareElConductivity();

  /**
   * Complex number: from a string to real,imaginary components
   */
  void sTc(std::string const stringCurrent, CFreal& currentRe, CFreal& currentIm);

private: //data

  /// socket for the vacuum electric field intensity (real and imaginary component)
  Framework::DataSocketSource<RealVector> socket_vacuumElFieldIntensity;

  /// socket for the vacuum electric field intensity in nodes (real and imaginary component)
  Framework::DataSocketSource<RealVector> socket_vacuumElFieldIntensityInNodes;

  /// socket for the vacuum electric field intensity in ghost cells (real and imaginary component)
  Framework::DataSocketSource<RealVector> socket_vacuumElFieldIntensityInGhostCells;

  /// socket for the time-averaged Joule heat source storage
  Framework::DataSocketSource<CFreal> socket_rmsJouleHeatSource;
  
  /// socket for the elCondfield storage
  Framework::DataSocketSource<CFreal> socket_elCondField;

  /// storage of volumes
  Framework::DataSocketSink<CFreal> socket_volumes;

  /// storage of nodes
  Framework::DataSocketSink < Framework::Node* , Framework::GLOBAL > socket_nodes;

  /// storage of states
  Framework::DataSocketSink < Framework::State* , Framework::GLOBAL > socket_states;

  /// storage of isOutward
  Framework::DataSocketSink<CFint> socket_isOutward;

  Framework::DataSocketSource<RealVector> socket_currentInCells;

  /// socket for the pressure
  Framework::DataSocketSource<CFreal> socket_pressure;
  
  /// thermodynamic library
  Common::SafePtr<Framework::PhysicalChemicalLibrary> m_library;

  /// builder of geometric entities
  Framework::GeometricEntityPool<Framework::StdTrsGeoBuilder> m_geoBuilder;

  /// source set
  Common::SafePtr<ST> m_srcTerm;
  
  /// Euler physical data
  RealVector m_physicalData;
  
  /// coefficient rescaling the electric current
  CFreal m_rescaleElCurrent;

  /// electric current amplitude - implicitly set to 1.0
  //CFreal m_elCurrent;
  RealVector m_elCurrent;
  std::string m_electricCurrent;
  //const CFuint REAL = 0;
  //const CFuint IMAGINARY = 1;

  /// array storing the electrical conductivity in the nodes
  RealVector m_radiusOfTheCell;

  /// desired torch power [kW]
  CFreal m_desiredPowerkW;

  /// number of coils
  CFuint m_nbCoils;

  /// radius of each circular coil
  std::vector<CFreal> m_radiusCoils;

  /// position on z-axis of each coil
  std::vector<CFreal> m_zPositionCoils;

  /// Save each iteration to a different Tecplot file (with suffix _iter#).
  bool  m_appendIter;

  /// Save each timestep to a different Tecplot file (with suffix _iter#).
  bool  m_appendTime;

  /// name of the Output File where to write the electric current
  std::string m_nameOutputFileElCurrent;

  /// Save each timestep to a different Tecplot file (with suffix _iter#).
  std::string m_toRun;
  
  /// flag to skip the preparation phase
  bool m_skipPreparation;

}; // end of class RMSJouleHeatSourceCoupling

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolumeICP

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "RMSJouleHeatSourceCoupling.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolumeICP_RMSJouleHeatSourceCoupling_hh
