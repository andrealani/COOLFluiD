#ifndef COOLFluiD_Numerics_Finite_VolumeICP_LorentzForceSourceTermComm_hh
#define COOLFluiD_Numerics_Finite_VolumeICP_LorentzForceSourceTermComm_hh

///////////////////////////////////////////////////////////////////////////

#include "Framework/DataProcessingData.hh"
#include "Framework/DataSocketSink.hh"
#include "Framework/CellTrsGeoBuilder.hh"
#include "Framework/GeometricEntityPool.hh"

//////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {  class PhysicalChemicalLibrary;  }

  namespace Numerics {

    namespace FiniteVolumeICP 
    {
      class VectorPotential;
    }
    
    namespace FiniteVolumeICP {


//////////////////////////////////////////////////////////////////////////

/**
 * This class implements the Lorentz force for the induction equations
 *
 * @author Vladimir Prokop
 * @author Emanuele Sartori
 *
 */
template <typename ST>
class LorentzForceSourceTerm : public Framework::DataProcessingCom {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor
   */
  LorentzForceSourceTerm(const std::string& name);

  /**
   * Default destructor
   */
  ~LorentzForceSourceTerm();

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
   * 
   * Needed for the new diamond shape computational stencil.
   * We use distance based values there. This function computes
   * distance weighted mean between 2 or 4 points.
   * 
   */
  CFreal distanceWeightedMean(CFreal value1, CFreal d1, CFreal value2, CFreal d2, CFreal value3 = 0., CFreal d3 = 0., CFreal value4 = 0., CFreal d4 = 0.);

private: //data

  /// array of temporary values
  std::vector<RealVector*> m_values;

  /// socket for the vacuum electric field intensity (real and imaginary component)
  Framework::DataSocketSink<RealVector> socket_vacuumElFieldIntensity;

  /// socket for the vacuum electric field intensity in nodes (real and imaginary component)
  Framework::DataSocketSink<RealVector> socket_vacuumElFieldIntensityInNodes;

  /// socket for the vacuum electric field intensity in ghost cells (real and imaginary component)
  Framework::DataSocketSink<RealVector> socket_vacuumElFieldIntensityInGhostCells;

  /// socket for the vacuum electric field intensity in ghost cells (real and imaginary component)
  Framework::DataSocketSink<CFreal> socket_elCondField;

  /// socket for the Lorentz Force (x and r component)
  Framework::DataSocketSource<RealVector> socket_LorentzForce;

  /// storage of volumes
  Framework::DataSocketSink<CFreal> socket_volumes;

  /// storage of nodes
  Framework::DataSocketSink < Framework::Node* , Framework::GLOBAL > socket_nodes;

  /// storage of states
  Framework::DataSocketSink < Framework::State* , Framework::GLOBAL > socket_states;

  /// storage of nstates (states in nodes)
  Framework::DataSocketSink<RealVector> socket_nstates;

  /// storage of normals
  Framework::DataSocketSink<CFreal> socket_normals;

  /// storage of isOutward
  Framework::DataSocketSink<CFint> socket_isOutward;

  /// socket for State's
  Framework::DataSocketSink<Framework::State*> socket_gstates;

  /// pointer to the physical-chemical library
  Common::SafePtr<Framework::PhysicalChemicalLibrary> m_library;

  /// builder of geometric entities
  Framework::GeometricEntityPool<Framework::StdTrsGeoBuilder> m_geoBuilder;

  /// builder of cells
  Common::SafePtr<Framework::GeometricEntityPool
		 <Framework::CellTrsGeoBuilder> > m_cellBuilder;
  
  /// source set
  Common::SafePtr<ST> m_srcTerm;
  
  /// array of temporary nodal states
  std::vector<RealVector*> m_states;

  /// Save each iteration to a different Tecplot file (with suffix _iter#).
  bool  m_appendIter;

  /// Save each timestep to a different Tecplot file (with suffix _iter#).
  bool  m_appendTime;

  /// name of the Output File where to write the Lorentz Force
  std::string m_nameOutputFileLorentzForce;

  // LorentzForce.FaceCenterComputationMethod:
  //  0: E in the middle of the face obtained from E in adjacent nodes 
  //  1: E in the middle of the face obtained with complete distance-base diamond-shape stencil (DEFAULT)
  //  2: E in the middle of the face obtained from E in cell centers
  CFuint m_FaceCenterComputationMethod;

  // LorentzForce.AverageInNodeStrategy:
  //  0: LorentzForce takes computed Ev coming from RMSJouleHeatSourceComm.cxx (DEFAULT)
  //  1: LorentzForce compute Ev usign distance-based average 
  //  2: LorentzForce computes Ev in the original way, usign volume-based average
  CFuint m_averageInNodeStrategy;

}; // end of class LorentzForceSourceTerm

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolumeICP

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolumeICP/LorentzForceSourceTermComm.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolumeICP_LorentzForceSourceTermComm_hh



