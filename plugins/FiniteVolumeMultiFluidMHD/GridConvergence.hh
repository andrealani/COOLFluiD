#ifndef COOLFluiD_Numerics_FiniteVolume_GridConvergence_hh
#define COOLFluiD_Numerics_FiniteVolume_GridConvergence_hh

///////////////////////////////////////////////////////////////////////////

#include "Framework/DataProcessingData.hh"
#include "Framework/DataSocketSink.hh"
#include "Framework/DataSocketSource.hh"
#include "Framework/CellTrsGeoBuilder.hh"
#include "Framework/GeometricEntityPool.hh"
#include "Framework/PhysicalConsts.hh"

//////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  
  namespace Framework {
    class PhysicalChemicalLibrary;
  }

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////

/**
 * This class compute some properties
 *  of the Sun atmosphere
 *
 * @author Alejandro Alvarez
 *
 */
class GridConvergence : public Framework::DataProcessingCom {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor
   */
  GridConvergence(const std::string& name);

  /**
   * Default destructor
   */
  ~GridConvergence();

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
   * Compute the analytical solution for the case of the isodensity vortex
   */
  void computeAnalyticalSolution();

  /**
   * Method writing the solution
   */
  void writeOutputFile();

  /**
   * Method constructing the file
   */
  boost::filesystem::path constructFilename();

  /**
  * Prepare the header of the output file for writing
  */
 void prepareOutputFile(std::ofstream& outputFile);

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

  /// save rate
  CFuint  m_saveRate;

  /// Save each iteration to a different Tecplot file (with suffix _iter#).
  bool  m_appendIter;

  /// Save each timestep to a different Tecplot file (with suffix _iter#).
  bool  m_appendTime;

  /// name of the Output File where to write the electric current
  std::string m_nameOutputFileError;

private: //function

private: //data

  /// storage of states
  Framework::DataSocketSink < Framework::State* , Framework::GLOBAL > socket_states;
  
  /// storage of ghost states
  /// this state list has ownership on the states
  Framework::DataSocketSink<Framework::State*> socket_gstates;
  
  /// storage of nstates (states in nodes)
  Framework::DataSocketSink<RealVector> socket_nstates;  
  
  /// storage of volumes
  Framework::DataSocketSink<CFreal> socket_volumes;  
  
  /// storage of nodes
  Framework::DataSocketSink < Framework::Node* , Framework::GLOBAL > socket_nodes;  
  
  /// storage of isOutward
  Framework::DataSocketSink<CFint> socket_isOutward; 
  
  /// storage of normals
  Framework::DataSocketSink<CFreal> socket_normals;  
  
  /// Bx theory solution
  Framework::DataSocketSource <CFreal> socket_BxTheory;

  /// By theory solution
  Framework::DataSocketSource <CFreal> socket_ByTheory;

  /// Ez theory solution
  Framework::DataSocketSource <CFreal> socket_EzTheory;

  /// Rho theory solution
  Framework::DataSocketSource <CFreal> socket_RhoTheory;
  
  /// Ux theory solution
  Framework::DataSocketSource <CFreal> socket_UxTheory;
  
  /// Uy theory solution
  Framework::DataSocketSource <CFreal> socket_UyTheory;

  /// T theory solution
  Framework::DataSocketSource <CFreal> socket_TTheory;
    
  /// pointer to the physical-chemical library In case in the future we take
  /// the cross-sections from library
  //Common::SafePtr<Framework::PhysicalChemicalLibrary> m_library;
  
  /// builder of geometric entities
  Framework::GeometricEntityPool<Framework::CellTrsGeoBuilder> m_geoBuilder;
  
  /// face builder
  Framework::GeometricEntityPool<Framework::FaceTrsGeoBuilder> m_faceBuilder;
  
  /// Error L1 Bx
  CFreal m_ErrorBxL1;

  /// Error L1 By
  CFreal m_ErrorByL1;

  /// Error L1 Ez
  CFreal m_ErrorEzL1;

  /// Error L1 Rho
  CFreal m_ErrorRhoL1;

  /// Error L1 Ux
  CFreal m_ErrorUxL1;
  
  /// Error L1 Uy
  CFreal m_ErrorUyL1;

  /// Error L1 T
  CFreal m_ErrorTL1;

  /// Error non-dimensional total
  CFreal m_ErrorTotalL1;

  /// Error L2 BX
  CFreal m_ErrorBxL2;

  /// Error L2 By
  CFreal m_ErrorByL2;

  /// Error L2 Ez
  CFreal m_ErrorEzL2;

  /// Error L2 Rho
  CFreal m_ErrorRhoL2;

  /// Error L2 Ux
  CFreal m_ErrorUxL2;

  /// Error L2 Uy
  CFreal m_ErrorUyL2;

  /// Error L2 T
  CFreal m_ErrorTL2;

  /// Error non-dimensional total
  CFreal m_ErrorTotalL2;
  
}; // end of class GridConvergence

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolumeArcJet

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolumeMultiFluidMHD_GridConvergence_hh



