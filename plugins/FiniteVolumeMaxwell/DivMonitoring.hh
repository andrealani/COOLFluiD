#ifndef COOLFluiD_Numerics_FiniteVolumeMaxwell_DivMonitoring_hh
#define COOLFluiD_Numerics_FiniteVolumeMaxwell_DivMonitoring_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/DataProcessingData.hh"
#include "Framework/DataSocketSink.hh"
#include "Environment/FileHandlerOutput.hh"

//////////////////////////////////////////////////////////////////////////////
 
namespace COOLFluiD {
  
  namespace Framework {
    class CellTrsGeoBuilder;
  }
  
  namespace Physics {
    namespace Maxwell {
      class MaxwellVarSet;
    }
  }
  
  namespace Numerics {
    
    namespace FiniteVolumeMaxwell {

//////////////////////////////////////////////////////////////////////////////

/**
 *
 * This class computes the divergence of magnetic and electric fields
 *
 * @author Alejandro Alvarez Laguna
 *
 */
class DivMonitoring : public Framework::DataProcessingCom {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */

  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor
   */
  DivMonitoring(const std::string& name);

  /**
   * Default destructor
   */
  ~DivMonitoring();

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
   * Compute the divergence of B and E
   */  
  void computeDivergence(); 
  
  /**
   * Compute the Analitical Solution of Testcases
   */  
  void computeAnaliticalSolution();  
  
  /**
   * Compute Cylindrical components of Magnetic Field
   */  
  void computeCylindricalComponents();   

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
  
  
   /**
   * Prepare the output file for writing
   */
  void prepareOutputFile(std::ofstream& outputFile);

  /**
   * Write the electric and magnetic in the output file
   */
  void writeOutputFile();

  /**
   * Construct the file path to which to write the Divergence term
   */
  boost::filesystem::path constructFilename();

//   /**
//    * Outputs the quality of the cells to a file
//    */
//   void printDivMonitoringToFile();
  
  std::string outputFileDivMonitoring() const
  {
    return m_tecplotFileDivMonitoring;
  }
  int outputFileDivMonitoring_SaveRate() const
  {
    return m_tecplotFileDivMonitoringSaveRate;
  }

private: //data

  /// socket for average B values in x-direction at the cell faces 
  Framework::DataSocketSource<CFreal> socket_avgBxFace;
  
  /// socket for average B values in y-direction at the cell faces 
  Framework::DataSocketSource<CFreal> socket_avgByFace;
  
  /// socket for average B values in z-direction at the cell faces 
  Framework::DataSocketSource<CFreal> socket_avgBzFace;
  
    /// socket for average E values in x-direction at the cell faces 
  Framework::DataSocketSource<CFreal> socket_avgExFace;
  
  /// socket for average E values in y-direction at the cell faces 
  Framework::DataSocketSource<CFreal> socket_avgEyFace;
  
  /// socket for average E values in z-direction at the cell faces 
  Framework::DataSocketSource<CFreal> socket_avgEzFace;
  
  /// storage of divB
  Framework::DataSocketSource <CFreal> socket_divB;
  
  /// storage of curlB
  Framework::DataSocketSource <CFreal> socket_curlB;  
  
  /// storage of divE
  Framework::DataSocketSource <CFreal> socket_divE;
  
  /// storage of theta
  Framework::DataSocketSource <CFreal> socket_theta;  
  
  /// storage of radial component of Magnetic field
  Framework::DataSocketSource <CFreal> socket_Bradial; 
  
  /// storage of azimuthal component of Magnetic field
  Framework::DataSocketSource <CFreal> socket_Btheta;   
  
  /// storage of azimuthal component of Magnetic field
  Framework::DataSocketSource <CFreal> socket_BthetaTheory;    
  
  /// storage of azimuthal component of Magnetic field
  Framework::DataSocketSource <CFreal> socket_ExPWTheory;  

  /// storage of azimuthal component of Magnetic field
  Framework::DataSocketSource <CFreal> socket_EyPWTheory;  

  /// storage of azimuthal component of Magnetic field
  Framework::DataSocketSource <CFreal> socket_EzPWTheory; 

  /// storage of azimuthal component of Magnetic field
  Framework::DataSocketSource <CFreal> socket_BxPWTheory;  

  /// storage of azimuthal component of Magnetic field
  Framework::DataSocketSource <CFreal> socket_ByPWTheory;  

  /// storage of azimuthal component of Magnetic field
  Framework::DataSocketSource <CFreal> socket_BzPWTheory; 


  /// storage of azimuthal component of Magnetic field
  Framework::DataSocketSource <CFreal> socket_Error;    
  
  /// storage of volumes
  Framework::DataSocketSink<CFreal> socket_volumes;
  /// storage of nodes
  Framework::DataSocketSink < Framework::Node* , Framework::GLOBAL > socket_nodes;

  /// storage of states
  Framework::DataSocketSink < Framework::State* , Framework::GLOBAL > socket_states;
  
  /// storage for State's
  Framework::DataSocketSink < Framework::State* > socket_gstates;
  
  /// storage of nstates (states in nodes)
  Framework::DataSocketSink<RealVector> socket_nstates;  
  
  /// storage of isOutward  
  Framework::DataSocketSink<CFint> socket_isOutward;
  
  /// storage of normals
  Framework::DataSocketSink<CFreal> socket_normals; 
  
 /// storage of ErrorExPW
 Framework::DataSocketSource <CFreal> socket_ErrorExPW;  

 /// storage of ErrorEyPW
 Framework::DataSocketSource <CFreal> socket_ErrorEyPW; 

 /// storage of ErrorEzPW
 Framework::DataSocketSource <CFreal> socket_ErrorEzPW;

 /// storage of ErrorBxPW
 Framework::DataSocketSource <CFreal> socket_ErrorBxPW;  

 /// storage of ErrorByPW
 Framework::DataSocketSource <CFreal> socket_ErrorByPW; 

 /// storage of ErrorBzPW
 Framework::DataSocketSource <CFreal> socket_ErrorBzPW;
 
//  /// storage of divBL2Norm
//  Framework::DataSocketSource <CFreal> socket_divBL2Norm;
//  
//  /// storage of divEL2Norm
//  Framework::DataSocketSource <CFreal> socket_divEL2Norm;
//  
//  /// storage of TEPWErrorL2Norm
//  Framework::DataSocketSource <CFreal> socket_TEPWErrorL2Norm; 
  
  /// builder of cells geometric entities
  Framework::GeometricEntityPool<Framework::CellTrsGeoBuilder> m_geoBuilder;  
  
  /// builder of face geometric entities
  Framework:: GeometricEntityPool<Framework::FaceTrsGeoBuilder> m_faceBuilder;
  
  /// update variable set
  Common::SafePtr<Physics::Maxwell::MaxwellVarSet> m_updateVarSet;
  
  ///store the Norm of divB
  RealVector m_divBL2Norm;

  ///store the Norm of divE
  RealVector m_divEL2Norm;

  ///store the Norm of Error in Transverse Electric Planar Wave
  RealVector m_TEPWErrorL2Norm;

  ///store the Norm of Error in Transverse Magnetic Planar Wave
  RealVector m_TMPWErrorL2Norm;
  
  ///store the Norm of Error in full Planar Wave
  RealVector m_fullPWErrorL2Norm;  

//  /// array storing the B Field divergence in the Cells
//  RealVector m_divB;  
//  
//  /// array storing the E Field divergence in the Cells
//  RealVector m_divE;    
    
  /// save rate
  CFuint  m_saveRate;

  /// Save each iteration to a different Tecplot file (with suffix _iter#).
  bool  m_appendIter;

  /// Save each timestep to a different Tecplot file (with suffix _iter#).
  bool  m_appendTime;

  /// name of the Output File where to write the electric current
  std::string m_nameOutputFileDivMonitoring;
  
  /// electric and Magnetic Field Divergence output file
  std::string m_tecplotFileDivMonitoring;
  int m_tecplotFileDivMonitoringSaveRate;
  
  /// Save each timestep to a different Tecplot file (with suffix _iter#).
  std::string m_toRun;  
  
}; // end of class DivMonitoring

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolumeMaxwell

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolumeMaxwell_DivMonitoring_hh
