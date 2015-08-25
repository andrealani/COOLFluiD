#ifndef COOLFluiD_Numerics_FiniteVolumeRadiation_Radiation_hh
#define COOLFluiD_Numerics_FiniteVolumeRadiation_Radiation_hh

///////////////////////////////////////////////////////////////////////////

#include "Framework/DataProcessingData.hh"
#include "Framework/DataSocketSink.hh"
#include "Framework/DataSocketSource.hh"
#include "Framework/CellTrsGeoBuilder.hh"
#include "Framework/GeometricEntityPool.hh"

#include "Common/OSystem.hh"

//////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  
  namespace Environment {
    class FileHandlerInput;
  }
  
  namespace Framework {
    class PhysicalChemicalLibrary;
  }

  namespace Numerics {

    namespace FiniteVolumeRadiation {

//////////////////////////////////////////////////////////////////////////

/**
 * This class compute the radiative heat transfer using a Finite Volume algorithm
 *
 * @author Alejandro Alvarez (C++ version of Alan Wray's algorithm)
 */
class Radiation : public Framework::DataProcessingCom {
public:

  /**
   * Constructor
   */
  Radiation(const std::string& name);

  /**
   * Default destructor
   */
  ~Radiation();
  
  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);  

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
   * Compute the directions depending on the option selected 
   */  
  void getDirections();
  
  /**
   * Compute the advance order depending on the option selected 
   */  
  void getAdvanceOrder(); 
  
  /**
   * Compute the advance order depending on the option selected 
   */  
  void getFieldOpacities(CFuint ib);
    
  /**
   * Compute the advance order depending on the option selected 
   */  
  void getFieldOpacities(const CFuint ib, const CFuint iCell);
  
  /**
   * Reads the binary file containing the opacities as function of the temperature, pressure
   * and wavelength  
   */
  void readOpacities();
  
  /**
   * Compute the order of the cells depending on the option selected (For the moment in the execute())
   */
  void findOrderOfAdvance();
  
  /**
   * Interpolates the values of the opacity tables
   */ 
  void tableInterpolate(CFreal T, CFreal p, CFuint ib, CFreal& val1, CFreal& val2);
  
  /**
   * Writes radial q and divQ to a file for the Sphere case 
   * (only if option RadialData is enabled)
   * 
   */
  void writeRadialData();

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
  
    /// storage of the the stage of the order of advance
  Framework::DataSocketSource <CFreal> socket_CellID;
  
  /// storage of the divq 
  Framework::DataSocketSource <CFreal> socket_divq;
  
  /// storage of the qx 
  Framework::DataSocketSource <CFreal> socket_qx;
  
  /// storage of the qy
  Framework::DataSocketSource <CFreal> socket_qy;
  
  /// storage of the qz
  Framework::DataSocketSource <CFreal> socket_qz;
//   
  /// storage of the radiative source table. Source[ib](it,ip)
  Framework::DataSocketSource <CFreal> socket_TempProfile; 
  
  /// pointer to the physical-chemical library
  Common::SafePtr<Framework::PhysicalChemicalLibrary> m_library; 
  
  /// builder of geometric entities
  Framework::GeometricEntityPool<Framework::CellTrsGeoBuilder> m_geoBuilder;
  
  /// face builder
  Framework::GeometricEntityPool<Framework::FaceTrsGeoBuilder> m_faceBuilder;
  
  /// temporary normal to the face
  RealVector m_normal; 
  
  /// Weights for the directions
  RealVector m_weight;
  
  /// Field source of oppacity table 
  RealVector m_fieldSource;
  
  /// Field Absorption of oppacity table used if exponential Method
  RealVector m_fieldAbsor;
  
  /// Field Absorption of oppacity table used if not Exponential Method
  RealVector  m_fieldAbSrcV;
  
  /// Field Absorption of oppacity table used if not Exponential Method
  RealVector  m_fieldAbV;  
  
  /// Exponent for the radiation of oppacity table
  RealVector m_In;
  
  /// Exponent for the radiation of oppacity table
  RealVector m_II;
  
  /// Opactities read from table. Stored as follows:
  /// kappa(p0,b0,T0),kappa(p0,b0,T1),...,kappa(p0,b0,Tn),
  /// kappa(p0,b1,T0),kappa(p0,b1,T1),...,kappa(p0,b1,Tn), .... 
  /// kappa(p0,bn,T0),kappa(p0,bn,T1),...,kappa(p0,bn,Tn), ....
  /// kappa(p1,b0,T0),kappa(p1,b0,T1),...,kappa(p1,b0,Tn), ...
  RealVector m_opacities;
  
  /// Radiative source reaf from table and stored the same way as the opacities
  RealVector m_radSource;
  
  /// storage of the temperatures of the opacity table. 
  RealVector m_Ttable;
  
  /// storage of the pressure of theopacity table
  RealVector m_Ptable; 
  
  /// Done status of a cell in a given direction at the end of a stage
  std::vector<bool> m_sdone;
  
  /// Done status of a cell in a given direction including those done in the current stage
  std::vector<bool> m_cdone;
  
  /// Directions 
  RealMatrix m_dirs;  
  
  /// Create the array advance_order, which contains, for each direction (1st index), the list of cells to be advanced.  The list is divided into "stages" that consist of 
  /// a set of cells that can be advanced in parallel.These sets are terminated by a negative entry.  E.g., if there are 8 cells and advance_order(1:8,1) = (1,2,-5,3,4,-8,6,-7), 
  /// then cells (1,2,5) can be done first, and are in fact the boundary cells for direction 1,
  /// then cells (3,4,8) can be done; finally cells (6,7) can be done to complete the sweep in direction 1.
  std::vector<std::vector<int> > m_advanceOrder;   
  
  /// Radiative flux
  RealMatrix m_q;
  
  /// Divergence of dif flux
  RealVector m_divq;
  
  /// Number of directions types
  CFuint m_nDirTypes;   
  
  /// user defined number of directions
  CFuint m_nDirs;
  
  /// option to print the radial q and divQ for the Sphere
  bool m_radialData;
  
  /// old algorithm just kept for comparison purposes
  bool m_oldAlgo;
  
  /// user defined number of points in the radial direction
  CFuint m_Nr;
  
  /// name of the temporary local directory where Parade is run
  boost::filesystem::path m_dirName;
  
  /// name of the .dat binary table with the opacities
  std::string m_binTabName;
  
  /// bool to write the table in a file
  bool m_writeToFile;
  
  /// File where the table is written
  std::string m_outTabName;
  
  /// Use exponential method
  bool m_useExponentialMethod;
  
  /// input file handle
  Common::SelfRegistPtr<Environment::FileHandlerInput> m_inFileHandle;
  
  /// opacities file
  boost::filesystem::path m_binTableFile;
  
  /// number of bins
  CFuint m_nbBins;
  
  /// number of Temperatures
  CFuint m_nbTemp;
  
  /// number of pressures
  CFuint m_nbPress;
  
  /// Radial average of q vector for a Sphere
  RealVector m_qrAv;
  
  /// Radial average of divQ for a Sphere
  RealVector m_divqAv;  
  
  /// constant pressure
  CFreal m_constantP;
  
  /// minimum temparature
  CFreal m_Tmin;
  
  /// maximum temperature
  CFreal m_Tmax;
  
  /// temperature interval
  CFreal m_deltaT;
  
  /// ID of pressure in the state vector
  CFuint m_PID;
  
  /// ID of temperature in the state vector
  CFuint m_TID;
  
}; // end of class Radiation

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolumeRadiation

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolumeRadiation_Radiation_hh



