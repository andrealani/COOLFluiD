#ifndef COOLFluiD_Physics_Parade_ParadeLibrary_hh
#define COOLFluiD_Physics_Parade_ParadeLibrary_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/RadiationLibrary.hh"
#include "MathTools/RealMatrix.hh"
#include "MathTools/RealVector.hh"
#include "Common/OSystem.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Environment {
    class FileHandlerInput;
    class FileHandlerOutput;
  }
  
  namespace Framework {
    class PhysicalChemicalLibrary;
  }
  
  namespace Parade {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents the interface for a ParadeLibrary.
 *
 * @author Andrea Lani
 * @author Alessandro Munafo'
 *
 */
class ParadeLibrary : public Framework::RadiationLibrary {
public:
  
  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);
  
  /**
   * Constructor without arguments
   */
  ParadeLibrary(const std::string& name);
  
  /**
   * Default destructor
   */
  ~ParadeLibrary();
  
  /**
   * Configures this configurable object.
   */
  virtual void configure ( Config::ConfigArgs& args );
  
  /**
   * Setups the data of the library
   */
  void setup();
  
  /**
   * Unsetups the data of the library
   */
  void unsetup();
  
  /// Get the number of loops neede to cover the whole spectrum
  virtual CFuint getWavLoopSize() const;
  
  /// Run the radiation code on the stagnation line
  virtual void runOnStagnationLine(Common::SafePtr<std::vector<CFuint> > stagnationLineCells,
				   Framework::ProxyDofIterator<CFreal>* pstates,
				   CFreal* qrad);
  
  /// Run the radiation code on a structured mesh
  virtual void runOnStructuredMesh(const std::vector<std::vector<CFuint>* >& meshByLine,
				   Framework::ProxyDofIterator<CFreal>* pstates,
				   CFreal* qrad);
  
  /// Compute radiative properties
  /// @param pstates    array fo flowfield solution (State's)
  /// @param data       matrix storing emission and absorption coefficients for each wavelength (row ID):
  ///                   data(index of spectral point, index of state*2)   = emission coeff
  ///                   data(index of spectral point, index of state*2+1) = absorption coeff
  /// @para, iWavRange  index of the loop over wavelengths ranges
  virtual void computeProperties(Framework::ProxyDofIterator<CFreal>* pstates,
				 RealMatrix& data, CFuint iWavRange);
  
private:
  
  /// Set up the library sequentially
  void setLibrarySequentially();
  
  /// update the range of wavelengths to use
  void updateWavRange(CFuint iWavRange);
  
  /// write the data (grid, temperatue, densities) corresponding to the local mesh
  void writeLocalData(Framework::ProxyDofIterator<CFreal>* pstates);
  
  /// read the radiative coefficients corresponding to the local mesh
  void readLocalRadCoeff(RealMatrix& data);
  
  /// run PARADE
  void runLibrary() const
  {
    std::string command = "./parade > outfile";
    Common::OSystem::getInstance().executeCommand(command);
  }
  
  /// run PARADE in parallel
  void runLibraryInParallel() const
  {
    std::string command = "cd " + m_paradeDir + " ; ./parade > outfile ; cd -";
    Common::OSystem::getInstance().executeCommand(command);
  }
  
  /// Read value from Parade input file
  template <typename T>
  void readValue(std::string& line, const std::string& name, T& value)
  {
    size_t posS = line.find(name);
    if (posS != std::string::npos) {
      line.erase(posS);
      Common::StringOps::trim(line);
      value = Common::StringOps::template from_str<T>(line);
    }
  }
  
  /// AL: what follows has to be VERIFIED
  
  /// write the stagnation line data
  void writeStagnationLineData(Common::SafePtr<std::vector<CFuint> > stagnationLineCells,
			       Framework::ProxyDofIterator<CFreal>* pstates);
  
  /// write all the data
  void writeAllData(const std::vector<std::vector<CFuint>* >& meshByLine,
		    Framework::ProxyDofIterator<CFreal>* pstates);
  
  /// compute Qrad source term for oprically thin case
  void computeQradEmDom(CFreal* qrad);
  
  /// read the ouput files for the stagnation line 
  void computeQradStagLine(CFreal* qrad);
  
  /// read the ouput files for the whole Mesh 
  void computeQradMesh(const std::vector<std::vector<CFuint>* >& meshByLine,
		       Framework::ProxyDofIterator<CFreal>* pstates,
		       CFreal* qrad);
  
  /// computhe the wall heat flux
  void computeRadWallHeatFlux(const std::vector<std::vector<CFuint>* >& meshByLine,
		              Framework::ProxyDofIterator<CFreal>* pstates);
  
  /// Planck function 
  CFreal Planck(const CFreal& lambda, const CFreal& T);
  
  /// Exponential integral of first order E1 
  CFreal E1(const CFreal& tau);
  
  /// Exponential integral of second order E2 
  CFreal E2(const CFreal& tau);
  
  /// Exponential integral of third order E3 
  CFreal E3(const CFreal& tau);
  
private: 
  
  /// input file handle
  Common::SelfRegistPtr<Environment::FileHandlerInput> m_inFileHandle;
  
  /// output file handle
  Common::SelfRegistPtr<Environment::FileHandlerOutput> m_outFileHandle;
  
  /// directory where Parade is launched
  std::string m_paradeDir; 
  
  /// path to the grid file
  boost::filesystem::path m_gridFile;
  
  /// path to the temperature file
  boost::filesystem::path m_tempFile;
  
  /// path to the densities file
  boost::filesystem::path m_densFile;
  
  /// path to the radiative properties file
  boost::filesystem::path m_radFile;
  
  /// thermodynamic library
  Common::SafePtr<Framework::PhysicalChemicalLibrary> m_library;
  
  /// number of points in the spectrum
  CFuint m_spectrumSize;
    
  /// roto-translational temperature ID
  CFuint m_trTempID;
  
  // electron temperature ID
  CFuint m_elTempID;
  
  // vibrational temperature ID
  CFuint m_vibTempID;
  
  /// array with wavelengths
  std::vector<CFreal> m_wavelengths;
  
  /// array with absorption coefficients
  RealMatrix _kNu;

  /// array with emission coefficients
  RealMatrix _emNu;

  /// array with optical thicknesses at given location and given wavelength
  RealMatrix _tauNu;

  /// array with source term at given location and given wavelength
  RealMatrix _sourceNu;
  
  /// array with molar masses
  RealVector m_mmasses;
  
  /// array with Avogadro number/molar masses
  RealVector m_avogadroOvMM;
    
  /// flag array to indicate molecular species 
  std::vector<bool> m_molecularSpecies;
  
  /// flag for emission dominated cases 
  bool flag_Em;
 
  /// flag for computing heat flux according to the tangent slab approximation 
  bool flag_HeatFlux;   
 
  /// Value at boundary 1 of temperature for solving the RTE by means of tangent slab method
  CFreal _Tb1;

  /// Value at boundary 2 of temperature for solving the RTE by means of tangent slab method
  CFreal _Tb2;

  /// Under relaxation factor to applied to the radiative source term. If not specified by the 
  //  user it is set to 1.0. 
  CFreal _uFac;
  
  /// minimum number density
  CFreal m_ndminFix;
  
  /// minimum temperature
  CFreal m_tminFix;
  
  /// name of the temporary local directory where Parade is run
  std::string m_localDirName;
  
  /// Reuse existing radiative data (requires the same number of processors as in the previous run).
  CFuint m_reuseProperties;
  
}; // end of class ParadeLibrary

//////////////////////////////////////////////////////////////////////////////

    } // namespace Parade

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_Parade_ParadeLibrary_hh
