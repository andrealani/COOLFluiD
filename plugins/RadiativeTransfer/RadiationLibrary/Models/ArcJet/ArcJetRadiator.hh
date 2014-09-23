#ifndef COOLFluiD_Physics_ArcJetRadiator_hh
#define COOLFluiD_Physics_ArcJetRadiator_hh

///////////////////////////////////////////////////////////////////////////
#include "RadiativeTransfer/RadiationLibrary/Radiator.hh"
#include "Framework/DataProcessingData.hh"
#include "Framework/DataSocketSink.hh"
#include "Framework/DataSocketSource.hh"
#include "Framework/CellTrsGeoBuilder.hh"
#include "Framework/GeometricEntityPool.hh"
#include "Framework/DofDataHandleIterator.hh"
#include "Framework/PhysicalChemicalLibrary.hh" 
#include "Common/OSystem.hh"
//////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  
  namespace Environment {
    class FileHandlerInput;
  }
  
 namespace RadiativeTransfer{ 

//////////////////////////////////////////////////////////////////////////

/**
 * This class represents the interface for a ArcJet Radiative Table
 *
 * @author Pedro Santos
 * @author Alejandro Laguna
 *
 */
class ArcJetRadiator : public Radiator {
public:
  
  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);
  
  /**
   * Constructor without arguments
   */
  ArcJetRadiator(const std::string& name);
  
  static std::string getClassName() { return "ArcJetRadiator"; }

  /**
   * Default destructor
   */
  ~ArcJetRadiator();
  
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
  
private:

  /// update the range of wavelengths to use
  void updateWavRange(CFreal wavMin, CFreal wavMax);
  
   /**
   * Reads the binary file containing the opacities as function of the temperature, 
   * pressure and wavelength  
   */
  void readOpacities();
  
  /**
   * Interpolates the values of the opacity tables
   */ 
  void tableInterpolate(CFreal T, CFreal p, CFuint ib, CFreal& val1, CFreal& val2);


 /**
   * Fill the m_data matrix with the interpolated values from the tables
   */ 
  void genData();

private: 
  
  //Vector containing the states id of the TRS;
  std::vector<CFuint> m_statesID;

  /// array storing absorption and emission coefficients
  /// wavelength       = m_radCoeff(local state ID, spectral point idx*3)
  /// emission coeff   = m_radCoeff(local state ID, spectral point idx*3+1)
  /// absorption coeff = m_radCoeff(local state ID, spectral point idx*3+2)
  RealMatrix m_data;

  /// name of the temporary local directory where Parade is run
  std::string m_localDirName;

  /// Opactities read from table. Stored as follows:
  // m_opacities[pressure][remperature][bin]
  std::vector< std::vector < std::vector< double > > >  m_opacities;
  
  /// Radiative source reaf from table and stored the same way as the opacities
  std::vector< std::vector < std::vector< double > > >  m_radSource;
  
  Common::SafePtr<Framework::DofDataHandleIterator<CFreal, Framework::State, Framework::GLOBAL> > m_pstates;
  std::auto_ptr<Framework::ProxyDofIterator<CFreal> > m_nstatesProxy;


  /// thermodynamic library
  Common::SafePtr<Framework::PhysicalChemicalLibrary> m_library;

  /// min, max and delta wavelenght for the current spectral loop
  CFreal m_wavMin, m_wavMax, m_dWav;

  ///Vector with the spectral loop powers of the states
  std::vector<CFreal> m_spectralLoopPowers;

  //vector with the comulative probability distributions
  std::vector<CFreal> m_cpdEms;

  /// input file handle
  Common::SelfRegistPtr<Environment::FileHandlerInput> m_inFileHandle;
  
  /// opacities file
  boost::filesystem::path m_binTableFile;
 
  /// name of the .dat binary table with the opacities
  std::string m_binTabName;

  /// storage of the temperatures of the opacity table. 
  std::vector<CFreal> m_Ttable;
  
  /// storage of the pressure of the opacity table
  std::vector<CFreal> m_Ptable; 
 
  /// storage of the wavelengths of the opacity table
  RealVector m_wavTable; 
 
  /// number of bins
  CFuint m_nbBins;
  
  /// number of Temperatures
  CFuint m_nbTemp;
  
  /// number of pressures
  CFuint m_nbPress;
 
public: // Radiator interface

  void setupSpectra(CFreal wavMin, CFreal wavMax);

  CFreal getEmission(CFreal lambda, RealVector &s_o);

  CFreal getAbsorption(CFreal lambda, RealVector &s_o);

  CFreal getSpectaLoopPower();

  void computeEmissionCPD();

  void getRandomEmission(CFreal &lambda, RealVector &s_o);
  
  void getData();

  inline void getSpectralIdxs(CFreal lambda, CFuint *idx);

  // ConfigObject interface
}; // end of class 

//////////////////////////////////////////////////////////////////////////////

    } // namespace Parade

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
#endif
