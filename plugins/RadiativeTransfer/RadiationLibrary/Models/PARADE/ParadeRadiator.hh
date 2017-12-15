#ifndef COOLFluiD_RadiativeTransfer_ParadeRadiator_hh
#define COOLFluiD_RadiativeTransfer_ParadeRadiator_hh

//////////////////////////////////////////////////////////////////////////////
#include "RadiativeTransfer/RadiationLibrary/Radiator.hh"
#include "MathTools/RealMatrix.hh"
#include "MathTools/RealVector.hh"
#include "Common/OSystem.hh"
#include "boost/filesystem.hpp"
#include "Framework/ProxyDofIterator.hh"
#include "Framework/DofDataHandleIterator.hh"
#include "Common/StringOps.hh"
#include "Framework/PhysicalChemicalLibrary.hh"
#include "Common/CFMap.hh"
///////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Environment {
    class FileHandlerInput;
    class FileHandlerOutput;
  }
  
  namespace RadiativeTransfer {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents the interface for a ParadeLibrary.
 *
 * @author Pedro Santos
 * @author Andrea Lani
 * @author Alessandro Munafo'
 *
 */
class ParadeRadiator : public Radiator {
public:
  
  /// Nested class providing some functions that can be called by both CPU and GPU
  class DeviceFunc {
  public:
    /// Default constructor
    HOST_DEVICE DeviceFunc() {}
    
    /// Destructor
    HOST_DEVICE ~DeviceFunc() {} 
    
    /// Compute all binned data corresponding to the given cell
    /// The template parameter prevents from having linking errors
    template <DeviceType DT>
    HOST_DEVICE void computeCellBins(const bool isEquil,
				     const CFuint nbPoints,
				     const CFuint i, 
				     const CFuint j, 
				     const CFuint nbBinsre, 
				     const CFuint testID,
				     const CFreal temp,
				     const CFreal dWavIn,
				     const CFreal* vctBins,  
				     CFreal* data,
				     CFreal* alpha_bin,
				     CFreal* emission_bin,
				     CFreal* B_binCurr);
    
    /// compute Planck function
    /// @param wav is wavelength in Angstrom
    HOST_DEVICE CFreal computePlanckFun(const CFreal h, 
					const CFreal c, 
					const CFreal k_b,
					const CFreal wav,
					const CFreal temp) const
    {
      return 2.*h*pow(c,2.)/pow(wav*1e-10,5.)*(1./(exp(h*c/(wav*1e-10)*k_b*temp))-1.);
    }
    
  };
  
  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);
  
  /**
   * Constructor without arguments
   */
  ParadeRadiator(const std::string& name);
  
  static std::string getClassName() { return "ParadeRadiator"; }

  /**
   * Default destructor
   */
  ~ParadeRadiator();
  
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
  
  void setupSpectra(CFreal wavMin, CFreal wavMax);
  
  CFreal getEmission(CFreal lambda, RealVector &s_o);
  
  CFreal getAbsorption(CFreal lambda, RealVector &s_o);
  
  CFreal getSpectraLoopPower();
  
  void computeEmissionCPD();
  
  void getRandomEmission(CFreal &lambda, RealVector &s_o);

  void getData();
  
  inline void getSpectralIdxs(CFreal lambda, CFuint& idx1, CFuint& idx2);

private:

  void setLibrarySequentially();
  
  /// update the range of wavelengths to use
  void updateWavRange(CFreal wavMin, CFreal wavMax);
  
  /// write the data (grid, temperatue, densities) corresponding to the local mesh
  void writeLocalData();
  
  /// read the radiative coefficients corresponding to the local mesh
  void readLocalRadCoeff();
  
  /// write the local mesh radiative coefficients to a ASCII file
  void writeLocalRadCoeffASCII(const CFuint nbCells);
  
  /// run PARADE
  void runLibrary() const
  {
    std::string command = "./parade > outfile";
    Common::OSystem::getInstance().executeCommand(command);
  }
  
  /// run PARADE in parallel
  void runLibraryInParallel() const
  { 
    CFLog(VERBOSE, "ParadeRadiator::runLibraryInParallel()\n");
    std::string command = "cd " + m_paradeDir.string() + " ; ./parade > outfile ; cd -";
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
  
  /// @return Flag telling if the process stores the full grid in parallel
  bool fullGridInProcess() const {return (m_namespace != "Default");}
  
  /// apply the binning method to reduce the spectral data
  void computeBinning();
  
  /// apply the banding method to reduce the spectral data
  void computeBanding();
    
  /// apply the binning/banding method to reduce the spectral data
  void computeBinningBanding();

  /// compute recv counts and dispacements for parallel communication
  void computeRecvCountsDispls(const CFuint totalNbCells, const CFuint sizeCoeff, 
			       CFuint& minSizeToSend, CFuint& maxSizeToSend,
			       std::vector<int>& recvCounts, std::vector<int>& displs);
  
 
  
  /// compute the averaged bins/bands
  virtual void computeAveragedBins
  (const CFuint nbBinsre, 
   const CFuint testID,
   Framework::LocalArray<CFreal>::TYPE& vctBins);
  
private: 
  
  /// Rank in the corresponding MPI group (namespace)
  CFuint m_rank;
  
  /// Number of ranks in the corresponding MPI group (namespace)
  CFuint m_nbProc;
  
  //Vector containing the states id of the TRS;
  std::vector<CFuint> m_statesID;

  /// input file handle
  Common::SelfRegistPtr<Environment::FileHandlerInput> m_inFileHandle;
  
  /// output file handle
  Common::SelfRegistPtr<Environment::FileHandlerOutput> m_outFileHandle;

  /// File where the table is written
  std::string m_outTabName;

  /// Bool expression to write the table in a file
  bool m_writePARADEToFile;
  
  /// directory where Parade is launched
  boost::filesystem::path m_paradeDir; 

  //name of the temporary local directory where Parade is run
  boost::filesystem::path m_dirName;
  
  /// path to the grid file
  boost::filesystem::path m_gridFile;
  
  /// path to the temperature file
  boost::filesystem::path m_tempFile;
  
  /// path to the densities file
  boost::filesystem::path m_densFile;
  
  /// path to the radiative properties file
  boost::filesystem::path m_radFile;
  
  /// roto-translational temperature ID
  CFuint m_trTempID;
  
  /// electronic temperature ID
  CFuint m_elTempID;
  
  /// vibrational temperature ID
  CFuint m_vibTempID;

  /// array storing absorption and emission coefficients
  /// wavelength       = m_radCoeff(local state ID, spectral point idx*3)
  /// emission coeff   = m_radCoeff(local state ID, spectral point idx*3+1)
  /// absorption coeff = m_radCoeff(local state ID, spectral point idx*3+2)
  Framework::LocalArray<CFreal>::TYPE m_data;
  
  /// vector storing the averaged absorption coefficient for each wavelength
  std::vector<CFreal> m_alphaav;
  
  /// vectors storing the matrix with the data for each bin and each cell, for the absorption, emission and source terms
  Framework::LocalArray<CFreal>::TYPE m_alpha_bin;
  Framework::LocalArray<CFreal>::TYPE m_emission_bin;
  Framework::LocalArray<CFreal>::TYPE m_B_bin;

  /// vector storing the averaged absorption coefficients data
  std::vector<CFreal> m_alpha_avbin;
  
  /// name of the temporary local directory where Parade is run
  std::string m_localDirName;
  
  /// namespace within which PARADE is run in parallel
  std::string m_namespace;
  
  /// Reuse existing radiative data (requires the same number of processors as in the previous run).
  bool m_reuseProperties;
  
  /// number of spectral points
  CFuint m_nbPoints;
  
  /// iterator for the state vector
  Common::SafePtr<Framework::DofDataHandleIterator<CFreal, Framework::State, Framework::GLOBAL> > m_pstates;
  
  /// thermodynamic library
  Common::SafePtr<Framework::PhysicalChemicalLibrary> m_library;
  
  /// minimum number density
  CFreal m_ndminFix;
  
  /// minimum temperature
  CFreal m_TminFix;
  
  /// bands' distribution
  std::string m_bandsDistr;
 
  /// array with molar masses
  RealVector m_mmasses;

  /// array with Avogadro number/molar masses
  RealVector m_avogadroOvMM;
  
  /// Path to Parade's binary files
  std::string m_libPath;
  
  /// min, max and delta wavelenght for the current spectral loop
  CFreal m_wavMin, m_wavMax, m_dWav;

  ///Vector with the spectral loop powers of the states
  std::vector<CFreal> m_spectralLoopPowers;
  
  /// vector with the comulative probability distributions
  std::vector<CFreal> m_cpdEms;
  
  /// flag telling whether the input flowfield is in LTE
  bool m_isLTE;

  bool m_Equilibrium;

  /// flag telling whether the input flowfield has only the mass fractions as data
  bool m_massfraction;

  bool m_TisLTE;
  
  /// flag telling whether the binning has to be applied
  bool m_binning;
  
  /// flag telling whether the banding has to be applied
  bool m_banding;
  
  /// flag telling to write the radiative coefficients to ASCII file 
  bool m_writeLocalRadCoeffASCII;
  
  /// flag telling to parallelize as much as possible to save memory
  bool m_saveMemory;
  
  /// flag array to indicate molecular species
  std::vector<bool> m_molecularSpecies;

  /// bool to write the table in a file
  bool m_writeHSNB;
  
}; // end of class ParadeRadiator

//////////////////////////////////////////////////////////////////////////////

template <DeviceType DT>
void ParadeRadiator::DeviceFunc::computeCellBins(const bool isEquil,
						 const CFuint nbPoints,
						 const CFuint i, 
						 const CFuint j, 
						 const CFuint nbBinsre, 
						 const CFuint testID,
						 const CFreal temp,
						 const CFreal dWavIn,
						 const CFreal* vctBins,  
						 CFreal* data,
						 //const CFreal* data,
						 CFreal* alpha_bin,
						 CFreal* emission_bin,
						 CFreal* B_binCurr)
{
  const CFuint nbCols = nbPoints*3;
  const CFreal emCoef = data[j*nbCols + i*3+1];
  
  // AL: just to reproduce the BUG
  CFreal& abCoef = data[j*nbCols + i*3+2];
  
  const CFreal h = 6.626070040e-34;     //SI units Js
  const CFuint c = 3e08; //SI units m/s
  const CFreal k_b = 1.3806485279e-23;  //SI units J/K
  CFreal Bs = 0.;
  
  if(!isEquil){ 
    if(abCoef=!0){ // this line is a clear BUG!!!!
      Bs = emCoef/abCoef;
    }
    else{
      Bs = computePlanckFun(h,c,k_b,data[j+nbCols+i*3],temp);
    }
  }
  else{
    Bs = computePlanckFun(h,c,k_b,data[j+nbCols+i*3],temp);
  }
  
  const CFuint idx0 = nbBinsre*j;
  
  /*
    alpha_bin[idx0]    = 0.;
    emission_bin[idx0] = 0.;
    B_binCurr[idx0]    = 0.;
  */
  const CFreal test   = data[j*nbCols + i*3+testID];
  
  for(CFuint k=1;k<nbBinsre;++k) {
    CFreal alpha_bincoeff = 0.;
    CFreal emission_bincoeff = 0.;
    CFreal B_bincoeff = 0.;
    
    if((test>=vctBins[k-1]) && (test<vctBins[k])) {
      // const CFreal dWav = dWavIn*1e-10;
      // B_bincoeff        = Bs*dWav;
      // alpha_bincoeff    = abCoef*B_bincoeff;
      // emission_bincoeff = emCoef*dWav;
      
      alpha_bincoeff = abCoef*Bs*dWavIn*1e-10;
      emission_bincoeff = emCoef*dWavIn*1e-10;
      B_bincoeff = Bs*dWavIn*1e-10;
    }
    else{
      B_bincoeff = 0.;
      alpha_bincoeff = 0.;
      emission_bincoeff = 0.;
    }
    
    const CFuint idx  = k + idx0;
    alpha_bin[idx]    += alpha_bincoeff;
    emission_bin[idx] += emission_bincoeff;
    B_binCurr[idx]    += B_bincoeff;
  }
}

//////////////////////////////////////////////////////////////////////////////
    
  } // namespace RadiativeTransfer

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_RadiativeTransfer_ParadeRadiator_hh
