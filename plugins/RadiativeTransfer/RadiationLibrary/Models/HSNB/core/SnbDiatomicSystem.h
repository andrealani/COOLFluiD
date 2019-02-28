#ifndef COOLFluiD_RadiativeTransfer_SNB_DIATOMIC_SYSTEM_H
#define COOLFluiD_RadiativeTransfer_SNB_DIATOMIC_SYSTEM_H

#include <string>
#include <cstdio>
#include <cmath>

#include "RadiativeTransfer/RadiationLibrary/Models/HSNB/core/ThermoData.h"
#include "RadiativeTransfer/RadiationLibrary/Models/HSNB/core/RadiativeSystem.h"
#include "RadiativeTransfer/RadiationLibrary/Models/HSNB/core/Operators.h"
#include "RadiativeTransfer/RadiationLibrary/Models/HSNB/core/SpeciesLoadData.h"
#include "RadiativeTransfer/Solvers/MonteCarlo/HSNBLocalParameterSet.hh"
#include "RadiativeTransfer/Solvers/MonteCarlo/HSNBPhotonTrace.hh"

class PhotonPath;

enum DopplerPDF
{
    EXPONENTIAL,
    TAILED_INVERSE_EXP
};
/*
enum NonUniformPath
{
    THIN,
    CURTIS_GODSON,
    LINDQUIST_SIMMONS
};
*/

class SnbDiatomicSystem : public RadiativeSystem<SnbDiatomicSystem>
{
public:

    /**
     * Constructor takes path to system name in the database.
     */
  SnbDiatomicSystem(COOLFluiD::RadiativeTransfer::SpeciesLoadData loadData,
		    const std::string& dpdf, const std::string& nup,
		    const COOLFluiD::RadiativeTransfer::ThermoData &thermo);
  
  /**
   * Copy constructor.
     */
  SnbDiatomicSystem(const SnbDiatomicSystem& system);
  
  /**
   * Destructor.
     */
    ~SnbDiatomicSystem();
    
    /**
     * Assignment operator.
     */
    SnbDiatomicSystem& operator = (SnbDiatomicSystem system) {
        swap(*this, system);
        return *this;
    }
    
    /**
     * Returns the species name associated with this radiative system.
     */
    std::string speciesName() const {
        return m_species;
    }
    
    /**
     * Returns the name of this system.
     */
    std::string systemName() const {
        return m_system;
    }
    
    size_t lowBand() const { return m_band1; }
    size_t highBand() const { return m_bandn; }
    //size_t nBands() const { return m_nbands; }
    int spectralGridSize() const { return m_nbands; }
    size_t nParams() const { return m_nparams; }
    double waveNumber(int i) const { return (m_band1+i)*1000+500; }
    
    /**
     * Initializes the field of local radiative properties
     * For thick systems: ku, etakappa, bd, bl
     * For thin systems: ku, eta
    */
    void setupLocalParameters(COOLFluiD::RadiativeTransfer::ThermoData &thermo);

    /**
     * Returns the local radiative property k of the band j
     * of the cell i
     */
    double getLocalParameter(const int& i, const int& j, const int& k) const;
    
    double wlod(double ku, double bl);
    double wdod(double ku, double bd);
    double wvod(double ku, double wd, double wl);


    double dwlod(double x, double r);
    double dwdod(double x, double r);
    void wquad(const double ku, const double bd, const double bl, double* const p_kbw);

    /// Returns the total emitted power of this mechanism in the cell.
    double emittedPower(int cellID,
			COOLFluiD::RadiativeTransfer::ThermoData &thermo);
    
    /// Computes the emission for each band
    void bandEmission(const int cellID,
		      COOLFluiD::RadiativeTransfer::ThermoData &thermo, double* const p_emis);
    
    /// Computes the transmissivity integral along a photon's path for optically thin systems
    double tau(const COOLFluiD::RadiativeTransfer::HSNBNonThickParameterSet& pathParams);

    /// Computes the transmissivity integral along a photon's path for optically thick systems
    double tau(const COOLFluiD::RadiativeTransfer::HSNBThickParameterSet &pathParams);

    double tauShortened(const COOLFluiD::RadiativeTransfer::HSNBThickParameterSet &pathParams, COOLFluiD::CFreal kappa0, COOLFluiD::CFreal distance0);


    /// Computes the optical thickness along the path including cells [0,ic).
    double opticalThickness(const COOLFluiD::RadiativeTransfer::HSNBNonThickParameterSet& pathParams);

    /// \brief Computes the optical thickness of a given path for an altered start distance
    ///
    /// To allow for computation of a slightly altered path where only the distance between the first two cells in
    /// the photon's path is altered. THe he cust
    double opticalThicknessShortened(const COOLFluiD::RadiativeTransfer::HSNBThickParameterSet &pathParams, COOLFluiD::CFreal kappa0, COOLFluiD::CFreal distance0);

    /// \brief Computes the optical thickness for a given photon path
    /// 
    /// \param pathParams set of parameters necessary
    double opticalThickness(const COOLFluiD::RadiativeTransfer::HSNBThickParameterSet &pathParams);

    ///
    /// \brief Add the parameters necessary for the computation of absorption to a photon path
    ///        If the values for \f$ \kappa and \beta \f$ have been precomputed they are loaded
    ///     from the mp_locparams lookup table. Otherwise they are interpolated in place for the
    ///     correct band.
    ///
    ///
    /// \param thermo struct handling access to the states
    /// \param nonThickParams Struct to store all parameters along a photons path to compute absorption
    /// \param cellDistance Distance that the photon has crossed from last to current cell
    /// \param localCellID id of the current cell
    /// \param sig Wavenumber of the photon's band
    ///
    void addStateParams
      (COOLFluiD::RadiativeTransfer::ThermoData &thermo,
       COOLFluiD::RadiativeTransfer::HSNBNonThickParameterSet& nonThickParams,
       COOLFluiD::CFreal cellDistance,
       COOLFluiD::CFuint localCellID, COOLFluiD::CFreal sig);
    
    ///
    /// \brief Add the parameters necessary for the computation of absorption to a photon path
    ///        If the values for \f$ \kappa and \beta \f$ have been precomputed they are loaded
    ///     from the mp_locparams lookup table. Otherwise they are interpolated in place for the
    ///     correct band. For optically thick systems \kappa and \beta have to be saved for all crossed cells
    ///     both for Doppler and Lorentz profiles
    ///
    ///
    /// \param thermo struct handling access to the states
    /// \param nonThickParams Struct to store all parameters along a photons path to compute absorption
    /// \param cellDistance Distance that the photon has crossed from last to current cell
    /// \param localCellID id of the current cell
    /// \param sig Wavenumber of the photon's band
    ///
    void addStateParams
      (COOLFluiD::RadiativeTransfer::ThermoData &thermo,
       COOLFluiD::RadiativeTransfer::HSNBThickParameterSet& thickParams,
       COOLFluiD::CFreal cellDistance,
       COOLFluiD::CFuint localCellID, COOLFluiD::CFreal sig);
    
    NonUniformPath getMechanismType() const;
    
    COOLFluiD::CFreal getBand(COOLFluiD::CFreal sig);
    
    COOLFluiD::CFreal getKappa
      (COOLFluiD::RadiativeTransfer::ThermoData &thermo,
       COOLFluiD::CFuint cellID, double sig);
    
    /// Returns the memory cost of storing all precomputed
    /// parameters in mp_locparams in byte (this equals the
    /// memory saving to be expected by setting usePrecomputedParamters
    /// to false)
    double getPrecomputedBufferByteSize() const;


    friend void swap(SnbDiatomicSystem&, SnbDiatomicSystem&);
    
private:


    ///
    /// \brief Checks whether the emission for all bands has been computed for the current emitting cell
    ///
    inline bool bandEmissionComputed(int forCellID);

    /**
     * Figures out what the max band range is for this band system.
     */
    void determineBandRange();

    ///
    /// \brief Checks whether a waveNumber lies within the band range availble in the database
    /// \param waveBand waveNumber to check availability for
    /// \return true if the wavenumber is in the available band range for this species.
    ///
    bool waveNumberIsInBandRange(COOLFluiD::CFuint waveBand);

    /**
     * Reads the temperature grid from the first band table.  This assumes that
     * the temperature grid is the same for all bands.
     */
    void loadTemperatureGrid();
    
    /**
     * Loads the parameter information for the given band.
     */
    void loadBandParameters(const size_t& iband);

    /**
     * Uses the given Tv and Tr to interpolate the parameters for each band from
     * the stored table data.
     */
    void getParameters(double Tv, double Tr, double* const p_params);


    ///
    /// \brief Get radiation parameters for optically thick parameters without using a lookup table
    /// Essentially uses the same interpolation routine as in setupLocalParameters only for one a single band
    ///
    /// \param thermo struct handling access to the states
    /// \param cellID id of the current cell
    /// \param band Index of the band to interpolate parameters for (using local indexing)
    ///
    void getThickParamsSingleBand(COOLFluiD::RadiativeTransfer::ThermoData &thermo, COOLFluiD::CFuint cellID, double band);

    ///
    /// \brief Get radiation parameters for optically thin parameters without using a lookup table
    /// Essentially uses the same interpolation routine as in setupLocalParameters only for one a single band
    ///
    /// \param thermo struct handling access to the states
    /// \param cellID id of the current cell
    /// \param band Index of the band to interpolate parameters for (using local indexing)
    ///
    void getThinParamsSingleBand(COOLFluiD::RadiativeTransfer::ThermoData &thermo, COOLFluiD::CFuint cellID, double band);

    std::string bandFilename(size_t band)
    {
        char filename [7];
        sprintf(filename, "%06zu", band*1000+500);
        return std::string(filename);
    } 


 
private:
    

    std::string m_directory;
    std::string m_species;
    std::string m_system;
    
    int m_species_index;
    DopplerPDF m_dpdf;
    NonUniformPath m_nup;

    bool m_test_qss;
    
    int m_band1;
    int m_bandn;
    int m_nbands;
    int m_nparams;
    int m_nv;
    int m_npoints;
    
    float*  mp_tv;
    float*  mp_tr;
    size_t* mp_iv;
    
    double* mp_params;
    double* mp_locparams;

    COOLFluiD::CFreal m_tempKu;
    COOLFluiD::CFreal m_tempBetaD;
    COOLFluiD::CFreal m_tempBetaL;
    COOLFluiD::CFint m_tempBand;

    bool m_usePrecomputedParams;
    bool localParamsSetup=false;
    COOLFluiD::CFuint m_paramCount=0;


    //Local variables to save memory allocations
    COOLFluiD::CFreal m_tv;
    COOLFluiD::CFreal m_tr;
    COOLFluiD::CFreal m_pa;
    COOLFluiD::CFreal mp_fac;

    float* m_pLower;
    size_t m_itv;
    double* m_ParamsCurrentCell;


    size_t m_nr0;
    size_t m_nr1;
    float* m_pMax;
    float* m_pMin;
    size_t m_itr;

    /*-------------------*/
    //Temporary variables for bilinear interpolation
    double m_tv1;
    double m_tv2;
    double m_tr1;
    double m_tr2;

    double m_area;
    double m_w11;
    double m_w12;
    double m_w21;
    double m_w22;

    size_t m_ndata;
    const double*  m_p11;
    const double*  m_p12;
    const double*  m_p21;
    const double*  m_p22;
    /*-------------------*/

    int m_currentParamIndex;

    /// Keeps track of the current cell that the band emission array has already been
    /// computed for to avoid redundant computation of this costly operation
    int m_currentBandEmissionCellIndex;

//    m_currentParamIndex

};

void swap(SnbDiatomicSystem& s1, SnbDiatomicSystem& s2);

// Quadritures
double equad(double x);
double hquad(double x);

#endif // SNB_DIATOMIC_SYSTEM_H

