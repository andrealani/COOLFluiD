#ifndef COOLFluiD_RadiativeTransfer_SNB_CO2_SYSTEM_H
#define COOLFluiD_RadiativeTransfer_SNB_CO2_SYSTEM_H

#include "RadiativeTransfer/RadiationLibrary/Models/HSNB/core/ThermoData.h"
#include "RadiativeTransfer/RadiationLibrary/Models/HSNB/core/Constants.h"
#include "RadiativeTransfer/RadiationLibrary/Models/HSNB/core/ThermoData.h"
#include "RadiativeTransfer/RadiationLibrary/Models/HSNB/core/RadiativeSystem.h"
#include "RadiativeTransfer/RadiationLibrary/Models/HSNB/core/Operators.h"
#include "Framework/DataHandle.hh"

#include "RadiativeTransfer/Solvers/MonteCarlo/HSNBLocalParameterSet.hh"

/*
enum NonUniformPath
{
    CURTIS_GODSON,
    LINDQUIST_SIMMONS
};
*/

const COOLFluiD::CFuint nbLocalParameters = 5;

class SnbCO2System : public RadiativeSystem<SnbCO2System>
{
public:

    /**
     * Constructor takes path to system name in the database.
     */
    SnbCO2System(const std::string& system, const std::string& nup, const COOLFluiD::RadiativeTransfer::ThermoData &thermo);
    
    bool waveNumberIsInBandRange(COOLFluiD::CFuint waveNumber);

    /**
     * Copy constructor.
     */
    SnbCO2System(const SnbCO2System& system);
    
    /**
     * Destructor.
     */
    ~SnbCO2System();
    
    /**
     * Assignment operator.
     */
    SnbCO2System& operator = (SnbCO2System system) {
        swap(*this, system);
        return *this;
    }
    
    size_t lowBand() const { return m_band1_down; }
    size_t highBand() const { return m_bandn_down; }
    int spectralGridSize() const { return m_nbands_down; }
    double waveNumber(int i) const { return (m_band1_down+i)*25.0+12.5; }

    void addStateParams(COOLFluiD::RadiativeTransfer::ThermoData &thermo, COOLFluiD::RadiativeTransfer::HSNBCO2ParameterSet& co2Params, COOLFluiD::CFreal cellDistance, COOLFluiD::CFuint localCellID, COOLFluiD::CFreal sig);
    
    /// Computes the absorption coefficient for a single wavenumber by interpolation
    void getAbsorptionSingleBand(COOLFluiD::RadiativeTransfer::ThermoData &thermo, COOLFluiD::CFuint cellID, double band, double& KD,double& KL, double& betaD, double& betaL);

    /// Computes the emission coefficient for a single wavenummber by interpolation
    ///
    /// Note that the interpolation is set up on every call, using setupParamInterpolation and
    /// getEmissionSingleBandFixedCell is much more efficient if the state remains unchanged between
    /// function calls.
    double getEmissionSingleBand(COOLFluiD::RadiativeTransfer::ThermoData& thermo,COOLFluiD::CFuint cellID, double band);

    /// Access parameter array to compute emission coefficient for a given state
    /// interpolation parameters have to be preset (weights w1, w2 e.g.)
    double getEmissionSingleBandFixedCell(COOLFluiD::CFuint band);

    /// Prepares interpolatio of parameters necessary for efficient in-place computations
    double setupParamInterpolation(COOLFluiD::RadiativeTransfer::ThermoData &thermo, COOLFluiD::CFuint cellID);

    size_t nParams() const { return m_nparams; }

    /**
     * Initializes the field of local radiative properties
    */
    void setupLocalParameters(COOLFluiD::RadiativeTransfer::ThermoData &thermo);

    /**
     * Returns the local radiative property k of the band j
     * of the cell i
     */
    double getLocalParameter(const int& i, const int& j, const int& k) const;
    
    double wlod(double kl, double bl);
    double wdod(double kd, double bd);
    double wvod(double kd, double kl, double wd, double wl);

    double dwlod(double x, double r);
    double dwdod(double x, double r);
    void wquad(const double kd, const double kl,
                 const double bd, const double bl, double* const p_kbw);

    /// Determines whether to use lorentz or doppler parameters.
    void setupProfileType(const std::vector<COOLFluiD::CFreal> &cellVolumes, COOLFluiD::RadiativeTransfer::ThermoData &thermo);
    
    /// Returns the total emitted power of this mechanism in the cell.
    double emittedPower(int cellID, COOLFluiD::RadiativeTransfer::ThermoData &thermo);

    /// Computes the emission for each band
    void bandEmission(const int cellID, COOLFluiD::RadiativeTransfer::ThermoData &thermo, double* const p_emis);

    /// Returns the optical thickness as precomputed in a photon path
    double opticalThickness(const COOLFluiD::RadiativeTransfer::HSNBCO2ParameterSet &pathParams);

    double tau(const COOLFluiD::RadiativeTransfer::HSNBCO2ParameterSet &pathParams);

    //Todo: Implementiere für eine Wavenumber mit bands_down ( im 25cm spectrum ), schreibe eigene emissionsroutine


    friend void swap(SnbCO2System&, SnbCO2System&);
    
private:
    
    /**
     * Figures out what the band range is for this band system.
     */
    void determineBandRange();

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
    void getParameters(double T, double* const p_params) const;

    std::string bandFilename(size_t band)
    {
        char filename [7];
        sprintf(filename, "%06zu", band*25);
        return std::string(filename);
    } 


 
private:
    size_t m_it;
    double m_eta;
    double m_t1;
    double m_t2;
    double m_w1;
    double m_w2;
    double* m_p1;
    double* m_p2;
    size_t m_nParamdata;

    std::string m_directory;
    
    int m_species_index;
    NonUniformPath m_nup;
    bool m_lorentz;

    COOLFluiD::CFreal m_tempKD;
    COOLFluiD::CFreal m_tempKL;

    COOLFluiD::CFreal m_tempBetaL;
    COOLFluiD::CFreal m_tempBetaD;

    COOLFluiD::CFreal m_tempKappa;


    COOLFluiD::CFuint m_ndata;
    COOLFluiD::CFuint m_nbLocalParams;
    
    size_t m_band1_up;
    size_t m_bandn_up;
    size_t m_nbands_up;
    size_t m_band1_down;
    size_t m_bandn_down;
    size_t m_nbands_down;
    size_t m_nparams;
    size_t m_npoints;
    
    float*  mp_t;
    double* mp_params;
    double* mp_locparams_up;
    double* mp_locparams_down;
};

void swap(SnbCO2System& s1, SnbCO2System& s2);

#endif // SNB_CO2_SYSTEM_H

