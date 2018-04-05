#ifndef COOLFluiD_RadiativeTransfer_SNB_CO2_SYSTEM_H
#define COOLFluiD_RadiativeTransfer_SNB_CO2_SYSTEM_H

#include <string>
#include <cstdio>
#include <cmath>
#include "RadiativeTransfer/RadiationLibrary/Models/HSNB/core/ThermoData.h"
#include "RadiativeTransfer/RadiationLibrary/Models/HSNB/core/Constants.h"
#include "RadiativeTransfer/RadiationLibrary/Models/HSNB/core/ThermoData.h"
#include "RadiativeTransfer/RadiationLibrary/Models/HSNB/core/RadiativeSystem.h"
#include "RadiativeTransfer/RadiationLibrary/Models/HSNB/core/Operators.h"
/*
enum NonUniformPath
{
    CURTIS_GODSON,
    LINDQUIST_SIMMONS
};
*/

using namespace COOLFluiD;
using namespace COOLFluiD::RadiativeTransfer;

class SnbCO2System : public RadiativeSystem<SnbCO2System>
{
public:

    /**
     * Constructor takes path to system name in the database.
     */
    SnbCO2System(const std::string& system, const std::string& nup, const ThermoData &thermo);
    
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
    
    size_t lowBand() const { return m_band1_up; }
    size_t highBand() const { return m_bandn_up; }
    int spectralGridSize() const { return m_nbands_up; }
    double waveNumber(int i) const { return (m_band1_up+i)*1000+500; }
    
    size_t nParams() const { return m_nparams; }

    /**
     * Initializes the field of local radiative properties
    */
    void setupLocalParameters(ThermoData &thermo);

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

    void downToUpSnb(const double* const p_down, double* const p_up);
 
private:
    
    std::string m_directory;
    
    int m_species_index;
    NonUniformPath m_nup;
    bool m_lorentz;
    
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

