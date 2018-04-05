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

using namespace COOLFluiD;
using namespace COOLFluiD::RadiativeTransfer;

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
    SnbDiatomicSystem(SpeciesLoadData loadData, const std::string& dpdf, const std::string& nup, const ThermoData &thermo);
    
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
    void setupLocalParameters(ThermoData &thermo);

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
    double emittedPower(int cellID, ThermoData &thermo);

    /// Computes the emission for each band
    void bandEmission(const int cellID, ThermoData &thermo, double* const p_emis);



    double tau(const HSNBNonThickParameterSet& pathParams);
    double tau(const HSNBThickParameterSet &pathParams);

    double tauShortened(const HSNBThickParameterSet &pathParams, CFreal kappa0, CFreal distance0);


    /// Computes the optical thickness along the path including cells [0,ic).
    double opticalThickness(const HSNBNonThickParameterSet& pathParams);

    double opticalThicknessShortened(const HSNBThickParameterSet &pathParams, CFreal kappa0, CFreal distance0);

    double opticalThickness(const HSNBThickParameterSet &pathParams);

    void addStateParams(ThermoData &thermo, HSNBNonThickParameterSet& nonThickParams, CFreal cellDistance, CFuint localCellID, CFreal sig);

    void addStateParams(HSNBThickParameterSet& thickParams, CFreal cellDistance, CFuint localCellID, CFreal sig);


    NonUniformPath getMechanismType() const;

    CFreal getBand(CFreal sig);



    friend void swap(SnbDiatomicSystem&, SnbDiatomicSystem&);
    
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
    void getParameters(double Tv, double Tr, double* const p_params) const;

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

    CFreal m_tempKu;
    CFreal m_tempBetaD;
    CFreal m_tempBetaL;
    CFint m_tempBand;

    bool localParamsSetup=false;
    CFuint paramCount=0;
};

void swap(SnbDiatomicSystem& s1, SnbDiatomicSystem& s2);

// Quadritures
double equad(double x);
double hquad(double x);

#endif // SNB_DIATOMIC_SYSTEM_H

