#ifndef COOLFluiD_RadiativeTransfer_SNB_CONTINUUM_SYSTEM_H
#define COOLFluiD_RadiativeTransfer_SNB_CONTINUUM_SYSTEM_H

#include "Environment/FileHandlerInput.hh"
#include "Environment/FileHandlerOutput.hh"

#include "RadiativeTransfer/RadiationLibrary/Models/HSNB/core/RadiativeSystem.h"
#include "RadiativeTransfer/RadiationLibrary/Models/HSNB/core/AtomicPartFunc.h"
#include "RadiativeTransfer/RadiationLibrary/Models/HSNB/core/Chineq.h"
#include "RadiativeTransfer/RadiationLibrary/Models/HSNB/core/Operators.h"
#include "RadiativeTransfer/RadiationLibrary/Models/HSNB/core/ThermoData.h"
#include "RadiativeTransfer/RadiationLibrary/Models/HSNB/core/SpeciesLoadData.h"
#include "RadiativeTransfer/Solvers/MonteCarlo/HSNBLocalParameterSet.hh"

class PhotonPath;

enum ContinuumSort
{
    BOUNDFREE,
    FREEFREE
};

class SnbContinuumSystem : public RadiativeSystem<SnbContinuumSystem>
{
public:

    /**
     * Constructor takes path to system name in the database.
     */
    SnbContinuumSystem(COOLFluiD::RadiativeTransfer::SpeciesLoadData loadData,
		       const COOLFluiD::RadiativeTransfer::ThermoData &thermo,
		       ContinuumSort sort = BOUNDFREE);
    
    /**
     * Copy constructor.
     */
    SnbContinuumSystem(const SnbContinuumSystem& system);
    
    /**
     * Destructor.
     */
    ~SnbContinuumSystem();
    
    /**
     * Assignment operator.
     */
    SnbContinuumSystem& operator = (SnbContinuumSystem system) {
        swap(*this, system);
        return *this;
    }
    
    /**
     * Returns the species name associated with this radiative system.
     */
    std::string speciesName() const {
        return m_species;
    }
    
    std::string productsName() const {
        return m_products;
    }

   /**
     * Returns the name of this system.
     */
    std::string systemName() const {
        return m_system;
    }

    Chineq* chi() const { return mp_chi; }    
    size_t lowBand() const { return m_band1; }
    size_t highBand() const { return m_bandn; }
    //size_t nBands() const { return m_nbands; }
    int spectralGridSize() const { return m_nbands; }
    size_t nParams() const { return m_nparams; }
    double waveNumber(int i) const { return (m_band1+i)*1000+500; }
    
    /**
     * Initializes the field of local radiative properties (ku, eta)
    */
    void setupLocalParameters(COOLFluiD::RadiativeTransfer::ThermoData &thermo);

    /**
     * Returns the local radiative property k of the band j
     * of the cell i
    */
    double getLocalParameter(const int& i, const int& j, const int& k) const;

    /**
     * @brief Exports the precomputed parameters into the file specified in exportPath
     *
     * This is useful for continua since the computation of parameters is costly.
     * @param exportPath export directory
     * @param thermo
     */
    void exportLocalParameters(boost::filesystem::path exportPath, COOLFluiD::RadiativeTransfer::ThermoData &thermo);

    void readLocalParameters(boost::filesystem::path importPath);




    /// Returns the total emitted power of this mechanism in the cell.
    double emittedPower(COOLFluiD::RadiativeTransfer::ThermoData &thermo, int i);

    double emittedPowerSaveMemory(int i);

    double localParamsInPlace(double& localEmissivity, int cellID, int band, COOLFluiD::RadiativeTransfer::ThermoData &thermo);

    /// Computes the emission for each band
    void bandEmission(COOLFluiD::RadiativeTransfer::ThermoData &thermo, int i, double* const p_emis);

    double opticalThickness(const PhotonPath& path, int ic, double sig);

    double opticalThickness
      (const COOLFluiD::RadiativeTransfer::HSNBNonThickParameterSet& pathParams);
    
    void addStateParams
      (COOLFluiD::RadiativeTransfer::ThermoData& thermo,
       COOLFluiD::RadiativeTransfer::HSNBNonThickParameterSet &nonThickParams,
       COOLFluiD::CFreal cellDistance, COOLFluiD::CFuint localCellID,
       COOLFluiD::CFreal sig);
    
    friend void swap(SnbContinuumSystem&, SnbContinuumSystem&);
    
private:

    void setupProducts();
    
    /**
     * Figures out what the band range is for this band system.
     */
    void determineBandRange();

    bool waveNumberIsInBandRange(COOLFluiD::CFuint waveNumber);

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
     * Uses the given T to interpolate the parameters for each band from
     * the stored table data.
     */
    void getParameters(double T, double* const p_params);


    void getParametersSingleBand(double T, COOLFluiD::CFuint band, double *p_params);


    double getThinAbsorptionSingleBand(COOLFluiD::RadiativeTransfer::ThermoData &thermo, COOLFluiD::CFuint cellID, double band);
    double getThinEmissionSingleBand(COOLFluiD::RadiativeTransfer::ThermoData &thermo, COOLFluiD::CFuint cellID, double band);

    /// Access parameter array to compute emission coefficient for a given state
    /// interpolation parameters have to be preset (weights w1, w2 e.g.)
    double getEmissionSingleBandFixedCell(COOLFluiD::RadiativeTransfer::ThermoData &thermo, COOLFluiD::CFuint band);

    void setupEmissionCoefficientsFixedCell(COOLFluiD::RadiativeTransfer::ThermoData &thermo, COOLFluiD::CFuint cellID);

    double setupParamInterpolation(COOLFluiD::RadiativeTransfer::ThermoData &thermo, COOLFluiD::CFuint cellID);

    double getParameterAt(double T, size_t band, size_t point) const;

       std::string bandFilename(size_t band)
    {
        char filename [7];
        sprintf(filename, "%06zu", band*1000+500);
        return std::string(filename);
    }
    
private:


       size_t m_it;
       double m_eta;
       double m_t1;
       double m_t2;
       double m_w1;
       double m_w2;
       const double* m_p1;
       const double* m_p2;
       size_t m_ndata;

       double *m_pBandParams;

       double m_nbParams;

       //Emission setup for a fixed cell
       double m_eion;
       double m_q;
       double m_lowT_fac;
       double m_pa;
       double m_pe;
       double m_chi_factor;

       // Electron species index
       size_t m_ie;

       // Local Paramset to be interpolated for the current cell
       double * m_pParamsCurrentCell;

       // Temporary storage to avoid reallocation of convenience variables
       double m_tr;
       double m_tv;

       // Number density array in current cell
       const double* m_nd;


       size_t m_nParamdata;

    COOLFluiD::Common::SelfRegistPtr<COOLFluiD::Environment::FileHandlerOutput> m_outFileHandle;
    COOLFluiD::Common::SelfRegistPtr<COOLFluiD::Environment::FileHandlerInput> m_inFileHandle;

    std::string m_directory;
    std::string m_datadir;
    std::string m_species;
    std::string m_products;
    std::string m_system;
    
    int m_species_index;
    ContinuumSort m_sort;
    Chineq* mp_chi;
    
    size_t m_band1;
    size_t m_bandn;
    size_t m_nbands;
    size_t m_nparams;
    size_t m_npoints;

    float* m_pLower;
    float*  mp_t;
    double* mp_params;
    double* mp_locparams;

    bool m_usePrecomputedParams;

    COOLFluiD::CFreal m_tempKu;
};

void swap(SnbContinuumSystem& s1, SnbContinuumSystem& s2);
#endif // SNB_CONTINUUM_SYSTEM_H

