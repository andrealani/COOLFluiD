#ifndef COOLFluiD_RadiativeTransfer_SNB_CONTINUUM_SYSTEM_H
#define COOLFluiD_RadiativeTransfer_SNB_CONTINUUM_SYSTEM_H

#include <string>
#include <cstdio>
#include <cmath>
#include <iostream>
#include <fstream>

#include "Environment/FileHandlerOutput.hh"
#include "Environment/FileHandlerInput.hh"

#include "Environment/SingleBehaviorFactory.hh"

#include "RadiativeTransfer/RadiationLibrary/Models/HSNB/core/RadiativeSystem.h"
#include "RadiativeTransfer/RadiationLibrary/Models/HSNB/core/AtomicPartFunc.h"
#include "RadiativeTransfer/RadiationLibrary/Models/HSNB/core/Chineq.h"
#include "RadiativeTransfer/RadiationLibrary/Models/HSNB/core/Operators.h"
#include "RadiativeTransfer/RadiationLibrary/Models/HSNB/core/ThermoData.h"
#include "RadiativeTransfer/RadiationLibrary/Models/HSNB/core/SpeciesLoadData.h"

#include "RadiativeTransfer/Solvers/MonteCarlo/HSNBLocalParameterSet.hh"

using namespace COOLFluiD::RadiativeTransfer;


class PhotonPath;

using namespace COOLFluiD;

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
    SnbContinuumSystem(SpeciesLoadData loadData, const ThermoData &thermo, ContinuumSort sort = BOUNDFREE);
    
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
    void setupLocalParameters(ThermoData &thermo);

    /**
     * Returns the local radiative property k of the band j
     * of the cell i
    */
    double getLocalParameter(const int& i, const int& j, const int& k) const;

    void exportLocalParameters(boost::filesystem::path exportPath, ThermoData &thermo);

    void readLocalParameters(boost::filesystem::path importPath);


    /// Returns the total emitted power of this mechanism in the cell.
    double emittedPower(int i);

    double emittedPowerSaveMemory(int i);

    double localParamsInPlace(double& localEmissivity, int cellID, int band, ThermoData &thermo);

    /// Computes the emission for each band
    void bandEmission(int i, double* const p_emis);




    double opticalThickness(const PhotonPath& path, int ic, double sig);

    double opticalThickness(const HSNBNonThickParameterSet& pathParams);

    void addStateParams(HSNBNonThickParameterSet &nonThickParams, CFreal cellDistance, CFuint localCellID, CFreal sig);


    
    friend void swap(SnbContinuumSystem&, SnbContinuumSystem&);
    
private:

    void setupProducts();
    
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
     * Uses the given T to interpolate the parameters for each band from
     * the stored table data.
     */
    void getParameters(double T, double* const p_params) const;

    double getParameterAt(double T, size_t band, size_t point) const;

       std::string bandFilename(size_t band)
    {
        char filename [7];
        sprintf(filename, "%06zu", band*1000+500);
        return std::string(filename);
    }
    
private:

    COOLFluiD::Common::SelfRegistPtr<COOLFluiD::Environment::FileHandlerOutput> m_outFileHandle;
    COOLFluiD::Common::SelfRegistPtr<COOLFluiD::Environment::FileHandlerInput> m_inFileHandle;

    std::string m_directory;
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
    
    float*  mp_t;
    double* mp_params;
    double* mp_locparams;

    CFreal m_tempKu;
};

void swap(SnbContinuumSystem& s1, SnbContinuumSystem& s2);
#endif // SNB_CONTINUUM_SYSTEM_H

