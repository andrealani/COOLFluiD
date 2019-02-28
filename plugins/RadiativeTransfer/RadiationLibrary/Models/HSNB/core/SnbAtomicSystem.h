#ifndef COOLFluiD_RadiativeTransfer_SNB_ATOMIC_SYSTEM_H
#define COOLFluiD_RadiativeTransfer_SNB_ATOMIC_SYSTEM_H

#include <string>
#include <cstdio>
#include <cmath>
#include <ctime>

#include "RadiativeTransfer/RadiationLibrary/Models/HSNB/core/RadiativeSystem.h"
#include "RadiativeTransfer/RadiationLibrary/Models/HSNB/core/Operators.h"
#include "RadiativeTransfer/RadiationLibrary/Models/HSNB/core/LblSpectralGrid.h"
#include "RadiativeTransfer/RadiationLibrary/Models/HSNB/core/LblAtomicSystem.h"
#include "RadiativeTransfer/RadiationLibrary/Models/HSNB/core/ThermoData.h"

class SnbAtomicSystem : public RadiativeSystem<SnbAtomicSystem>
{
public:

    /**
     * Constructor gets the name of the systems to consider.
     */
    SnbAtomicSystem(const std::string& spectradir, 
                    const std::string& computespectra,
                    const std::string& grid_type,
                    const std::vector<std::string>& systemlist, COOLFluiD::RadiativeTransfer::ThermoData &thermo);
    
    /**
     * Copy constructor.
     */
    SnbAtomicSystem(const SnbAtomicSystem& sys)
        : m_systems(sys.m_systems),
          m_specdir(sys.m_specdir),
          mp_grid(
              sys.mp_grid == NULL ? NULL : new LblSpectralGrid(*(sys.mp_grid))),
          mp_bounds(sys.mp_bounds == NULL ? NULL : new int [sys.m_nbands+1]),
          m_nloclblparams(sys.m_nloclblparams),
          mp_loclblparams(
              sys.mp_loclblparams == NULL ? NULL :
              new double [sys.m_nloclblparams]),
          m_nlocparams(sys.m_nlocparams),
          mp_locparams(
              sys.mp_locparams == NULL ? NULL : new double [sys.m_nlocparams]),
          m_species_index(sys.m_species_index),
          m_nsystems(sys.m_nsystems),
          m_band1(sys.m_band1),
          m_bandn(sys.m_bandn),
          m_nbands(sys.m_nbands),
          m_computelbl(sys.m_computelbl),
          m_spectra_time(sys.m_spectra_time)
    {
        std::copy(sys.mp_bounds, sys.mp_bounds+m_nbands+1, mp_bounds);
        std::copy(
            sys.mp_loclblparams, sys.mp_loclblparams+m_nloclblparams,
            mp_loclblparams);
        std::copy(
            sys.mp_locparams, sys.mp_locparams+m_nlocparams,
            mp_locparams);
    }

    /**
     * Destructor.
     */
    ~SnbAtomicSystem();
    
    /**
     * Assignment operator.
     */
    SnbAtomicSystem& operator= (SnbAtomicSystem sys)
    {
        swap(*this, sys);
        return *this;
    }

    clock_t getSpectraCpuClocks() const { return m_spectra_time; }

    size_t lowBand() const { return m_band1; }
    size_t highBand() const { return m_bandn; }
    int spectralGridSize() const { return m_nbands; }
    double waveNumber(int i) const { return (m_band1+i)*1000+500; }

    void setupLocalParameters(COOLFluiD::RadiativeTransfer::ThermoData &thermo);

    double getLocalParameter(const int& i, const int& j, const int& k) const;


    friend void swap(SnbAtomicSystem&, SnbAtomicSystem&);
    
private:
    
    /**
     * Figures out what the band range is for this band system.
     */
    void determineBandRange();
    void readLblSpectra(int index, double* const p_spectra);
    void writeLblSpectra(int index, const double* const p_spectra);
    void computeLblSpectra(int index, COOLFluiD::RadiativeTransfer::ThermoData &thermo,
             double* const p_spectra);
    void computeLblSpectraEq(int index, COOLFluiD::RadiativeTransfer::ThermoData &thermo,
             double* const p_spectra);
    void lblToSnb(double* const p_lbl, double* const p_snb);
    void lblToSnbInit();

 
private:
    
    std::vector<LblAtomicSystem> m_systems;
    std::string m_specdir;

    LblSpectralGrid* mp_grid;
    int* mp_bounds;

    size_t  m_nloclblparams;
    double* mp_loclblparams;
    size_t  m_nlocparams;
    double* mp_locparams;

    std::vector<int> m_species_index;
    size_t m_nsystems;
    
    size_t m_band1;
    size_t m_bandn;
    size_t m_nbands;

    bool m_computelbl;
    clock_t m_spectra_time;
};

#endif // SNB_ATOMIC_SYSTEM_H
