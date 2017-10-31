#include <iostream>
#include <iomanip>
#include <fstream>

#include "RadiativeTransfer/RadiationLibrary/Models/HSNB/core/FieldData.h"
#include "RadiativeTransfer/RadiationLibrary/Models/HSNB/core/LineOfSight.h"
#include "RadiativeTransfer/RadiationLibrary/Models/HSNB/core/SnbAtomicSystem.h"
#include "RadiativeTransfer/RadiationLibrary/Models/HSNB/core/SnbDiatomicSystem.h"
#include "RadiativeTransfer/RadiationLibrary/Models/HSNB/core/SnbCO2System.h"
#include "RadiativeTransfer/RadiationLibrary/Models/HSNB/core/SnbContinuumSystem.h"
#include "RadiativeTransfer/RadiationLibrary/Models/HSNB/core/Operators.h"
#include "RadiativeTransfer/RadiationLibrary/Models/HSNB/core/StringUtils.h"

/**
 * @class RteSolution
 * @brief Abstract interface for all RTE solvers.
 */
class RteSolution
{

public:

  RteSolution(const FieldData& field,
              SnbAtomicSystem* atoms, 
              std::vector<SnbDiatomicSystem>& diatomics,
              std::vector<SnbContinuumSystem>& continua,
              SnbCO2System* co2);

  ~RteSolution();

  int spectralGridSize() const { return m_nbands; } 

  double spectralFluxCpuTime() const {
      return (double) m_spec_flux_time / CLOCKS_PER_SEC;
  }

  double spectraCpuTime() const {
      return (double) m_spectra_time / CLOCKS_PER_SEC;
  }

  double eqInt(double nu,double t);

  void computeChemicalSourceTerms(double* const p_chemSource, const double* const p_specUnu);

  void computeEnergySourceTerm(double* const p_energySource, const double* const p_specFlux);

  void computeOmegaSourceTerm(double* const p_omegaSource, const double* const p_specUnu);

  void computeFluxField(double* const p_specFlux, double* const p_specUnu);

  void computeIntensity(bool by_system = false);

  void computeTau();

  void intensityCheck();

private:

  void computePath(
      const LineOfSight& los, double* const p_specInt, bool by_system = false);

  void setupBC(double* const p_wallInt, double twall, bool by_system = false);

  void writeFieldResults(const std::string& field_name, const double* const p_field);

  void writeResults(const std::string& field_name, const double* const p_field, bool by_system = false);

  void determineBandRange();

private:

  FieldData m_field;
  SnbAtomicSystem* m_atoms;
  std::vector<SnbDiatomicSystem> m_diatomics;
  std::vector<SnbContinuumSystem> m_continua;
  SnbCO2System* m_co2;
  
  size_t m_bmin, m_bmax;
  int m_nbands;
  clock_t m_spectra_time;
  clock_t m_spec_flux_time;
};
