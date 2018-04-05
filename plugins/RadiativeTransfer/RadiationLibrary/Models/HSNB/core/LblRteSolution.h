#include <iostream>
#include <iomanip>
#include <fstream>


#include "RadiativeTransfer/RadiationLibrary/Models/HSNB/core/LblSpectralGrid.h"
#include "RadiativeTransfer/RadiationLibrary/Models/HSNB/core/FieldData.h"
#include "RadiativeTransfer/RadiationLibrary/Models/HSNB/core/Operators.h"
#include "RadiativeTransfer/RadiationLibrary/Models/HSNB/core/StringUtils.h"

/**
 * @class LblRteSolution
 * @brief Abstract interface for all RTE solvers.
 */
class LblRteSolution
{

public:

  LblRteSolution(const FieldData& field,
                 const LblSpectralGrid& lblgrid,
                 double* const spectrum);

  ~LblRteSolution();

  double eqInt(double nu,double t);

/*
  void computeChemicalSourceTerms(double* const p_chemSource, const double* const p_specUnu);

  void computeOmegaSourceTerm(double* const p_omegaSource, const double* const p_specUnu);
*/

  void computeEnergySourceTerm(double* const p_energySource, const double* const p_specFlux);

  void computeFluxField(double* const p_specFlux);

private:

  double E1(double x);

  double E2(double x);

  double E3(double x);

  void setupBC(double* const p_wallInt, double twall);

  void writeFieldResults(const std::string& field_name, const double* const p_field);

  void writeResults(const std::string& field_name, const double* const p_field);


private:

  FieldData m_field;
  LblSpectralGrid m_grid;
  double* mp_spectrum;
  
};
