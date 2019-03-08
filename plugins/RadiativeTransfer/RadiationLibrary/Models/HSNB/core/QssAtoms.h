#ifndef COOLFluiD_RadiativeTransfer_QSS_ATOMS_H
#define COOLFluiD_RadiativeTransfer_QSS_ATOMS_H

#include "RadiativeTransfer/RadiationLibrary/Models/HSNB/core/AtomicPartFunc.h"
#include <string>

class QssAtoms
{
public:
    QssAtoms(const std::string& name_species, const std::string datadir);

    ~QssAtoms();

    // Computes corrections for each group of Johnston's QSS model
    void computeCorrections(COOLFluiD::RadiativeTransfer::ThermoData &thermo);

    // Returns the correction to apply to the population of each 
    // electronic state belonging to the group g
    double getCorrection(int const g) {return mp_fac[g];}

private:

    void fitJohnston(double const Te, double const Ne,
                     double* const p_r0, double* const p_r1);

    void readJohnstonParameters();

private:

    size_t m_ngroups;

    double* mp_r0_params;
    double* mp_r1_params;
    double* mp_fac;

    std::string m_species;
    std::string m_datadir;
};

#endif // QSS_ATOMS_H
