#ifndef COOLFluiD_RadiativeTransfer_QSS_MOLECULES_H
#define COOLFluiD_RadiativeTransfer_QSS_MOLECULES_H

#include "RadiativeTransfer/RadiationLibrary/Models/HSNB/core/AtomicPartFunc.h"
#include "RadiativeTransfer/RadiationLibrary/Models/HSNB/core/MolecularPartFunc.h"
#include "RadiativeTransfer/RadiationLibrary/Models/HSNB/core/MolecularElecPartFunc.h"
#include <string>

enum QssModel
{
    JOHNSTON_FIT,
    DISS_EQUIL,
    CN_VIOLET,
    NONE
};

class QssMolecules
{
public:
    QssMolecules(const std::string& name_system, const std::string& name_species, const std::string datadir);

    ~QssMolecules();

    // This function computes a correction to provide to the spontaneous emission coefficient
    double levelCorrection(COOLFluiD::RadiativeTransfer::ThermoData &thermo);

private:

    double levelCorrectCNViolet(COOLFluiD::RadiativeTransfer::ThermoData &thermo);

    // Chemical equilibrium between the molecule and the dissociated atom
    // is assumed to compute the correction
    double levelCorrectionDissEquil(COOLFluiD::RadiativeTransfer::ThermoData &thermo);

    // Johnston's fit are used to compute non-equilibrium population
    double levelCorrectionJohnston(const COOLFluiD::RadiativeTransfer::ThermoData& thermo);

    double fitJohnston(double const Te, double const Ne);

    void readJohnstonParameters();

private:

    QssModel m_qss_model;

    int m_upper_energy_index;
    char m_upper_level_char;
    double* mp_johnston_params;


    std::string m_datadir;
    std::string m_system;
    std::string m_species;
};


#endif // QSS_MOLECULES_H
