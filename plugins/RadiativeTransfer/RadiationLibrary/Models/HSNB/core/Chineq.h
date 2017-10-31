#ifndef COOLFluiD_RadiativeTransfer_CHINEQ_H
#define COOLFluiD_RadiativeTransfer_CHINEQ_H

#include "RadiativeTransfer/RadiationLibrary/Models/HSNB/core/AtomicPartFunc.h"
#include "RadiativeTransfer/RadiationLibrary/Models/HSNB/core/MolecularPartFunc.h"
//#include "RadiativeTransfer/RadiationLibrary/Models/HSNB/core/FieldData.h"
#include <string>
#include "RadiativeTransfer/RadiationLibrary/Models/HSNB/core/ThermoData.h"

using namespace COOLFluiD::RadiativeTransfer;


class Chineq
{
public:
    Chineq(const std::string& name_system, const std::string& name_species)
        : m_system(name_system), m_species(name_species)
    { }

    double computeChi(ThermoData &thermo);

    double concNegIon(ThermoData& thermo);

    double ionizationEnergy();

private:

    double atomicPhotoionization(ThermoData &thermo);

    double molecularPhotoionization(ThermoData &thermo);

    double oxygenPhotodissociation(ThermoData &thermo);

private:

    std::string m_system;
    std::string m_species;
};


#endif // CHINEQ_H
