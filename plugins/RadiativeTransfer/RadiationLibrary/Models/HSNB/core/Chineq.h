#ifndef COOLFluiD_RadiativeTransfer_CHINEQ_H
#define COOLFluiD_RadiativeTransfer_CHINEQ_H

#include "RadiativeTransfer/RadiationLibrary/Models/HSNB/core/AtomicPartFunc.h"
#include "RadiativeTransfer/RadiationLibrary/Models/HSNB/core/MolecularPartFunc.h"
#include "RadiativeTransfer/RadiationLibrary/Models/HSNB/core/ThermoData.h"

class Chineq
{
public:
 Chineq(const std::string& name_system,
	const std::string& name_species, std::string datadir)
   : m_system(name_system), m_species(name_species), m_datadir(datadir)
  { }
  
  double computeChi(COOLFluiD::RadiativeTransfer::ThermoData &thermo);
  
  double concNegIon(COOLFluiD::RadiativeTransfer::ThermoData& thermo);
  
  double ionizationEnergy();
  
private:
  
  double atomicPhotoionization(COOLFluiD::RadiativeTransfer::ThermoData &thermo);
  
  double molecularPhotoionization(COOLFluiD::RadiativeTransfer::ThermoData &thermo);
  
  double oxygenPhotodissociation(COOLFluiD::RadiativeTransfer::ThermoData &thermo);
  
private:

    std::string m_system;
    std::string m_species;
    std::string m_datadir;
};


#endif // CHINEQ_H
