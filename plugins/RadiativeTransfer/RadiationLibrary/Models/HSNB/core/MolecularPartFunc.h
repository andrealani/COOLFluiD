#ifndef COOLFluiD_RadiativeTransfer_MOLECULAR_PART_FUNC_H
#define COOLFluiD_RadiativeTransfer_MOLECULAR_PART_FUNC_H

//#include "RadiativeTransfer/RadiationLibrary/Models/HSNB/core/FieldData.h"
#include <cstdio>
#include <cmath>
#include <string>
#include "RadiativeTransfer/RadiationLibrary/Models/HSNB/core/ThermoData.h"

class MolecularPartFunc
{
public:

    // Constructor
    MolecularPartFunc(const std::string& name, const std::string datadir);

    // Destructor
    ~MolecularPartFunc();

    double Q(const double& tr, const double& tv);

    double debye(COOLFluiD::RadiativeTransfer::ThermoData& thermo);

private:

    void loadTemperatureGrid();
    void loadParameters();

private:

    std::string m_name;
    std::string m_datadir;

    float*  mp_tv;
    float*  mp_tr;
    double* mp_params;

    const int m_ntv;
    const int m_ntr;

};


#endif // MOLECULAR_PART_FUNC_H
