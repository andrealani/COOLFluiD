#ifndef COOLFluiD_RadiativeTransfer_MOLECULAR_PART_FUNC_H
#define COOLFluiD_RadiativeTransfer_MOLECULAR_PART_FUNC_H

//#include "RadiativeTransfer/RadiationLibrary/Models/HSNB/core/FieldData.h"
#include <cstdio>
#include <cmath>
#include <string>
#include "RadiativeTransfer/RadiationLibrary/Models/HSNB/core/ThermoData.h"

using namespace COOLFluiD::RadiativeTransfer;

class MolecularPartFunc
{
public:

    // Constructor
    MolecularPartFunc(const std::string& name);

    // Destructor
    ~MolecularPartFunc();

    double Q(const double& tr, const double& tv);

    double debye(ThermoData& thermo);

private:

    void loadTemperatureGrid();
    void loadParameters();

private:

    std::string m_name;
    float*  mp_tv;
    float*  mp_tr;
    double* mp_params;

    const int m_ntv;
    const int m_ntr;

};


#endif // MOLECULAR_PART_FUNC_H
