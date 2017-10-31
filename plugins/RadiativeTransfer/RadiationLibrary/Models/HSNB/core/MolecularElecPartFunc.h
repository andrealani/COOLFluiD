#ifndef COOLFluiD_RadiativeTransfer_MOLECULAR_ELEC_PART_FUNC_H
#define COOLFluiD_RadiativeTransfer_MOLECULAR_ELEC_PART_FUNC_H

//#include "RadiativeTransfer/RadiationLibrary/Models/HSNB/core/FieldData.h"
#include <cstdio>
#include <cmath>
#include <string>

class MolecularElecPartFunc
{
public:

    // Constructor
    MolecularElecPartFunc(const std::string& name);

    // Destructor
    ~MolecularElecPartFunc();

    // Returns the electronic partition function of the i-th electronic
    // level at temperatures (tr,tv)
    double Qel(const int i, const double& tr, const double& tv);

    // Returns the energy of the i-th electronic level
    double Eel(const int i) const { return mp_energy[i]; }

private:

    void loadEnergyLevels();
    void loadTemperatureGrid();
    void loadParameters();

private:

    std::string m_name;
    float*  mp_tv;
    float*  mp_tr;
    double* mp_params;
    double* mp_energy;

    size_t m_ntv;
    size_t m_ntr;
    size_t m_elec;

};

#endif // MOLECULAR_PART_FUNC_H
