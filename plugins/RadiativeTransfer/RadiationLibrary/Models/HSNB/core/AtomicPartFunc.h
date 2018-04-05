#ifndef COOLFluiD_RadiativeTransfer_ATOMIC_PART_FUNC_H
#define COOLFluiD_RadiativeTransfer_ATOMIC_PART_FUNC_H

#include "RadiativeTransfer/RadiationLibrary/Models/HSNB/core/PartitionFunction.h"
#include <string>
#include <vector>

class AtomicPartFunc : public PartitionFunction
{
public:
    AtomicPartFunc(const std::string& name);

    virtual ~AtomicPartFunc() {}

    // Computes the partition function at the given thermodynamic state
    double Q(ThermoData &thermo);

    // Computes the partition function for electro negative ions
    double Qneg(const double& t); 

    // Computes the Debye length at the given thermodynamic state
    double debye(ThermoData &thermo);

private:

    std::string m_name;

    struct Level {
        double E;
        double g;
        int l;
    };

    std::vector<Level> m_levels;

    struct CompareLevels {
        bool operator() (const Level& l1, const Level& l2) const {
            return (l1.l < l2.l);
        }
    };
};


#endif // ATOMIC_PART_FUNC_H
