#ifndef COOLFluiD_RadiativeTransfer_HYDROGEN_PART_FUNC_H
#define COOLFluiD_RadiativeTransfer_HYDROGEN_PART_FUNC_H

#include "RadiativeTransfer/RadiationLibrary/Models/HSNB/core/ThermoData.h"
#include "RadiativeTransfer/RadiationLibrary/Models/HSNB/core/PartitionFunction.h"

/**
 * Provides the simple hydrogen partition function.
 */
class HydrogenPartFunc : public PartitionFunction
{
public:

    /// Empty constructor.
    HydrogenPartFunc() {}

    /// Empty destructor.
    virtual ~HydrogenPartFunc() {}

    /**
     * Returns partition function for the Hydrogen atom.
     */
    double Q(ThermoData& thermo)
    {
        const double RY_over_KB = 157887.693241;
        const double Tel = thermo.Tel();
        double q = 1.0, i2;

        for (int i = 2; i < 41; ++i) {
            i2 = (double)i*i;
            q += i2*std::exp(RY_over_KB*(1.0/i2-1.0)/Tel);
        }

        return 2.0*q;
    }
};

#endif // HYDROGEN_PART_FUNC_H
