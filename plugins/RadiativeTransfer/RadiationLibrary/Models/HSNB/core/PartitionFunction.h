#ifndef COOLFluiD_RadiativeTransfer_PARTITION_FUNCTION_H
#define COOLFluiD_RadiativeTransfer_PARTITION_FUNCTION_H

#include "RadiativeTransfer/RadiationLibrary/Models/HSNB/core/ThermoData.h"

using namespace COOLFluiD::RadiativeTransfer;

/**
 * Abstract base class for all partition function classes.
 */
class PartitionFunction
{
public:

    /// Empty constructor.
    PartitionFunction() {}

    /// Empty destructor.
    virtual ~PartitionFunction() {}

    /**
     * Returns partition function for the associated species and thermodynamic
     * state.
     */
    virtual double Q(ThermoData& thermo) = 0;
};

#endif // PARTITION_FUNCTION_H
