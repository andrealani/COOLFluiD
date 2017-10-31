#ifndef COOLFluiD_RadiativeTransfer_LORENTZ_MODEL_H
#define COOLFluiD_RadiativeTransfer_LORENTZ_MODEL_H

#include "RadiativeTransfer/RadiationLibrary/Models/HSNB/core/LineData.h"
#include "RadiativeTransfer/RadiationLibrary/Models/HSNB/core/ThermoData.h"
#include <vector>

using namespace COOLFluiD::RadiativeTransfer;

/**
 * Abstract base class for all Lorentz broadening models.
 */
class LorentzModel
{
public:
    /// Constructor
    LorentzModel() {};

    /// Destructor
    virtual ~LorentzModel() {};

    /**
     * Updates the line data vector with the correct Lorentz HWHM.  Default
     * behavior is to set the HWHM values to zero.
     */
    virtual void update(const ThermoData& thermo) { }

    /**
     * Returns the Lorentz HWHM of line i.
     */
    virtual double hwhm(int i) { return 0.0; }

    virtual void setHWHM(std::vector<LineData>& lines) = 0;
};

#endif // LORENTZ_MODEL_H
