#ifndef COOLFluiD_RadiativeTransfer_RADIATIVE_SYSTEM_H
#define COOLFluiD_RadiativeTransfer_RADIATIVE_SYSTEM_H

#include <string>

#include "RadiativeTransfer/RadiationLibrary/Models/HSNB/core/Operators.h"

enum NonUniformPath
{
    THIN,
    CURTIS_GODSON,
    LINDQUIST_SIMMONS
};
/**
 * @class RadiativeSystem
 * @brief Provides static inheritance interface for all radiative systems.
 */
template <typename SysType>
class RadiativeSystem
{
public:
    
    /**
     * Returns the species name associated with this radiative system.
     */
    std::string speciesName() const {
        return static_cast<const SysType&>(*this).speciesName();
    }
    
    /**
     * Returns the name of this system.
     */
    std::string systemName() const {
        return static_cast<const SysType&>(*this).systemName();
    }
    
    /**
     * Returns the size of the spectral grid.
     */
    int spectralGridSize() const {
        return static_cast<const SysType&>(*this).spectralGridSize();
    }
    
    /**
     * Returns the wavenumber associated with index i in the spectral grid.
     */
    double waveNumber(int i) const {
        return static_cast<const SysType&>(*this).waveNumber(i);
    }


private:

    std::string bandFilename(size_t band) {
        return static_cast<const SysType&>(*this).bandFilemame(band);
    }
};

#endif
