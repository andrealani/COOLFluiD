#ifndef COOLFluiD_RadiativeTransfer_LORENTZ_HYDROGEN_H
#define COOLFluiD_RadiativeTransfer_LORENTZ_HYDROGEN_H

#include "RadiativeTransfer/RadiationLibrary/Models/HSNB/core/Constants.h"
#include "RadiativeTransfer/RadiationLibrary/Models/HSNB/core/LorentzModel.h"
#include <vector>

#include "RadiativeTransfer/RadiationLibrary/Models/HSNB/core/ThermoData.h"

using namespace COOLFluiD::RadiativeTransfer;


/**
 * Implements a simple Lorentz model for Hydrogen.
 */
class LorentzHydrogen : public LorentzModel
{
public:
    /// Constructor
    LorentzHydrogen(
        const std::string& species,
        const std::vector< std::pair<int, int> >& levels) :
        m_ne23(0.0),
        m_vcs_facs(levels.size(), 8065.73*1.05e-15)
    {
        int n, m;
        for (int i = 0; i < levels.size(); ++i) {
            n = levels[i].second;
            m = levels[i].first;

            m_vcs_facs[i] *= (n == m+1 ? 0.642 : 1.0)*(n*n - m*m);
        }
    }

    /// Destructor
    virtual ~LorentzHydrogen() {}

    /// Updates thermo dependent data.
    void update(ThermoData& thermo) {
        m_ne23 = std::pow(1.0e-6*thermo.N(thermo.speciesIndex("e-")), 2.0/3.0);
    }

    /**
     * Returns the Lorentz HWHM of line i.
     */
    double hwhm(int i) {
        return m_ne23*m_vcs_facs[i];
    }

    void setHWHM(std::vector<LineData>& lines)
    {
        for (int i = 0; i < lines.size(); ++i)
            lines[i].gaml = m_ne23*m_vcs_facs[i];
    }

private:

    double m_ne23;
    std::vector<double> m_vcs_facs;
};

#endif // LORENTZ_HYDROGEN_H
