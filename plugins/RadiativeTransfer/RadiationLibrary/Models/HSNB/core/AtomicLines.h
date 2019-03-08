#ifndef COOLFluiD_RadiativeTransfer_ATOMIC_LINES_H
#define COOLFluiD_RadiativeTransfer_ATOMIC_LINES_H

#include "RadiativeTransfer/RadiationLibrary/Models/HSNB/core/LineData.h"
#include "RadiativeTransfer/RadiationLibrary/Models/HSNB/core/LorentzModel.h"
#include "RadiativeTransfer/RadiationLibrary/Models/HSNB/core/PartitionFunction.h"
#include "RadiativeTransfer/RadiationLibrary/Models/HSNB/core/HSNBSharedPtr.h"
#include "RadiativeTransfer/RadiationLibrary/Models/HSNB/core/ThermoData.h"

#include "RadiativeTransfer/Solvers/MonteCarlo/HSNBPhotonTrace.hh"
#include <string>
#include <vector>

class LblSpectralGrid;
class PhotonPath;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

namespace RadiativeTransfer {

//////////////////////////////////////////////////////////////////////////////

/**
 * Represents a collection of atomic lines, and provides functions for computing
 * spectral coefficients.
 */
class AtomicLines
{
public:


    /// Constructor
    AtomicLines(const std::string baseDir, const std::string& atom, const ThermoData &thermo);

    /// Returns number of lines considered for the atom.
    int nLines() const {
        return m_line_data.size();
    }

    /// Returns the name of the atom.
    const std::string& name() const {
        return m_atom;
    }

    /// Returns the current line data vector (from last update)
    const std::vector<LineData>& lineData() const {
        return m_line_data;
    }

    /// Returns the line data for all transitions considered for the atom.
    const std::vector<LineData>& lineData(const ThermoData& thermo);

    /// Computes spectrum based on spectrum type.
    void spectrum(ThermoData &thermo, const LblSpectralGrid& grid,
        double* const p_spectrum);

    bool wavenumberIsValid(double sig);

    /// Computes the optically thin emission of each line.
    void thinLineEmission(ThermoData &thermo, double * const p_emis, CFreal &sumEmission, bool printDebug=false);

    double opticalThickness(ThermoData &thermo, const PhotonPath& path, int ic, double sig);

    double opticalThickness(HSNBAtomicParameterSet &atomParams);

    void addStateParams(ThermoData &thermo, HSNBAtomicParameterSet& atomParams, CFreal cellDistance, CFuint localCellID, CFreal sig);

    /// Initializes the field of local radiative properties.
    void setupLocalParameters(ThermoData &thermo);

    /// Saves out the line data
    void save(const std::string& file_name) const;

private:

    /// Updates the Doppler broadening parameters.
    void updateDopplerHWHM(double T);

    /// Loads line data from a NIST data file with TOPBASE indices.
    void loadNistData();

    /// Loads line data from a special file for H lines.
    void loadHData();

private:
    std::string m_datadir;

    std::string m_atom;   /// Name of atom
    int m_species_index;  /// Species index in thermodynamic database
    double m_mw;          /// Species molecular weight
    bool m_qss;           /// Whether or not to use a qss model

    /// Line data for all transitions considered for this atom
    std::vector<LineData> m_line_data;

    /// Lorentz broadening parameters
    HSNBSharedPtr<LorentzModel> mp_lorentz;

    /// Atomic partition function
    HSNBSharedPtr<PartitionFunction> mp_Q;

    /// Costly parameters, it is useful to store these once.
    struct LocalParams {
        double Q;
        double N;
    };
    std::vector<LocalParams> m_loc_params;

    // Parameters used in opticalThickness
    double m_opt_thick;
    int    m_last_loc;
    double m_last_sig;

    // Range of wavenumbers to be considered as centering sigc
    CFreal m_hiSigc;
    CFreal m_lowSigc;

    CFuint m_hiBand;
    CFuint m_loBand;

    CFuint m_tempLocalCellID;

    CFreal m_Q;
    CFreal m_ni;
    CFreal m_f2;

    double m_f1;
};

#endif // ATOMIC_LINES_H


}
}


