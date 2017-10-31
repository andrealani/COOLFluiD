#include "AtomicLines.h"
#include "AtomicPartFunc.h"
#include "Constants.h"
#include "ThermoData.h"
#include "HydrogenPartFunc.h"
#include "LblSpectralGrid.h"
#include "LorentzHydrogen.h"
#include "LorentzTopbase.h"
#include "Voigt.h"
#include "PhotonPath.h"
#include "ThermoData.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

namespace RadiativeTransfer {

//////////////////////////////////////////////////////////////////////////////

//==============================================================================

AtomicLines::AtomicLines(const std::string baseDir,
    const std::string& name, const ThermoData& thermo) :
    m_atom(name),
    m_species_index(thermo.speciesIndex(name)),
    m_mw(0.0),
    m_qss(false),
    m_opt_thick(0.0),
    m_last_loc(0),
    m_last_sig(0.0)
{
    m_datadir=baseDir+"/data";

    // Check if QSS is toggled for this atom
    if (m_atom[m_atom.size()-1] == '*') {
        m_qss = true;
        m_atom.erase(m_atom.size()-1);
    }

    // Check if this atom is present in the thermodynamic database and load
    // species specific constants
    if (m_species_index < 0) {
        std::cout << "Error loading " << m_atom << ": species is not loaded!"
                  << std::endl;
        std::exit(1);
    }

    m_mw = thermo[m_species_index].mw;

    // Load line data and Lorentz model
    if (m_atom == "H")
        loadHData();  // special case for H atoms
    else
        loadNistData();

    // Initialize non-Boltzmann correction factors
    for (int i = 0; i < nLines(); ++i) {
        m_line_data[i].nonbolt_l = 1.0;
        m_line_data[i].nonbolt_u = 1.0;
    }
}

//==============================================================================

void AtomicLines::updateDopplerHWHM(double T)
{
    // Temperature factor
    const double fac = std::sqrt(2.0*RU*T*LN2/(m_mw*C0*C0));

    for (int i = 0; i < nLines(); ++i)
        m_line_data[i].gamd = fac * m_line_data[i].sigc;
}

//==============================================================================

const std::vector<LineData>& AtomicLines::lineData(
    const ThermoData& thermo)
{
    // Update the Doppler broadening parameters
    updateDopplerHWHM(thermo.Th());

    // Update the Lorentz broadening parameters
    mp_lorentz->update(thermo);
    mp_lorentz->setHWHM(m_line_data);
    //for (int i = 0; i < nLines(); ++i)
    //    m_line_data[i].gaml = mp_lorentz->hwhm(i);

    // Update the non-Boltzmann correction factors
    // TBD
    if (m_qss)
        std::cout << "QSS not implemented for atomic lines!" << std::endl;

    return m_line_data;
}

//==============================================================================

// Simple class used in the Voigt evaluation in AtomicLines::spectrum()
class PlusTimesEq2
{
public:
    PlusTimesEq2(const double& a, const double& b) : m_a(a), m_b(b) {}
    void operator () (double* const y, const double& x, const int& i) const {
        y[i*2]   += m_a*x;
        y[i*2+1] += m_b*x;
    }
private:
    const double& m_a;
    const double& m_b;
};

//==============================================================================

void AtomicLines::spectrum(ThermoData& thermo, const LblSpectralGrid& grid,
    double* const p_spectrum)
{
    // Update the line data
    lineData(thermo);

    // Some constants
    const double ni = thermo.N(m_species_index);
    const double f1 = QE*QE*HP*ni/(2.0*ME*EPS0*mp_Q->Q(thermo)); // [W/sr]
    const double f2 = HP*C0*100.0/(KB*thermo.Tel()); // [cm]
    const double f3 = 1.0e-4/(2.0*HP*C0*C0);

    // Initialize spectra with 0s
    std::fill(p_spectrum, p_spectrum + 2*grid.size(), 0.0);

    // Compute line dependent terms
    for (int i = 0; i < nLines(); ++i) {
        const LineData& line = m_line_data[i];
        //if (line.sigc > grid.max() + 10.0*line.gaml) continue;

        // eta / sigma in [W cm-2 sr-1]
        double eta = f1*line.sigc*line.sigc*std::exp(line.lnglflu-f2*line.Eu);

        // Apply line profile
        // eta_ul / sigma * f(sigma - sigma_ul)
        // eta_ul / sigma * f(sigma - sigma_ul) * theta_l / theta_u
        Voigt::drayson(
            grid, line.gamd, line.gaml, line.sigc, p_spectrum,
            PlusTimesEq2(eta*line.nonbolt_u, eta*line.nonbolt_l));
    }

    double sig;
    int ieta, ikap;
    double f;

    // Finalize spectrum with grid only dependent terms
    for (int i = 0, ie = 0, ik = 1; i < grid.size(); ++i, ie+=2, ik+=2) {
        sig = grid[i];

        // Split exponential to avoid overflow at high wavenumbers
        f = exp(0.5*sig*f2);
        p_spectrum[ik] = f3*f*(f*p_spectrum[ik]-p_spectrum[ie]/f)/(sig*sig);
        p_spectrum[ie] *= sig;
    }
}

//==============================================================================

void AtomicLines::thinLineEmission(ThermoData& thermo, double* const p_emis, CFreal& sumEmission, bool printDebug)
{
    // Update the line data
    lineData(thermo);

    // Some constants
    const double ni = thermo.N(m_species_index);
    const double f1 = QE*QE*HP*ni/(2.0*ME*EPS0*mp_Q->Q(thermo)); // [W/sr]
    const double f2 = HP*C0*100.0/(KB*thermo.Tel()); // [cm]

    // Compute line dependent terms
    for (int i = 0; i < nLines(); ++i) {
        const LineData& line = m_line_data[i];
        const double sigc = line.sigc;
        p_emis[i] = f1*sigc*sigc*sigc*std::exp(line.lnglflu-f2*line.Eu);

        if (printDebug) {
            CFLog(INFO, "AtomicLines::thinLineEmission => p_emis[" << i << "]=" << p_emis[i] << " \n");
            CFLog(INFO, "AtomicLines::thinLineEmission => f1=" << f1 << ", sigc=" << sigc << " \n");
            CFLog(INFO, "AtomicLines::thinLineEmission => ni=" << ni << ", mp_Q->Q(thermo)=" << mp_Q->Q(thermo) << " \n");
            CFLog(INFO, "AtomicLines::thinLineEmission => " << m_atom<<", m_species_index=" << m_species_index << ", mp_Q->Q(thermo)=" << mp_Q->Q(thermo) << " \n");
        }

        sumEmission+=p_emis[i];
    }
}

//==============================================================================

//double AtomicLines::tau(
//    const FieldData& field, const LineOfSight& los, double sig)
//{
//   return std::exp(-opticalThickness(field, los, sig));
//}

//==============================================================================

//double AtomicLines::opticalThickness(
//    const FieldData& field, const LineOfSight& los, double sig)
//{
//    ThermoData& thermo = ThermoData::getInstance();
//
//    double tau = 0.0;
//
//    // Loop over each cell in the los
//    for (int i = 0; i < los.nPath(); ++i) {
//        // Update line data for current point
//        thermo.setState(
//            field.cell(i).Tr(), field.cell(i).Tv(), field.cell(i).P(),
//            field.cell(i).X());
//        lineData(thermo);
//
//        const double ni = thermo.N(m_species_index);
//        const double f1 = QE*QE*HP*ni/(2.0*ME*EPS0*mp_Q->Q(thermo)); // [W/sr]
//        const double f2 = HP*C0*100.0/(KB*thermo.Tel()); // [cm]
//        const double f3 = 1.0e-4/(2.0*HP*C0*C0);
//
//        // Split exponential to avoid overflow at high wavenumbers
//        double f = exp(0.5*sig*f2);
//
//        // Loop over each line and compute the absorption coefficient at sig
//        double kappa = 0.0;
//        for (int j = 0; j < nLines(); ++j) {
//            // Get line data
//            const LineData& line = m_line_data[j];
//            //if (std::abs(line.sigc - sig) > 10*(line.gamd+line.gaml)) continue;
//            // Voigt
//            double v = Voigt::humlicek(line.gaml, line.gamd, sig-line.sigc);
//            // eta / sigma in [W cm-2 sr-1]
//            double eta = f1*line.sigc*line.sigc*std::exp(line.lnglflu-f2*line.Eu);
//            kappa += f3*v*eta*f*(f*line.nonbolt_l/line.nonbolt_u-line.nonbolt_u/f)/(sig*sig);
//        }
//
//        // Compute optical length
//        tau += kappa*los.length(i);
//    }
//
//    return tau;
//}

double AtomicLines::opticalThickness(ThermoData &thermo, const PhotonPath& path, int ic, double sig)
{
    assert(ic <= path.nCells());
    const double f1 = 0.25e-4*QE*QE/(C0*C0*ME*EPS0*sig*sig);

    // If the point in which we want the optical thickness is less than the last
    // point we computed, start from the beginning.  In addition, if a -1 is
    // given for ic, then compute the whole path from the beginning.
    if (ic < m_last_loc || sig != m_last_sig) {
        if (ic == -1) ic = path.nCells();
        m_opt_thick = 0.0;
        m_last_loc = 0;
        m_last_sig = sig;
    }


    double Q, ni, f2, f, sum1, sum2;

    // Loop over each cell in the los
    for (int i = m_last_loc; i < ic; ++i) {
        // Update line data for current point
        int id = path.cellID(i);
        thermo.setState(i);
        lineData(thermo);

        f2 = HP*C0*100.0/(KB*thermo.Tel()); // [cm]

        double sum1 = 0, sum2 = 0, f;
        for (int j = 0; j < nLines(); ++j) {
            const LineData& line = m_line_data[j];
            f = Voigt::humlicek(line.gaml, line.gamd, sig-line.sigc);
            f *= line.sigc*line.sigc*std::exp(line.lnglflu-f2*line.Eu);
            sum1 += f*line.nonbolt_l;
            sum2 += f*line.nonbolt_u;
        }

        // Split exponential to avoid overflow at high wavenumbers
        f = exp(0.5*sig*f2);
        Q  = m_loc_params[id].Q;
        ni = m_loc_params[id].N;
        m_opt_thick += f1*ni/Q*f*(f*sum1 - sum2/f)*path.cellDistance(i);
    }

    m_last_loc = ic;
    return m_opt_thick;
}

double AtomicLines::opticalThickness(HSNBAtomicParameterSet& atomParams)
{
    //CFLog(VERBOSE, "AtomicLines::opticalThickness => opticalThickness =" << atomParams.optThick<< "\n" );
    return atomParams.optThick;
}

void AtomicLines::addStateParams(ThermoData &thermo,HSNBAtomicParameterSet &atomParams, CFreal cellDistance, CFuint localCellID, CFreal sig)
{


    const double f1 = 0.25e-4*QE*QE/(C0*C0*ME*EPS0*sig*sig);

    double Q, ni, f2;

    // Assume thermoState has been set outside of this function

    lineData(thermo);

    f2 = HP*C0*100.0/(KB*thermo.Tel()); // [cm]

    double sum1 = 0, sum2 = 0, f;
    for (int j = 0; j < nLines(); ++j) {
        const LineData& line = m_line_data[j];
        f = Voigt::humlicek(line.gaml, line.gamd, sig-line.sigc);
        f *= line.sigc*line.sigc*std::exp(line.lnglflu-f2*line.Eu);
        sum1 += f*line.nonbolt_l;
        sum2 += f*line.nonbolt_u;
    }

    // Split exponential to avoid overflow at high wavenumbers
    f = exp(0.5*sig*f2);
    Q  = m_loc_params[localCellID].Q;
    ni = m_loc_params[localCellID].N;

//    CFLog(VERBOSE, "AtomicLines::addStateParams => Q= " << Q << ", ni="<< ni <<" \n");

    atomParams.optThick+=f1*ni/Q*f*(f*sum1 - sum2/f)*cellDistance;

    if (std::isnan(atomParams.optThick)) {
        std::cout << "Computed singular optical thickness, printing telemetry: \n";
        std::cout << "Enum: f1=" << f1 << ", sig=" << sig << ", ni=m_loc_params[" << localCellID << "].N=" << ni << "\n";
        std::cout << "Denum: f2=" << f2 << ", thermo.Tel()=" << thermo.Tel() << ", f=" << f << ", sum1=" << sum1
                  << ", sum2=" << sum2 << " cellDistance=" << cellDistance << "\n";
    }

//    CFLog(VERBOSE, "AtomicLines::addStateParams => f1=" << f1 << ", f="<< f << ", f2=" << f2 << ", sum1=" << sum1 << ", sum2=" <<sum2 << "\n");
//    std::cout << "Enum: f1=" << f1 << ", sig=" << sig << ", ni=m_loc_params[" << localCellID << "].N=" << ni << "\n";
//    std::cout << "Denum: f2=" << f2 << ", thermo.Tel()=" << thermo.Tel() << ", f=" << f << ", sum1=" << sum1
//              << ", sum2=" << sum2 << " cellDistance=" << cellDistance << "\n";
    //    //DELETE LATER
//    CFLog(VERBOSE, "AtomicLines::addStateParams => optThick= " << f1*ni/Q*f*(f*sum1 - sum2/f)*cellDistance << " \n");
}



void AtomicLines::setupLocalParameters(ThermoData &thermo)
{
    m_loc_params.resize(thermo.nCells());

    for (int i = 0; i < thermo.nCells(); ++i) {
        thermo.setState(i);
        m_loc_params[i].Q = mp_Q->Q(thermo);      
        m_loc_params[i].N = thermo.N(m_species_index);
    }
}

//==============================================================================

void AtomicLines::loadNistData()
{
    // First open the NIST database file
    std::string database = m_datadir + "/lbl_data/lines/atoms/" + m_atom + ".nist";
    std::ifstream file(database.c_str());

//    std::cout << "Loading line data for atom " + m_atom << std::endl;

    if (!file.is_open()) {
        std::cout << "Could not open line database file " << database << "!"
                  << std::endl;
        exit(1);
    }

    // Read the first line for the number of atomic lines in the database
    int nlines;
    file >> nlines;

    // Allocate storage for line data
    m_line_data.resize(nlines);
    std::vector<int> topbase_index(nlines);

    // Read in the raw line data
    char buff[1024];

    for (size_t i = 0; i < nlines; ++i) {
        file >> topbase_index[i]
             >> m_line_data[i].sigc    // [A]
             >> m_line_data[i].lnglflu
             >> m_line_data[i].Eu;      // [cm-1]
        file.getline(buff, 1024, '\n');
    }

    // Close the file
    file.close();

    // Some double checking
    for (size_t i = nlines-2; i--; )
        assert(m_line_data[i].sigc <= m_line_data[i+1].sigc);

    // Convert line parameters to correct units
    const double ln10 = std::log(10.0);
    for (int i = 0; i < nlines; ++i) {
        // Shift TOPBASE index
        topbase_index[i] = std::max(0, topbase_index[i] - 1);

        // wavelength in A -> wavenumber in cm-1
        m_line_data[i].sigc = 1.0e8 / m_line_data[i].sigc;

        // log10(gl*flu) -> ln(gl*flu)
        m_line_data[i].lnglflu *= ln10;
    }

    // Load the Lorentz broadening model
    mp_lorentz = HSNBSharedPtr<LorentzModel>(
        new LorentzTopbase(m_atom, topbase_index));

    // Load the partition function
    mp_Q = HSNBSharedPtr<PartitionFunction>(new AtomicPartFunc(m_atom));
}

//==============================================================================

void AtomicLines::loadHData()
{
    // First open the NIST database file

    std::string database = m_datadir + "/lbl_data/lines/atoms/H.nist";
    std::ifstream file(database.c_str());

    if (!file.is_open()) {
        std::cout << "Could not open line database file " << database << "!"
                  << std::endl;
        exit(1);
    }

    // Read the first line for the number of atomic lines in the database
    int nlines;
    file >> nlines;

    // Allocate storage for line data
    m_line_data.resize(nlines);
    std::vector< std::pair<int, int> > transitions(nlines);

    // Read in the raw line data
    char buff[1024];

    for (size_t i = 0; i < nlines; ++i) {
        file >> m_line_data[i].sigc    // [A]
             >> m_line_data[i].lnglflu
             >> m_line_data[i].Eu      // [cm-1]
             >> transitions[i].first
             >> transitions[i].second;
        file.getline(buff, 1024, '\n');
    }

    // Close the file
    file.close();

    // Convert line parameters to correct units
    const double ln10 = std::log(10.0);
    for (int i = 0; i < nlines; ++i) {
        // wavelength in A -> wavenumber in cm-1
        m_line_data[i].sigc = 1.0e8 / m_line_data[i].sigc;

        // log10(gl*flu) -> ln(gl*flu)
        m_line_data[i].lnglflu *= ln10;

        //std::cout << m_line_data[i].sigc << " " << m_line_data[i].lnglflu << std::endl;
    }

    // Load the Lorentz broadening model
    mp_lorentz = HSNBSharedPtr<LorentzModel>(
        new LorentzHydrogen(m_atom, transitions));

    // Load the partition function
    mp_Q = HSNBSharedPtr<PartitionFunction>(new HydrogenPartFunc());
}

//==============================================================================

void AtomicLines::save(const std::string& file_name) const
{
    static const double ln10 = std::log(10.0);
    std::ofstream f(file_name.c_str());

    f << nLines() << "\n";
    for (int i = 0; i < nLines(); ++i) {
        f << std::setw(4) << dynamic_cast<const LorentzTopbase*>(
                mp_lorentz.operator->())->topbaseID(i) + 1;
        f << std::fixed << std::setprecision(4);
        f << std::setw(15) << 1.0e8 / m_line_data[i].sigc;
        f << std::fixed << std::setprecision(3);
        f << std::setw(9) << m_line_data[i].lnglflu / ln10;
        f << std::fixed << std::setprecision(4);
        f << std::setw(13) << m_line_data[i].Eu;
        f << "\n";
    }

    f.close();
}

}
}

