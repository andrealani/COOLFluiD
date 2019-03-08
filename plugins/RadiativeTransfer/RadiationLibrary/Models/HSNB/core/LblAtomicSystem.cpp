#include "LblAtomicSystem.h"
#include "LblSpectralGrid.h"
#include "Constants.h"
#include "Voigt.h"

#include <string>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <cassert>

using namespace std;
using namespace COOLFluiD::RadiativeTransfer;

void swap(LblAtomicSystem& s1, LblAtomicSystem& s2)
{
    std::swap(s1.m_name, s2.m_name);
    std::swap(s1.m_index, s2.m_index);
    std::swap(s1.m_nlines, s2.m_nlines);
    std::swap(s1.mp_line_data, s2.mp_line_data);
    std::swap(s1.m_lorentz_type, s2.m_lorentz_type);
    std::swap(s1.m_nglor, s2.m_nglor);
    std::swap(s1.m_ntopbase, s2.m_ntopbase);
    std::swap(s1.mp_lorentz_data, s2.mp_lorentz_data);
    std::swap(s1.m_part_func, s2.m_part_func);
    std::swap(s1.m_test_qss, s2.m_test_qss);
}

//==============================================================================

LblAtomicSystem::LblAtomicSystem(const ThermoData& thermo, const std::string& name, const std::string datadir, const bool& qss)
    : m_name(name), m_index(thermo.speciesIndex(name)), m_nlines(0),
      m_part_func(name,datadir), m_test_qss(qss), m_datadir(datadir)
{
    cout << "Loading atomic system: " << m_name << endl;

    loadLineData();
    loadLorentzData();
    if (m_test_qss) {
        loadQssInterface();
        cout << "Non-Boltzmann corrections computed from Johnston's QSS model" << endl;
    }
}

//==============================================================================

LblAtomicSystem::~LblAtomicSystem()
{
    if (mp_line_data != NULL) delete [] mp_line_data;
    if (mp_lorentz_data != NULL) delete [] mp_lorentz_data;
}

//==============================================================================

class PlusTimesEq
{
public:
    PlusTimesEq(const double& a) : m_a(a) {}
    void operator () (double* const y, const double& x, const int& i) const {
        y[i] += m_a*x;
    }
private:
    const double& m_a;
};

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

class PlusTimesEq3
{
public:
    PlusTimesEq3(const double& a) : m_a(a) {}
    void operator () (double* const y, const double& x, const int& i) const {
        y[2*i] += m_a*x;
    }
private:
    const double& m_a;
};

void LblAtomicSystem::spectra(
    ThermoData& thermo, const LblSpectralGrid& grid, double* const p_spectra,
    const SpectralDataType type)
{
	// Thermodynamic state variables
	const double Th  = thermo.Th();
	const double Tel = thermo.Tel();
	const double Te  = thermo.Te();
	const double mw  = thermo[m_index].mw; // [kg/mol]
	const double ni  = thermo.N(m_index); // [m-3]

	cout << "Computing spectra for " << m_name << ": ";

	// Skip this system if the number density is way too low to matter
	if (ni < 1000.0) {
		cout << "number density too small, skipping..." << endl;
		return;
	}

	// Compute the partition function and constant factors
	const double Q = m_part_func.Q(thermo);
	cout << "N = " << ni << ", Q = " << Q << endl;

	const double fac1 = QE*QE*HP*ni/(2.0*ME*EPS0*Q); // [W/sr]
	const double fac2 = HP*C0*100.0/(KB*Tel); // [cm]
	const double fac3 = 2.0*HP*C0*C0*10000.0; // [W cm^2 sr-1]

	// Constant factors for doppler and lorentz line-widths
	const double doppler_fac = std::sqrt(2.0*RU*Th*LN2/(mw*C0*C0));

	double lorentz_fac [11];
	double Tcal = min(max(Th/1000.0, 2.0), 25.0);
	double omega = Tcal - (int) Tcal;
	int ioni = (int) Tcal - 2;
	if ((int) Tcal == 25) {
		omega = 1.0;
		ioni  = 22;
	}

	switch (m_lorentz_type) {
	case NEUTRAL:
		lorentz_fac[8] = 0.5*std::pow(Th/2000.0, 0.3);
		lorentz_fac[0] = lorentz_fac[8]*thermo.N(thermo.speciesIndex("O"));
		lorentz_fac[1] = lorentz_fac[8]*thermo.N(thermo.speciesIndex("N"));
		lorentz_fac[2] = lorentz_fac[8]*thermo.N(thermo.speciesIndex("O2"));
		lorentz_fac[3] = lorentz_fac[8]*thermo.N(thermo.speciesIndex("N2"));
		lorentz_fac[4] = lorentz_fac[8]*thermo.N(thermo.speciesIndex("NO"));
		lorentz_fac[5] = lorentz_fac[8]*thermo.N(thermo.speciesIndex("C"));
		lorentz_fac[6] = lorentz_fac[8]*thermo.N(thermo.speciesIndex("CO"));
		lorentz_fac[7] = lorentz_fac[8]*thermo.N(thermo.speciesIndex("CO2"));
		lorentz_fac[8] = 0.0; // ignoring this factor
		lorentz_fac[10] = 0.5e-22*thermo.N(thermo.speciesIndex("e-"));
		lorentz_fac[9] = lorentz_fac[10]*std::pow(Te/2000.0, 1.0/6.0);
		break;
	case SINGLE_ION:
		lorentz_fac[0] = 1.0e-22*thermo.N(thermo.speciesIndex("e-"));
		break;
	default:
		cout << "this lorentz type is not supported yet!" << endl;
		//exit(1);
	}

    double sigc, sigc3, eta, kappa, gamd, gaml;
    double *pg;

    // Compute Desequilibrium coefficients
    double* p_cu = new double[m_nlines];
    double* p_cl = new double[m_nlines];
    fill(p_cu, p_cu+m_nlines, 1.0);
    fill(p_cl, p_cl+m_nlines, 1.0);

    if (m_test_qss) {
        QssAtoms qss(m_name, m_datadir);
        qss.computeCorrections(thermo);
        for (int i = 0; i < m_nlines; i++) {
            p_cu[i] = qss.getCorrection(mp_line_data[i].Gqssu);
            p_cl[i] = qss.getCorrection(mp_line_data[i].Gqssl);
        }
        //for (int i=0; i<6; i++)
        //    cout << "GROUP " << i << " - " << qss.getCorrection(i) << endl;
    }

    #pragma omp parallel for private(sigc, sigc3, eta, kappa, gamd, gaml, pg)
    for (int i = 0; i < m_nlines; ++i) {
        sigc = mp_line_data[i].sigma; // center of the line [cm-1]
        sigc3 = sigc*sigc*sigc; // [cm-3]

        // Compute eta and kappa
        eta = // [W cm-3 sr-1]
            fac1*sigc3*std::exp(mp_line_data[i].lngf-fac2*mp_line_data[i].Eu)*p_cu[i];
        kappa = // [cm-2]
            eta/(fac3*sigc3)*(std::exp(fac2*sigc)*p_cl[i]/p_cu[i]-1.0);
        // Warning: huge value of fac2*sigc (high wavenumer, low temperature) => infinite exponential factor 
        if (std::isnan(kappa) || std::isinf(kappa)) kappa =0.0;
        //cout << "line #: " << i+1 << endl;
        //cout << "fac1 = " << fac1 << ", fac2 = " << fac2 << ", sigc3 = " << sigc3 << endl;
        //cout << "eta = " << eta << ", kappa = " << kappa << endl;

        // Compute the Doppler and Lorentz HWHM in [cm-1]
        gamd = sigc*doppler_fac;
        gaml = 0.0;
        pg = mp_lorentz_data + m_nglor*mp_line_data[i].topbase;

        switch (m_lorentz_type) {
        case NEUTRAL:
            for (int k = 0; k < 11; ++k)
                gaml += std::abs(pg[k])*lorentz_fac[k];
            break;
        case SINGLE_ION:
            gaml = lorentz_fac[0]*(pg[ioni]*(1.0-omega)+pg[ioni+1]*omega);
            break;
        default:
            cout << "Warning: only neutral or single ion case supported..." << endl;
        }

        //cout << "gamd = " << gamd << ", gaml = " << gaml << endl;

        // Multiply by line-shape and add to the spectra
        switch (type) {
        case ETA:
            Voigt::drayson(grid, gamd, gaml, sigc, p_spectra, PlusTimesEq(eta));
            break;
        case KAPPA:
            Voigt::drayson(grid, gamd, gaml, sigc, p_spectra, PlusTimesEq(kappa));
            break;
        case ETA_KAPPA:
            Voigt::drayson(grid, gamd, gaml, sigc, p_spectra, PlusTimesEq2(eta, kappa));
            break;
        }

        //cout << endl;
    }

    delete[] p_cu;
    delete[] p_cl;
}

//==============================================================================

// This way of computing spectra ensure that equilibrium is retrieved
void LblAtomicSystem::spectraEq(ThermoData &thermo, const LblSpectralGrid& grid, double* const p_spectra)
{
	// Thermodynamic state variables
	const double Th  = thermo.Th();
	const double Tel = thermo.Tel();
	const double Te  = thermo.Te();
	const double mw  = thermo[m_index].mw; // [kg/mol]
	const double ni  = thermo.N(m_index); // [m-3]

	cout << "Computing spectra for " << m_name << ": ";

	// Skip this system if the number density is way too low to matter
	if (ni < 1000.0) {
		cout << "number density too small, skipping..." << endl;
		return;
	}

	// Compute the partition function and constant factors
	const double Q = m_part_func.Q(thermo);
	cout << "N = " << ni << ", Q = " << Q << endl;

	const double fac1 = QE*QE*HP*ni/(2.0*ME*EPS0*Q); // [W/sr]
	const double fac2 = HP*C0*100.0/(KB*Tel); // [cm]
	const double fac3 = 2.0*HP*C0*C0*10000.0; // [W cm^2 sr-1]

	// Constant factors for doppler and lorentz line-widths
	const double doppler_fac = std::sqrt(2.0*RU*Th*LN2/(mw*C0*C0));

	double lorentz_fac [11];
	double Tcal = min(max(Th/1000.0, 2.0), 25.0);
	double omega = Tcal - (int) Tcal;
	int ioni = (int) Tcal - 2;
	if ((int) Tcal == 25) {
		omega = 1.0;
		ioni  = 22;
	}

	switch (m_lorentz_type) {
	case NEUTRAL:
		lorentz_fac[8] = 0.5*std::pow(Th/2000.0, 0.3);
		lorentz_fac[0] = lorentz_fac[8]*thermo.N(thermo.speciesIndex("O"));
		lorentz_fac[1] = lorentz_fac[8]*thermo.N(thermo.speciesIndex("N"));
		lorentz_fac[2] = lorentz_fac[8]*thermo.N(thermo.speciesIndex("O2"));
		lorentz_fac[3] = lorentz_fac[8]*thermo.N(thermo.speciesIndex("N2"));
		lorentz_fac[4] = lorentz_fac[8]*thermo.N(thermo.speciesIndex("NO"));
		lorentz_fac[5] = lorentz_fac[8]*thermo.N(thermo.speciesIndex("C"));
		lorentz_fac[6] = lorentz_fac[8]*thermo.N(thermo.speciesIndex("CO"));
		lorentz_fac[7] = lorentz_fac[8]*thermo.N(thermo.speciesIndex("CO2"));
		lorentz_fac[8] = 0.0; // ignoring this factor
		lorentz_fac[10] = 0.5e-22*thermo.N(thermo.speciesIndex("e-"));
		lorentz_fac[9] = lorentz_fac[10]*std::pow(Te/2000.0, 1.0/6.0);
		break;
	case SINGLE_ION:
		lorentz_fac[0] = 1.0e-22*thermo.N(thermo.speciesIndex("e-"));
		break;
	default:
		cout << "This Lorentz type is not supported yet!" << endl;
		//exit(1);
	}

    double sigc, eta, gamd, gaml;
    double *pg;

    #pragma omp parallel for private(sigc, eta, gamd, gaml, pg)
    for (int i = 0; i < m_nlines; ++i) {
        sigc = mp_line_data[i].sigma; // center of the line [cm-1]
        
        // Compute eta 
        assert(mp_line_data[i].lngf-fac2*mp_line_data[i].Eu < MAXLOG);
        eta = // [W cm-2 sr-1]
            fac1*sigc*sigc*std::exp(mp_line_data[i].lngf-fac2*mp_line_data[i].Eu);

        // Compute the Doppler and Lorentz HWHM in [cm-1]
        gamd = sigc*doppler_fac;
        gaml = 0.0;
        pg = mp_lorentz_data + m_nglor*mp_line_data[i].topbase;

        switch (m_lorentz_type) {
        case NEUTRAL:
            for (int k = 0; k < 11; ++k)
                gaml += std::abs(pg[k])*lorentz_fac[k];
            break;
        case SINGLE_ION:
            gaml = lorentz_fac[0]*(pg[ioni]*(1.0-omega)+pg[ioni+1]*omega);
            break;
        default:
            cout << "This Lorentz type is not supported yet!" << endl;
        }

        // Multiply by line-shape and add to the spectra
        Voigt::drayson(grid, gamd, gaml, sigc, p_spectra, PlusTimesEq3(eta));

    }

}

void LblAtomicSystem::getLineData(ThermoData& thermo, std::vector<AtomicLineData>& line_data)
{
    // Thermodynamic state variables
    const double Th  = thermo.Th();
    const double Tel = thermo.Tel();
    const double Te  = thermo.Te();
    const double mw  = thermo[m_index].mw; // [kg/mol]
    const double ni  = thermo.N(m_index); // [m-3]

    cout << "Computing line data for " << m_name << ": ";

    // Skip this system if the number density is way too low to matter
    if (ni < 1000.0) return;

    // Compute the partition function and constant factors
    const double Q = m_part_func.Q(thermo);
    cout << "N = " << ni << ", Q = " << Q << endl;

    const double fac1 = QE*QE*HP*ni/(2.0*ME*EPS0*Q); // [W/sr]
    const double fac2 = HP*C0*100.0/(KB*Tel); // [cm]
    const double fac3 = 2.0*HP*C0*C0*10000.0; // [W cm^2 sr-1]

    // Constant factors for doppler and lorentz line-widths
    const double doppler_fac = std::sqrt(2.0*RU*Th*LN2/(mw*C0*C0));

    double lorentz_fac [11];
    double Tcal = min(max(Th/1000.0, 2.0), 25.0);
    double omega = Tcal - (int) Tcal;
    int ioni = (int) Tcal - 2;
    if ((int) Tcal == 25) {
        omega = 1.0;
        ioni  = 22;
    }

    switch (m_lorentz_type) {
    case NEUTRAL:
        lorentz_fac[8] = 0.5*std::pow(Th/2000.0, 0.3);
        lorentz_fac[0] = lorentz_fac[8]*thermo.N(thermo.speciesIndex("O"));
        lorentz_fac[1] = lorentz_fac[8]*thermo.N(thermo.speciesIndex("N"));
        lorentz_fac[2] = lorentz_fac[8]*thermo.N(thermo.speciesIndex("O2"));
        lorentz_fac[3] = lorentz_fac[8]*thermo.N(thermo.speciesIndex("N2"));
        lorentz_fac[4] = lorentz_fac[8]*thermo.N(thermo.speciesIndex("NO"));
        lorentz_fac[5] = lorentz_fac[8]*thermo.N(thermo.speciesIndex("C"));
        lorentz_fac[6] = lorentz_fac[8]*thermo.N(thermo.speciesIndex("CO"));
        lorentz_fac[7] = lorentz_fac[8]*thermo.N(thermo.speciesIndex("CO2"));
        lorentz_fac[8] = 0.0; // ignoring this factor
        lorentz_fac[10] = 0.5e-22*thermo.N(thermo.speciesIndex("e-"));
        lorentz_fac[9] = lorentz_fac[10]*std::pow(Te/2000.0, 1.0/6.0);
        break;
    case SINGLE_ION:
        lorentz_fac[0] = 1.0e-22*thermo.N(thermo.speciesIndex("e-"));
        break;
    default:
        cout << "this lorentz type is not supported yet!" << endl;
        //exit(1);
    }

    for (int i = 0; i < m_nlines; ++i) {
        double sigc = mp_line_data[i].sigma; // center of the line [cm-1]

        // Compute eta
        double eta = // [W cm-2 sr-1]
            fac1*sigc*sigc*std::exp(mp_line_data[i].lngf-fac2*mp_line_data[i].Eu);

        // Compute the Doppler and Lorentz HWHM in [cm-1]
        double gamd = sigc*doppler_fac;
        double gaml = 0.0;
        double* pg = mp_lorentz_data + m_nglor*mp_line_data[i].topbase;

        switch (m_lorentz_type) {
        case NEUTRAL:
            for (int k = 0; k < 11; ++k)
                gaml += std::abs(pg[k])*lorentz_fac[k];
            break;
        case SINGLE_ION:
            gaml = lorentz_fac[0]*(pg[ioni]*(1.0-omega)+pg[ioni+1]*omega);
            break;
        default:
            cout << "This Lorentz type is not supported yet!" << endl;
        }

        // Add the line data to the list
        line_data.push_back(AtomicLineData(gaml, gamd, sigc, eta));
    }
}

//==============================================================================

void LblAtomicSystem::loadLineData()
{
    // First open the NIST database file
//    std::string datadir = std::string(std::getenv("HTGR_DATA_DIRECTORY"));
    string database = m_datadir + "/lbl_data/lines/atoms/" + m_name + ".nist";
    ifstream file(database.c_str());
    
    if (!file.is_open()) {
        cout << "Could not open line database file " << database << "!" << endl;
        exit(1);
    }
    
    // Read the first line for the number of atomic lines in the database
    file >> m_nlines;
    
    // Allocate storage for line data
    mp_line_data = new LineData [m_nlines];
    
    // Read in the raw line data
    char buff[1024];
    
    m_ntopbase = 0;
    for (size_t i = 0; i < m_nlines; ++i) {
        file >> mp_line_data[i].topbase
             >> mp_line_data[i].sigma
             >> mp_line_data[i].lngf
             >> mp_line_data[i].Eu;
        file.getline(buff, 1024, '\n');
        m_ntopbase = std::max(mp_line_data[i].topbase, m_ntopbase);
    }
    
    // Close the file
    file.close();
    
    // Convert line parameters to correct units
    const double ln10 = std::log(10.0e0);
    for (int i = 0; i < m_nlines; ++i) {
        LineData& data = mp_line_data[i];
        data.topbase = (data.topbase > 0 ? data.topbase-1 : 0); // shift index
        data.sigma = 1.0e8/data.sigma; // wavelength in A -> wavenumber in cm-1
        data.lngf *= ln10; // log10(gl*flu) -> ln(gl*flu)
    }
}

void LblAtomicSystem::loadLorentzData()
{
    if (m_name == "H") return;

    // Open the Lorentz line-width database for this atom
//    std::string datadir = std::string(std::getenv("HTGR_DATA_DIRECTORY"));
    string database = m_datadir + "/lbl_data/widths/atoms/" + m_name + ".lor";
    ifstream file(database.c_str());
    
    if (!file.is_open()) {
        cout << "Could not open database file " << database << "!" << endl;
        exit(1);
    }
    
    // First read the type of lorentz data
    int lorentz_type;
    file >> lorentz_type;
    m_lorentz_type = static_cast<LorentzType>(lorentz_type);

    // Allocate storage for Lorentz parameters
    switch (m_lorentz_type) {
    case NEUTRAL:
        m_nglor = 11; break;
    case SINGLE_ION:
        m_nglor = 24; break;
    case HYDROGEN:
        m_nglor = 2; break;
    default:
        cout << "Lorentz type not implemented!" << endl;
        exit(1);
    }
    mp_lorentz_data = new double [m_ntopbase*m_nglor];
    
    // Read in the lorentz data
    double* pg;
    size_t i; file >> i;
    while (i <= m_ntopbase && !file.eof()) {
        pg = mp_lorentz_data + (i-1)*m_nglor;
        for (int k = 0; k < m_nglor; ++k)
            file >> pg[k];
        file >> i;
    }
    
    // Close Lorentz file
    file.close();
}

void LblAtomicSystem::loadQssInterface()
{
    // Open the index file interfacing NIST transitions with the groups defining
    // by Johnston in his QSS model
//    std::string datadir = std::string(std::getenv("HTGR_DATA_DIRECTORY"));
    string database = m_datadir + "/qss/atoms/" + m_name + ".tr.qss";
    ifstream file(database.c_str());
    
    if (!file.is_open()) {
        cout << "Could not open line database file " << database << "!" << endl;
        exit(1);
    }

    // Read the indices
    for (size_t i = 0; i < m_nlines; ++i)
        file >> mp_line_data[i].Gqssu
             >> mp_line_data[i].Gqssl;
    
    // Close the file
    file.close();
}

//==============================================================================
