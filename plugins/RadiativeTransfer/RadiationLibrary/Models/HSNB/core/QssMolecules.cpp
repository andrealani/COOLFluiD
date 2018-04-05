#include "QssMolecules.h"
#include <iostream>
#include <fstream>
#include <cstdlib>

using namespace std;

QssMolecules::QssMolecules(const std::string& name_system, const std::string& name_species)
        : m_system(name_system), m_species(name_species), mp_johnston_params(NULL),
          m_upper_energy_index(-1), m_upper_level_char('Z')
{

    // Systems for which Johnston's QSS model is used for computing
    // non-Boltmann electronic population
    if ( (m_system == "N2_1+")
       | (m_system == "N2_2+")
       | (m_system == "N2+_1-")
       | (m_system == "N2+_2-")
       | (m_system == "N2+M")) {
        m_qss_model = JOHNSTON_FIT;
        mp_johnston_params = new double [16];

        // Setup the features for the upper energy level
        // of the transition of the given system
        if (m_system == "N2_1+") {
            m_upper_level_char = 'B';
            m_upper_energy_index = 2;
        } else if (m_system == "N2_2+") {
            m_upper_level_char = 'C';
            m_upper_energy_index = 10;
        } else if (m_system == "N2+M") {
            m_upper_level_char = 'A';
            m_upper_energy_index = 1;
        } else if (m_system == "N2+_1-") {
            m_upper_level_char = 'B';
            m_upper_energy_index = 2;
        } else if (m_system == "N2+_2-") {
            m_upper_level_char = 'C';
            m_upper_energy_index = 12;
        }
        readJohnstonParameters();

    // Systems for which the upper energy level of the transition is above 
    // the dissociation limit (N2-VUV). An equilibrium assumption with the
    // associated atom is used to populate the level.
    } else if ( (m_system == "N2BH1")
              | (m_system == "N2BH2")
              | (m_system == "N2CY")
              | (m_system == "N2WJ")
              | (m_system == "N2W")) {
        m_qss_model = DISS_EQUIL;

    // CN violet
    } else if (m_system == "CNvio") {
        m_qss_model = CN_VIOLET;
    } else {
        m_qss_model = NONE;
    }
    //cout << "QSS model for system " << m_system << " - " << m_qss_model << endl;

}

QssMolecules::~QssMolecules()
{
    if (mp_johnston_params)
        delete [] mp_johnston_params;
}

double QssMolecules::levelCorrection(ThermoData& thermo)
{

    switch (m_qss_model) {

        // Johnston's QSS model is used for computing
        // non-Boltmann electronic population
        case JOHNSTON_FIT:
            return levelCorrectionJohnston(thermo); 
        break;

        // An equilibrium assumption with the associated atom 
        // is used to populate the level (N2-VUV)
        case DISS_EQUIL:
            return levelCorrectionDissEquil(thermo); 
        break;

        case CN_VIOLET:
            return levelCorrectCNViolet(thermo);
        break;

        case NONE:
            return 1.0; 
        break;

    }
}

double QssMolecules::levelCorrectCNViolet(ThermoData& thermo)
{
    double Tel = thermo.Tel();
    double Ta  = std::sqrt(thermo.Tr()*thermo.Tv());
    double f0 = 1.1*std::exp(37049.4/Tel);

    double n = 0.0;
    for (int i = 0; i < thermo.nSpecies(); ++i)
        n += thermo.N(i);
    double ne = n - thermo.N(thermo.speciesIndex("e-"));
    double nh = n - ne;
    double a  = (1.8e5*nh + 6.25e8*ne)/NA*std::sqrt(Ta)*std::exp(-37049.4/Ta);
    double f1 = a / (1.0/62.5e-9 + a*f0);

    double f2 = 2.5*std::exp(-13293.65/Tel);
    double f3 = f0*(1.0 + f2) + 7.0;

    //return (f3*(1-f2/(1+f2))/((1+f1)/f1-f2/(1+f2)));
    return f1*f0;
}

double QssMolecules::levelCorrectionDissEquil(ThermoData& thermo)
{

    double tr = thermo.Tr();
    double tv = thermo.Tv();
    const size_t index_at = thermo.speciesIndex("N");
    const size_t index_mol = thermo.speciesIndex("N2");
    AtomicPartFunc qat("N");
    MolecularPartFunc qmol("N2");

    double q_ratio = qmol.Q(tv,tv) / pow(qat.Q(thermo),2.0);
    double n_ratio = std::pow( thermo.N()[index_at], 2.0)
                   / thermo.N()[index_mol];

    double mn = thermo[index_at].mw / (2*NA); // reduced mass of dissociation products
    double xi = std::pow((TWOPI * mn * KB * tr )/ (HP * HP), 1.5);

    const double ediss = 78711 * HP * C0 * 100 ; // cm-1 => J

    return min(max(n_ratio * q_ratio * std::exp(ediss/(KB*tv)) / xi, 0.0), 10.0);

}

double QssMolecules::levelCorrectionJohnston(const ThermoData& thermo)
{
    const int ie = thermo.speciesIndex("e-");
    double tr = thermo.Tr();
    double tv = thermo.Tv();
    double ne = thermo.N(ie);

    MolecularPartFunc part(m_species);
    MolecularElecPartFunc part_elec(m_species);
    double n_equil = part_elec.Qel(m_upper_energy_index,tr,tv)
                   * exp(-part_elec.Eel(m_upper_energy_index)*HP*C0*100/KB/tv)
                   / part.Q(tr,tv);

    double n_johnston = n_equil;
    if (ie >= 0)
        n_johnston = fitJohnston(tv, ne);
    //cout << "NJ " << tr << " " << tv << " " << ne << " " << n_johnston << " " << n_equil << " " << min(max(n_johnston / n_equil, 0.0), 1.0) << endl;

    return min(max(n_johnston / n_equil, 0.0), 10.0);
}

double QssMolecules::fitJohnston(double const Te, double const Ne)
{

    // Clip the Te and Ne ranges
    double T = max( min(Te, 14000.0), 7000.0);
    double N = max( min(Ne*1.0e-6,  1.0e16), 1.0e14); // part/cm3

    double T2 = T*T;
    double T3 = T2*T;
    double lnN = log(N);
    double lnN2 = lnN*lnN;
    double lnN3 = lnN2*lnN;

    // Compute the fit for population according to Johnston
    // Appendix C, p214 of his thesis
    double fac1 = (mp_johnston_params[ 0]*T3
                +  mp_johnston_params[ 1]*T2
                +  mp_johnston_params[ 2]*T
                +  mp_johnston_params[ 3]) * lnN3;
    double fac2 = (mp_johnston_params[ 4]*T3
                +  mp_johnston_params[ 5]*T2
                +  mp_johnston_params[ 6]*T
                +  mp_johnston_params[ 7]) * lnN2;
    double fac3 = (mp_johnston_params[ 8]*T3
                +  mp_johnston_params[ 9]*T2
                +  mp_johnston_params[10]*T
                +  mp_johnston_params[11]) * lnN;
    double fac4 = (mp_johnston_params[12]*T3
                +  mp_johnston_params[13]*T2
                +  mp_johnston_params[14]*T
                +  mp_johnston_params[15]);
    return fac1 + fac2 + fac3 + fac4;

}

void QssMolecules::readJohnstonParameters()
{

    // Open the file where parameters are stored
    std::string datadir = std::string(std::getenv("HTGR_DATA_DIRECTORY"));
    std::string database = datadir + "/qss/molecules/" + m_species + "_" + m_upper_level_char + ".qss";
    std::ifstream file(database.c_str());

    // Read the parameters
    for (int i = 0; i < 16; ++i)
        file >> mp_johnston_params[i];
    file.close();

}
