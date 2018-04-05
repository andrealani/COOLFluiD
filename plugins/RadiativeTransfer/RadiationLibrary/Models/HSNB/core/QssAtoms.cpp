#include "QssAtoms.h"
#include <limits>
#include <string>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>

using namespace std;

QssAtoms::QssAtoms(const std::string& name_species)
        : m_species(name_species)
{
    if (m_species == "N" || m_species == "O") {
        cout << "QSS model for " << m_species << " atomic populations" << endl;
    } else {
        cout << "QSS model is not implemented for species " << m_species << endl;
        exit(1);
    }
    
    readJohnstonParameters();
}

QssAtoms::~QssAtoms()
{
    delete [] mp_r0_params;
    delete [] mp_r1_params;
    delete [] mp_fac;
}

void QssAtoms::computeCorrections(ThermoData& thermo)
{
    // Thermodynamic data
    const size_t index = thermo.speciesIndex(m_species);
    const size_t index_ion = thermo.speciesIndex(m_species + "+");
    const size_t index_elec = thermo.speciesIndex("e-");
    double Tel = thermo.Tel();
    double N = thermo.N(index);
    double Nion = thermo.N(index_ion);
    double Ne = thermo.N(index_elec);
    AtomicPartFunc part(m_species), part_ion(m_species + "+");
    double Qat = part.Q(thermo), Qion = part_ion.Q(thermo);
    double xi = std::pow((TWOPI * ME * KB * Tel )/ (HP * HP), 1.5);
    double eion;
    if (m_species == "N") {
        eion = 117225 * HP * C0 * 100.0 ; // cm-1 => J 
    } else if (m_species == "O") {
        eion = 109837 * HP * C0 * 100.0 ; // cm-1 => J
    }

    // Computes population ratio betwwen Saha-Boltzmann and Boltzmann
    // (equivalent to the chineq factor for photoionization)
    double chi = Nion*Ne*exp(eion/(KB*Tel))*Qat/(N*Qion*xi*2.0);

    // Computes Johnston parameters
    double* p_r0 = new double [m_ngroups];
    double* p_r1 = new double [m_ngroups];
    fitJohnston(Tel, Ne, p_r0, p_r1);

    // Computes the factor to apply to the equilibrium population 
    // of each electronic state belonging to the group j
    mp_fac[0] = 1.0;
    for (int j=0; j<m_ngroups; j++)
        mp_fac[j+1] = min(p_r0[j]*chi + p_r1[j], 1.0);

}

void QssAtoms::fitJohnston(double const Te, double const Ne,
                           double* const p_r0, double* const p_r1)
{
    // Clip the Te and Ne ranges
    double T = max( min(Te, 14000.0), 7000.0);
    double N = max( min(Ne*1.0e-6,  1.0e16), 1.0e14); // part/cm3

    double T2 = T*T;
    double T3 = T2*T;
    double lnN = log(N);
    double lnN2 = lnN*lnN;
    double lnN3 = lnN2*lnN;

    double fac1, fac2, fac3, fac4;

    // Compute the fit for population according to Johnston
    // Appendix C, p212-214 of his thesis
    for (int group=0; group < m_ngroups; group++ ) {

        fac1 = (mp_r0_params[ 0+group*16]*T3
             +  mp_r0_params[ 1+group*16]*T2
             +  mp_r0_params[ 2+group*16]*T
             +  mp_r0_params[ 3+group*16]) * lnN3;
        fac2 = (mp_r0_params[ 4+group*16]*T3
             +  mp_r0_params[ 5+group*16]*T2
             +  mp_r0_params[ 6+group*16]*T
             +  mp_r0_params[ 7+group*16]) * lnN2;
        fac3 = (mp_r0_params[ 8+group*16]*T3
             +  mp_r0_params[ 9+group*16]*T2
             +  mp_r0_params[10+group*16]*T
             +  mp_r0_params[11+group*16]) * lnN;
        fac4 = (mp_r0_params[12+group*16]*T3
             +  mp_r0_params[13+group*16]*T2
             +  mp_r0_params[14+group*16]*T
             +  mp_r0_params[15+group*16]);
        p_r0[group] = fac1 + fac2 + fac3 + fac4;

        fac1 = (mp_r1_params[ 0+group*16]*T3
             +  mp_r1_params[ 1+group*16]*T2
             +  mp_r1_params[ 2+group*16]*T
             +  mp_r1_params[ 3+group*16]) * lnN3;
        fac2 = (mp_r1_params[ 4+group*16]*T3
             +  mp_r1_params[ 5+group*16]*T2
             +  mp_r1_params[ 6+group*16]*T
             +  mp_r1_params[ 7+group*16]) * lnN2;
        fac3 = (mp_r1_params[ 8+group*16]*T3
             +  mp_r1_params[ 9+group*16]*T2
             +  mp_r1_params[10+group*16]*T
             +  mp_r1_params[11+group*16]) * lnN;
        fac4 = (mp_r1_params[12+group*16]*T3
             +  mp_r1_params[13+group*16]*T2
             +  mp_r1_params[14+group*16]*T
             +  mp_r1_params[15+group*16]);
        p_r1[group] = fac1 + fac2 + fac3 + fac4;

    }

}

void QssAtoms::readJohnstonParameters()
{

    m_ngroups = 5;
    std::string datadir = std::string(std::getenv("HTGR_DATA_DIRECTORY"));
    std::string database;
    std::ifstream file;
    std::ostringstream ss;
    mp_r0_params = new double [16*m_ngroups];
    mp_r1_params = new double [16*m_ngroups];
    mp_fac = new double [m_ngroups+1];
    fill(mp_fac, mp_fac+m_ngroups+1, 1.0);

    for (int group=0; group < m_ngroups; group++) {
        // Open the file where parameters of the group are stored
        ss.str("");
        ss << group+1;
        std::string database = datadir + "/qss/atoms/" + m_species + "_" + ss.str() + ".qss";
        file.open(database.c_str());

        // Read r0 and r1 parameters
        for (int i = 0; i < 16; ++i)
            file >> mp_r0_params[group*16+i];
        for (int i = 0; i < 16; ++i)
            file >> mp_r1_params[group*16+i];
        file.close();
    }

}
