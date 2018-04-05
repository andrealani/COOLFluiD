#include "Chineq.h"
#include <iostream>
#include <fstream>
#include <cstdlib>

using namespace std;

double Chineq::concNegIon(ThermoData &thermo)
{

    const size_t index_neutral = thermo.speciesIndex(m_species.substr(0,1));
    const size_t index_elec = thermo.speciesIndex("e-");

    AtomicPartFunc qneg(m_species), qat(m_species.substr(0,1));
    double q_ratio = qneg.Qneg(thermo.Tel()) / qat.Q(thermo) ;

    double xi = std::pow((TWOPI * ME * KB * thermo.Tel() )/ (HP * HP), 1.5);
    const double eion  = ionizationEnergy();

    return thermo.N()[index_elec] * thermo.N()[index_neutral]
           * q_ratio * std::exp(eion/(KB*thermo.Tel())) / (2*xi);

}

double Chineq::computeChi(ThermoData& thermo)
{

    if (m_system == "N_bf" || m_system == "O_bf" || m_system == "C_bf" || m_system == "H_bf") {
        return atomicPhotoionization(thermo);
    } else if (m_system == "N2_bf" || m_system == "NO_bf" || m_system == "O2_bf") {
        return molecularPhotoionization(thermo);
    } else if (m_system == "O2_bf_SR") {
        return oxygenPhotodissociation(thermo);
    } else {
        return 1.0;
    }

}

double Chineq::atomicPhotoionization(ThermoData& thermo)
{
    double tel = thermo.Tel();

    const size_t index_atom = thermo.speciesIndex(m_species);
    const size_t index_ion = thermo.speciesIndex(m_species + "+");
    const size_t index_elec = thermo.speciesIndex("e-");

    double q_ratio;
    if (m_species == "H")
        q_ratio = 2.0;
    else
        q_ratio = AtomicPartFunc(m_species).Q(thermo) /
            AtomicPartFunc(m_species + "+").Q(thermo);

    //AtomicPartFunc qat(m_species), qion(m_species + "+");
    //double q_ratio = qat.Q(thermo)
    //               / qion.Q(thermo);
  
    double n_ratio = thermo.N()[index_elec] * thermo.N()[index_ion] 
                   / thermo.N()[index_atom];

    double xi = std::pow((TWOPI * ME * KB * tel )/ (HP * HP), 1.5);

    const double eion  = ionizationEnergy();
    // Correction of the ionization energy from the Debye ionization lowering
    const double alpha = thermo[index_atom].charge;

    double api;
    if (m_species == "H")
        api = 0.0;
    else
        api = (alpha+1.0)*QE*QE/(FOURPI*EPS0*AtomicPartFunc(m_species).debye(thermo));
    //const double rhod  = (qat.debye(thermo);
    //const double api   = (alpha+1.0)*QE*QE/(FOURPI*EPS0*rhod);

    //return n_ratio * q_ratio * std::exp((eion-api)/(KB*tv)) / (2*xi);
    // WARNING: the exponential factor has been splitted in two contribution
    // exp((eion-api)/kTv) = exp(eion/2kTv) * exp (eion/2kTv - api/kTv)
    return n_ratio * q_ratio * std::exp(eion/(2*KB*tel)-api/(KB*tel)) / (2*xi);
                
}

double Chineq::molecularPhotoionization(ThermoData& thermo)
{
    double tv = thermo.Tv();

    const size_t index_mol = thermo.speciesIndex(m_species);
    const size_t index_ion = thermo.speciesIndex(m_species + "+");
    const size_t index_elec = thermo.speciesIndex("e-");

    MolecularPartFunc qmol(m_species), qion(m_species + "+");
    double q_ratio = qmol.Q(tv, tv) / qion.Q(tv, tv);
 
    double n_ratio = thermo.N()[index_elec] * thermo.N()[index_ion] 
                   / std::max(1000.0, thermo.N()[index_mol]);

    double xi = std::pow((TWOPI * ME * KB * tv )/ (HP * HP), 1.5); 

    const double eion = ionizationEnergy();
    // Correction of the ionization energy from the Debye ionization lowering
    const double alpha = thermo[index_mol].charge;
    const double rhod  = qmol.debye(thermo);
    const double api   = (alpha+1.0)*QE*QE/(FOURPI*EPS0*rhod);

    //return n_ratio * q_ratio * std::exp((eion-api)/(KB*tv)) / (2*xi);
    // WARNING: the exponential factor has been splitted in two contribution
    // exp((eion-api)/kTv) = exp(eion/2kTv) * exp (eion/2kTv - api/kTv)
    return n_ratio * q_ratio * std::exp(eion/(2*KB*tv)-api/(KB*tv)) / (2*xi);

}

double Chineq::oxygenPhotodissociation(ThermoData& thermo)
{
    double tr = thermo.Tr();
    double tv = thermo.Tv();

    const size_t index_at = thermo.speciesIndex("O");
    const size_t index_mol = thermo.speciesIndex("O2");

    double qat  = AtomicPartFunc("O").Q(thermo);
    double qmol = MolecularPartFunc("O2").Q(tr, tr);

    double q_ratio = qmol / (qat * qat);

    double n_ratio = thermo.N()[index_at] * thermo.N()[index_at]
                   / thermo.N()[index_mol];

    double Mo = thermo[index_at].mw; // reduced mass of dissociation products
    double xi = TWOPI * Mo * KB * tr / (NA * HP * HP);
    xi = std::sqrt(xi)*xi;

    const double ediss = 41260. * HP * C0 * 100 ; // cm-1 => J
    const double eO1D  = 15867. * HP * C0 * 100 ; // cm-1 => J

    return n_ratio * q_ratio * std::exp((ediss+eO1D)/(KB*tr)-eO1D/(KB*tv)) / xi;
}

double Chineq::ionizationEnergy()
{

  // Ionization energies converted in J
  if (m_species == "N") {
    return 117225 * HP * C0 * 100.0 ; // cm-1 => J 
  } else if (m_species == "O") {
    return 109837 * HP * C0 * 100.0 ; // cm-1 => J
  } else if (m_species == "N-") {
    return 0.1 * QE ; // eV => J
  } else if (m_species == "O-") {
    return 17004.8122 * KB ; // K => J
  } else if (m_species == "N2") {
    return 15.581 * QE ; // eV => J
  } else if (m_species == "O2") {
    return 12.07 * QE ; // eV => J
  } else if (m_species == "NO") {
    return 9.2643 * QE ; // eV => J
  } else if (m_species == "C2") {
    return 11.4 * QE ; // eV => J source: NIST Chemical Webbook
  } else if (m_species == "CO") {
    return 14.014 * QE ; // eV => J source: NIST Chemical Webbook
  } else if (m_species == "CN") {
    return 13.598 * QE ; // eV => J source: NIST Chemical Webbook
  } else if (m_species == "C") {
    return 11.2603 * QE ; // eV => J source: NIST Chemical Webbook
  } else if (m_species == "H") {
    return EIONH;
  } else if (m_species == "H2") {
    return 1488.0 / 1000.0 /  NA; // kJ/mol -> J/atom
  } else {
    cout << "ERROR: Missing ionization energy for " << m_species << "!" << endl;
    return 0;
  }
 
}
