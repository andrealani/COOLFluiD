#include "AtomicPartFunc.h"
#include "ThermoData.h"
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cassert>
#include <algorithm>

using namespace std;



AtomicPartFunc::AtomicPartFunc(const std::string& name) :
    m_name(name)
{
    // Open the partition function data file
    std::string datadir = std::string(std::getenv("HTGR_DATA_DIRECTORY"));
    string database = datadir + "/lbl_data/part/atoms/" + m_name + ".part";
    ifstream file(database.c_str());

    if (!file.is_open()) {
        cout << "Could not open database file " << database << "!" << endl;
        exit(1);
    }

    // Load the levels
    int nlevels;
    file >> nlevels;

    assert(nlevels > 0);
    m_levels.resize(nlevels);

    for (int i = 0; i < nlevels; ++i) {
        file >> m_levels[i].E;
        file >> m_levels[i].g;
        file >> m_levels[i].l;
    }

    file.close();

    // Order the levels in terms of the limit
    std::sort(m_levels.begin(), m_levels.end(), CompareLevels());
}

double AtomicPartFunc::debye(ThermoData& thermo)
{
    double th = thermo.Th();
    double te = thermo.Te();
    double* numberDensity = thermo.N();

    double den = 1.0e-50;

    size_t ie = thermo.speciesIndex("e-");
    if (ie < thermo.nSpecies()) {
        den = thermo[ie].charge;
        den *= den*numberDensity[ie]*th/te;
    }

    for (size_t i = 0; i < thermo.nSpecies(); ++i) {
        if (i == ie) continue;
        const SpeciesData& sp = thermo[i];
        den += sp.charge*sp.charge*numberDensity[i];
    }

    return std::sqrt(th*EPS0*KB/den)/QE;
}

double AtomicPartFunc::Q(ThermoData& thermo)
{

    const size_t index = thermo.speciesIndex(m_name);

    // Compute neff and nlim
    const double alpha = thermo[index].charge;
    const double rhod  = debye(thermo);
    const double api   = (alpha+1.0)*QE*QE/(FOURPI*EPS0*rhod);
    const double neff  = std::min((alpha+1.0)*std::sqrt(EIONH/api),100.0);
    const int    nlim  = (int) floor(neff);
    const double fac   = -EIONH/(KB*thermo.Tel());

    // Loop over each level and compute the partition function
    double Q1 = 0.0;
    double Q2 = 0.0;
    
    // Note we could speed this up by finding the index that we would break on
    // before the loop and only loop up to that point, then we could implement
    // this with no branching.  If this is still a bottle-neck, then could
    // consider making a table lookup instead.  Should setup some testing before
    // moving forward to make sure we don't break it.
    double E;
    for (int i = 0; i < m_levels.size(); ++i) {
        if (m_levels[i].l > nlim + 1) break;
        E = m_levels[i].g*std::exp(m_levels[i].E*fac);
        Q2 += E;
        if (m_levels[i].l <= nlim)
            Q1 += E;
    }

    //cout << "Q = " << (Q1 + (neff-(double)nlim)*(Q2-Q1)) << endl;
    return (Q1 + (neff-(double)nlim)*(Q2-Q1));
}

double AtomicPartFunc::Qneg(const double& t)
{

  if (m_name == "N-") {
     return 9.0+5.0*std::exp(-1.333*QE/(KB*t))
               +1.0*std::exp(-2.771*QE/(KB*t));
  } else if (m_name == "O-") {
     return 6.0+6.0*std::exp(-2.000*QE/(KB*t));
  } else if (m_name == "C-") {
     return 4.0+10.*std::exp(-1.23*QE/(KB*t));
  } else {
    cout << "ERROR Q X-" << endl; 
    return 0.0;
  }

}
