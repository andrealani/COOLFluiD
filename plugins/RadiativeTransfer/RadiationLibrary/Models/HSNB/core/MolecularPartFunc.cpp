#include "MolecularPartFunc.h"

#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <cstdlib>
#include <cmath>

#include "StringUtils.h"

using namespace std;

MolecularPartFunc::MolecularPartFunc(const std::string& name)
    : m_name(name), m_ntr(76), m_ntv(76)
{

    mp_tv = new float [m_ntv];
    mp_tr = new float [m_ntv];
    mp_params = new double [m_ntv*m_ntr];
    
    loadTemperatureGrid(); 
    loadParameters();

}

MolecularPartFunc::~MolecularPartFunc()
{
    delete [] mp_tr;
    delete [] mp_tv;
    delete [] mp_params;
}

double MolecularPartFunc::debye(ThermoData& thermo)
{
    double den = 1.0e-50;

    size_t ie = thermo.speciesIndex("e-");
    if (ie < thermo.nSpecies()) {
        den = thermo[ie].charge;
        den *= den*thermo.N()[ie]*thermo.Th()/thermo.Te();
    }

    for (size_t i = 0; i < thermo.nSpecies(); ++i) {
        if (i == ie) continue;
        const SpeciesData& sp = thermo[i];
        den += sp.charge*sp.charge*thermo.N()[i];
    }

    return std::sqrt(thermo.Th()*EPS0*KB/den)/QE;
    
}

double MolecularPartFunc::Q(const double& tr, const double& tv)
{
    // Step 1: Determine the indices for Tv and Tr which either bound the
    // lower left corner of the the quad in which (Tv,Tr) falls or the
    // correct lower left corner which will be used to extrapolate to (Tv,Tr)
    // if it falls out of the grid boundaries.
    size_t itv=0, itr=0;
    while (tr >= mp_tr[itr] & itr<m_ntr-2) 
        itr++;
    while (tv >= mp_tv[itv] & itv<m_ntv-2) 
        itv++;

    // Step 2: Compute the bilinear interpolation
    const double tr1 = mp_tr[itr];
    const double tr2 = mp_tr[itr+1];
    const double tv1 = mp_tv[itv];
    const double tv2 = mp_tv[itv+1];

    const double area = (tv2-tv1)*(tr2-tr1);
    const double w11 = (tv2-tv)*(tr2-tr)/area;
    const double w12 = (tv2-tv)*(tr-tr1)/area;
    const double w21 = (tv-tv1)*(tr2-tr)/area;
    const double w22 = (tv-tv1)*(tr-tr1)/area;

    const double p11 = mp_params[ itv   *m_ntr+itr  ];
    const double p12 = mp_params[ itv   *m_ntr+itr+1];
    const double p21 = mp_params[(itv+1)*m_ntr+itr  ];
    const double p22 = mp_params[(itv+1)*m_ntr+itr+1];

    return w11*p11 + w12*p12 + w21*p21 + w22*p22;

}

void MolecularPartFunc::loadTemperatureGrid()
{

    std::string datadir = std::string(std::getenv("HTGR_DATA_DIRECTORY"));
    std::string database = datadir + "/lbl_data/part/molecules/" + m_name + ".part";
    std::ifstream file(database.c_str());

    // Read the temperature grid, same for Tr and Tv
    double temp;
    for (int i=0; i<m_ntr; i++) {
       mp_tr[i]=0.0;
       file >> temp;
       file >> mp_tr[i];
       file >> temp;
       mp_tv[i] = mp_tr[i];
    }

   file.close();

}

void MolecularPartFunc::loadParameters()
{
    // Open the file where parameters are stored
    std::string datadir = std::string(std::getenv("HTGR_DATA_DIRECTORY"));
    std::string database = datadir + "/lbl_data/part/molecules/" + m_name + ".part";
    std::ifstream file(database.c_str());
    double temp;
    
    // Read the parameter information
    for (int i = 0; i < m_ntr*m_ntv; ++i) {
        file >> temp;
        file >> temp;
        file >> mp_params[i];
    }
    
    file.close();
}


