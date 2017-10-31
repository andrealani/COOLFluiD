#include "MolecularElecPartFunc.h"

#include <limits>
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

MolecularElecPartFunc::MolecularElecPartFunc(const std::string& name)
    : m_name(name), m_ntr(76), m_ntv(76)
{
    loadEnergyLevels();

    mp_tv = new float [m_ntv];
    mp_tr = new float [m_ntv];
    mp_params = new double [m_ntv*m_ntr*m_elec];
    
    loadTemperatureGrid(); 
    loadParameters();
}

MolecularElecPartFunc::~MolecularElecPartFunc()
{
    delete [] mp_tr;
    delete [] mp_tv;
    delete [] mp_params;
    delete [] mp_energy;
}

double MolecularElecPartFunc::Qel(const int i, const double& tr, const double& tv)
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

    const double p11 = mp_params[( itv   *m_ntr+itr  )*m_elec+i];
    const double p12 = mp_params[( itv   *m_ntr+itr+1)*m_elec+i];
    const double p21 = mp_params[((itv+1)*m_ntr+itr  )*m_elec+i];
    const double p22 = mp_params[((itv+1)*m_ntr+itr+1)*m_elec+i];

    return w11*p11 + w12*p12 + w21*p21 + w22*p22;
}

void MolecularElecPartFunc::loadTemperatureGrid()
{
    std::string datadir = std::string(std::getenv("HTGR_DATA_DIRECTORY"));
    std::string database = datadir + "/lbl_data/part/molecules/" + m_name + "_elec.part";
    std::ifstream file(database.c_str());

    // Read the temperature grid, same for Tr and Tv
    double tmp;
    string line;
    for (int i=0; i<m_ntr; i++) {
       file >> tmp >> mp_tr[i];
       file.ignore(numeric_limits<int>::max(), '\n');
       mp_tv[i] = mp_tr[i];
    }
    file.close();
}

void MolecularElecPartFunc::loadParameters()
{
    // Open the file where parameters are stored
    std::string datadir = std::string(std::getenv("HTGR_DATA_DIRECTORY"));
    std::string database = datadir + "/lbl_data/part/molecules/" + m_name + "_elec.part";
    std::ifstream file(database.c_str());
    double temp;
    
    // Read the parameter information
    for (int i = 0; i < m_ntr*m_ntv; ++i) {
        file >> temp;
        file >> temp;
        for (int k = 0; k < m_elec; ++k)
            file >> mp_params[i*m_elec + k];
    }
    
    file.close();
}

void MolecularElecPartFunc::loadEnergyLevels()
{
    // Open the file where energy levels are stored
    std::string datadir = std::string(std::getenv("HTGR_DATA_DIRECTORY"));
    std::string database = datadir + "/lbl_data/part/molecules/" + m_name + ".don";
    std::ifstream file(database.c_str());

    // Read the energy levels
    file >> m_elec;
    mp_energy = new double [m_elec];
    for (int i = 0; i < m_elec; ++i)
        file >> mp_energy[i];

    file.close();
}
