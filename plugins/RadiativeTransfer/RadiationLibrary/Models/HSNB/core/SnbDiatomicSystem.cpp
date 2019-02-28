#include "SnbDiatomicSystem.h"
#include "QssMolecules.h"
#include "RadiativeTransfer/RadiationLibrary/Models/HSNB/core/PhotonPath.h"

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
using namespace COOLFluiD;
using namespace COOLFluiD::RadiativeTransfer;

SnbDiatomicSystem::SnbDiatomicSystem(SpeciesLoadData loadData,
				     const std::string& dpdf,
				     const std::string& nup,
				     const ThermoData& thermo)
  : m_directory(loadData.baseDirectory)
{

    m_usePrecomputedParams=true;
    m_paramCount = 0;

    // Add the full path to the default data directory

    m_species = loadData.speciesName;
    m_system  = loadData.systemName;
    m_directory += "/data/modele_msbe/";

    m_directory += m_species + "/" + m_system;


    m_test_qss = false;
    if(m_system[m_system.size()-1] == '*') {
        m_test_qss = true;
        m_system.erase(m_system.size()-1);
        m_directory.erase(m_directory.size()-1);
    }

    m_species_index = thermo.speciesIndex(m_species);

    CFLog(INFO, "SnbDiatomicSystem::SnbDiatomicSystem => m_species=" << this->speciesName() << ", m_system=" << this->m_system
          << ", Index=" << m_species_index<< "\n");

    if (m_species_index < 0) {
        cout << "Error loading SNB system '" << m_system << "'! "
             << "The species isn't loaded for this system." << endl;
        exit(1);
    }

    if (dpdf == "E")
        m_dpdf = EXPONENTIAL;
    else
        m_dpdf = TAILED_INVERSE_EXP;
    
    // Determine the min and max band
    determineBandRange();


    if (thermo.sigMin()!=-1) {
       CFuint lowWav=round(thermo.sigMin() / 1000.0);
       CFuint hiWav=round(thermo.sigMax() / 1000.0);

       if (lowWav<hiWav) {
          if (waveNumberIsInBandRange(lowWav)) {
               m_band1=lowWav;
          }


          if (waveNumberIsInBandRange(hiWav)) {
               m_bandn=hiWav;
          }

          // If the custom spectrum is not contained by the database information skip species
          if (m_band1>hiWav || m_bandn<lowWav) {
               m_band1=m_bandn+1;
          }

          m_nbands = m_bandn - m_band1 + 1;
//          cout <<"Diatomic species: " << speciesName() <<",lowWav " << lowWav << ", hiWav " << hiWav <<", m_band1 " << m_band1
//               << ", m_bandn " << m_bandn << ", m_nbands: " << m_nbands << "\n";
       }



       m_nbands = m_bandn - m_band1 + 1;


    }


    // Check if the band range could be loaded
    if (m_nbands < 0) {
        cout << "Could not find SNB system '" << m_system << "'" << endl;
        exit(1);
    }

    // Load the temperature grid and determine number of parameters
    loadTemperatureGrid();

    if (nup == "T")
        m_nup = THIN;
    else if (nup == "LS")
        m_nup = LINDQUIST_SIMMONS;
    else
        m_nup = CURTIS_GODSON;
    // Force the non uniform path treatment to thin for thin systems
    if(m_nparams == 2) m_nup = THIN;
 
//    cout << "Loading diatomic system " << m_system
//              << " (" << m_band1 << " - " << m_bandn << ") "
//              << m_nparams << ", Doppler PDF = "
//              << (m_dpdf == EXPONENTIAL ? "E" : "H")
//              << ", Non uniform path approximation = "
//              << (m_nup == THIN ? "T" : (m_nup == CURTIS_GODSON ? "CG" : "LS" ))
//              << ", QSS model = "
//              << (m_test_qss == true ? "YES" : "NO") << endl;
    
    // Allocate storage for the parameter information
    mp_params = new double [m_npoints * m_nbands * m_nparams];
    mp_locparams = NULL;
    
    // Now load the parameter information
    #pragma omp parallel for
    for (size_t i = 0; i < m_nbands; ++i)
        loadBandParameters(i);
}

SnbDiatomicSystem::SnbDiatomicSystem(const SnbDiatomicSystem& system)
    : m_directory(system.m_directory),
      m_species(system.m_species),
      m_system(system.m_system),
      m_species_index(system.m_species_index),
      m_dpdf(system.m_dpdf),
      m_nup(system.m_nup),
      m_test_qss(system.m_test_qss),
      m_band1(system.m_band1),
      m_bandn(system.m_bandn),
      m_nbands(system.m_nbands),
      m_nparams(system.m_nparams),
      m_nv(system.m_nv),
      m_npoints(system.m_npoints),
      mp_tv(system.mp_tv == NULL ? NULL : new float [m_nv]),
      mp_tr(system.mp_tr == NULL ? NULL : new float [system.mp_iv[1]]),
      mp_iv(system.mp_iv == NULL ? NULL : new size_t [m_nv+1]),
      mp_params(system.mp_params == NULL ? NULL :
          new double [m_npoints*m_nbands*m_nparams]),
      mp_locparams(NULL)
{
    copy(system.mp_tv, system.mp_tv+m_nv, mp_tv);
    copy(system.mp_tr, system.mp_tr+system.mp_iv[1], mp_tr);
    copy(system.mp_iv, system.mp_iv+m_nv+1, mp_iv);
    copy(
        system.mp_params, system.mp_params+m_npoints*m_nbands*m_nparams,
        mp_params);
}

void swap(SnbDiatomicSystem& s1, SnbDiatomicSystem& s2)
{
    std::swap(s1.m_directory, s2.m_directory);
    std::swap(s1.m_species, s2.m_species);
    std::swap(s1.m_system, s2.m_system);
    std::swap(s1.m_species_index, s2.m_species_index);
    std::swap(s1.m_dpdf, s2.m_dpdf);
    std::swap(s1.m_nup, s2.m_nup);
    std::swap(s1.m_test_qss, s2.m_test_qss);
    std::swap(s1.m_band1, s2.m_band1);
    std::swap(s1.m_bandn, s2.m_bandn);
    std::swap(s1.m_nbands, s2.m_nbands);
    std::swap(s1.m_nparams, s2.m_nparams);
    std::swap(s1.m_nv, s2.m_nv);
    std::swap(s1.m_npoints, s2.m_npoints);
    std::swap(s1.mp_tv, s2.mp_tv);
    std::swap(s1.mp_tr, s2.mp_tr);
    std::swap(s1.mp_iv, s2.mp_iv);
    std::swap(s1.mp_params, s2.mp_params);
    std::swap(s1.mp_locparams, s2.mp_locparams);
}

SnbDiatomicSystem::~SnbDiatomicSystem()
{
    delete [] mp_tr;
    delete [] mp_tv;
    delete [] mp_iv;
    delete [] mp_params;


    if (mp_locparams)
        delete [] mp_locparams;
}

void SnbDiatomicSystem::getParameters(
    double Tv, double Tr, double* const p_params)
{
    // Step 1: Determine the indices for Tv and Tr which either bound the
    // lower left corner of the the quad in which (Tv,Tr) falls or the
    // correct lower left corner which will be used to extrapolate to (Tv,Tr)
    // if it falls out of the grid boundaries.

    // Index for Tv:  0 <= itv < nv-2
    m_pLower = std::lower_bound(mp_tv, mp_tv+m_nv, (float)Tv);

//    size_t m_itv;
    if (m_pLower == mp_tv+m_nv)
        m_itv = m_nv-2;
    else
        m_itv = (mp_tv == m_pLower ? 0 : std::distance(mp_tv, m_pLower)-1);
    
    // Index for Tr:  2 <= itr < min(nr(itv), nr(itv+1), the lower right
    // corner must also be available, otherwise itr is adjusted (itr
    // measures the distance from the maximum, thus the smaller it is, the
    // higher the temperature)
    m_nr0 = mp_iv[m_itv+1]-mp_iv[m_itv];
    m_nr1 = mp_iv[m_itv+2]-mp_iv[m_itv+1];
    m_pMax = mp_tr + mp_iv[1];
    m_pMin = m_pMax - m_nr0;
    m_pLower = std::lower_bound(m_pMin, m_pMax, (float)Tr);

    if (m_pLower == m_pMax)
        m_itr = 2;
    else
        m_itr = std::min(
            m_pLower == m_pMin ? m_nr0 : std::distance(m_pLower, m_pMax)+1, m_nr1);
    
    // Step 2: Compute the bilinear interpolation
    m_tv1 = mp_tv[m_itv];
    m_tv2 = mp_tv[m_itv+1];
    m_tr1 = mp_tr[mp_iv[1]-m_itr];
    m_tr2 = mp_tr[mp_iv[1]-m_itr+1];
    
//    cout << m_system << endl;
//    cout << tv1 << ", " << tv2 << endl;
//    cout << tr1 << ", " << tr2 << endl;
//    cout << Tv  << ", " << Tr  << endl;
//    cout << itv << ", " << itr << endl;

    m_area = (m_tv2-m_tv1)*(m_tr2-m_tr1);
//    cout << "area = " << area << endl;
    m_w11 = (m_tv2-Tv)*(m_tr2-Tr)/m_area;
    m_w12 = (m_tv2-Tv)*(Tr-m_tr1)/m_area;
    m_w21 = (Tv-m_tv1)*(m_tr2-Tr)/m_area;
    m_w22 = (Tv-m_tv1)*(Tr-m_tr1)/m_area;
    
    m_ndata = m_nparams*m_nbands;
    m_p11 = mp_params + (mp_iv[m_itv]+(m_nr0-m_itr))*m_ndata;
    m_p12 = m_p11 + m_ndata;
    m_p21 = mp_params + (mp_iv[m_itv+1]+(m_nr1-m_itr))*m_ndata;
    m_p22 = m_p21 + m_ndata;

    for (int i = 0; i < m_ndata; ++i) {
        p_params[i] = m_w11*m_p11[i] + m_w12*m_p12[i] + m_w21*m_p21[i] + m_w22*m_p22[i];
        if (p_params[i]<0.0)
            p_params[i]=0.0;
    }
}



void SnbDiatomicSystem::getThickParamsSingleBand(ThermoData &thermo,CFuint cellID, double band)
{
    thermo.setState(cellID);
    mp_fac=1.0;
//    m_nparams=4;

    m_tr  = thermo.Tr();
    m_tv  = thermo.Tv();

    m_pa  = thermo.N(m_species_index) * KB * m_tr;
    m_pLower = std::lower_bound(mp_tv, mp_tv+m_nv, (float)m_tv);

    if (m_pLower == mp_tv+m_nv)
        m_itv = m_nv-2;
    else
        m_itv = (mp_tv == m_pLower ? 0 : std::distance(mp_tv, m_pLower)-1);

    m_nr0 = mp_iv[m_itv+1]-mp_iv[m_itv];
    m_nr1 = mp_iv[m_itv+2]-mp_iv[m_itv+1];

    m_pMax = mp_tr + mp_iv[1];
    m_pMin = m_pMax - m_nr0;
    m_pLower = std::lower_bound(m_pMin, m_pMax, (float)m_tr);

    if (m_pLower == m_pMax)
        m_itr = 2;
    else
        m_itr = std::min(
            m_pLower == m_pMin ? m_nr0 : std::distance(m_pLower, m_pMax)+1, m_nr1);

    // Step 2: Compute the bilinear interpolation
    m_tv1 = mp_tv[m_itv];
    m_tv2 = mp_tv[m_itv+1];
    m_tr1 = mp_tr[mp_iv[1]-m_itr];
    m_tr2 = mp_tr[mp_iv[1]-m_itr+1];

    m_area = (m_tv2-m_tv1)*(m_tr2-m_tr1);
//    cout << "area = " << area << endl;
    m_w11 = (m_tv2-m_tv)*(m_tr2-m_tr)/m_area;
    m_w12 = (m_tv2-m_tv)*(m_tr-m_tr1)/m_area;
    m_w21 = (m_tv-m_tv1)*(m_tr2-m_tr)/m_area;
    m_w22 = (m_tv-m_tv1)*(m_tr-m_tr1)/m_area;

    m_ndata = m_nparams*m_nbands;


    m_p11 = mp_params + (mp_iv[m_itv]+(m_nr0-m_itr))*m_ndata;
    m_p12 = m_p11 + m_ndata;
    m_p21 = mp_params + (mp_iv[m_itv+1]+(m_nr1-m_itr))*m_ndata;
    m_p22 = m_p21 + m_ndata;

    m_currentParamIndex = m_nparams*band;

    //Has to be multiplied with cellDistance
    m_tempKu = (m_w11*m_p11[m_currentParamIndex] + m_w12*m_p12[m_currentParamIndex] + m_w21*m_p21[m_currentParamIndex] + m_w22*m_p22[m_currentParamIndex]) * m_pa;

    //        mp_locparams[i*ndata+b*4+1] = w11*p11[i+2] + w12*p12[i+2] + w21*p21[i+2] + w22*p22[i+2] * mp_fac;
    m_tempBetaD = m_w11*m_p11[m_currentParamIndex+3] + m_w12*m_p12[m_currentParamIndex+3] + m_w21*m_p21[m_currentParamIndex+3] + m_w22*m_p22[m_currentParamIndex+3];



    if (m_nparams < 6)
        m_tempBetaL = (m_w11*m_p11[m_currentParamIndex+4] + m_w12*m_p12[m_currentParamIndex+4] + m_w21*m_p21[m_currentParamIndex+4] + m_w22*m_p22[m_currentParamIndex+4])*m_pa/100000.0;
    else
        m_tempBetaL = (m_w11*m_p11[m_currentParamIndex+4] + m_w12*m_p12[m_currentParamIndex+4] + m_w21*m_p21[m_currentParamIndex+4] + m_w22*m_p22[m_currentParamIndex+4])*m_pa
                + m_w11*m_p11[m_currentParamIndex+5] + m_w12*m_p12[m_currentParamIndex+5] + m_w21*m_p21[m_currentParamIndex+5] + m_w22*m_p22[m_currentParamIndex+5];


    if (m_tempKu<0.0) {
        m_tempKu=0.0;
    }

    if (m_tempBetaD<0.0) {
        m_tempBetaD=0.0;
    }

    if (m_tempBetaL<0.0) {
        m_tempBetaL=0.0;
    }

}

void SnbDiatomicSystem::getThinParamsSingleBand(ThermoData &thermo,CFuint cellID, double band)
{
    thermo.setState(cellID);
//    m_nparams=2;

    m_tr  = thermo.Tr();
    m_tv  = thermo.Tv();
    m_pa  = thermo.N(m_species_index) * KB * m_tr;

    // Step 1: Determine the indices for Tv and Tr which either bound the
    // lower left corner of the the quad in which (Tv,Tr) falls or the
    // correct lower left corner which will be used to extrapolate to (Tv,Tr)
    // if it falls out of the grid boundaries.

    // Index for Tv:  0 <= itv < nv-2
    m_pLower = std::lower_bound(mp_tv, mp_tv+m_nv, (float)m_tv);

//    size_t m_itv;
    if (m_pLower == mp_tv+m_nv)
        m_itv = m_nv-2;
    else
        m_itv = (mp_tv == m_pLower ? 0 : std::distance(mp_tv, m_pLower)-1);

    // Index for Tr:  2 <= itr < min(nr(itv), nr(itv+1), the lower right
    // corner must also be available, otherwise itr is adjusted (itr
    // measures the distance from the maximum, thus the smaller it is, the
    // higher the temperature)
    m_nr0 = mp_iv[m_itv+1]-mp_iv[m_itv];
    m_nr1 = mp_iv[m_itv+2]-mp_iv[m_itv+1];
    m_pMax = mp_tr + mp_iv[1];
    m_pMin = m_pMax - m_nr0;
    m_pLower = std::lower_bound(m_pMin, m_pMax, (float)m_tr);

    if (m_pLower == m_pMax)
        m_itr = 2;
    else
        m_itr = std::min(
            m_pLower == m_pMin ? m_nr0 : std::distance(m_pLower, m_pMax)+1, m_nr1);

    // Step 2: Compute the bilinear interpolation
    m_tv1 = mp_tv[m_itv];
    m_tv2 = mp_tv[m_itv+1];
    m_tr1 = mp_tr[mp_iv[1]-m_itr];
    m_tr2 = mp_tr[mp_iv[1]-m_itr+1];


    m_area = (m_tv2-m_tv1)*(m_tr2-m_tr1);
//    cout << "area = " << area << endl;
    m_w11 = (m_tv2-m_tv)*(m_tr2-m_tr)/m_area;
    m_w12 = (m_tv2-m_tv)*(m_tr-m_tr1)/m_area;
    m_w21 = (m_tv-m_tv1)*(m_tr2-m_tr)/m_area;
    m_w22 = (m_tv-m_tv1)*(m_tr-m_tr1)/m_area;

    m_ndata = m_nparams*m_nbands;
    m_p11 = mp_params + (mp_iv[m_itv]+(m_nr0-m_itr))*m_ndata;
    m_p12 = m_p11 + m_ndata;
    m_p21 = mp_params + (mp_iv[m_itv+1]+(m_nr1-m_itr))*m_ndata;
    m_p22 = m_p21 + m_ndata;


    //
    m_currentParamIndex = m_nparams*band;
    //Still has to be multiplied with cellDistance!
    m_tempKu = (m_w11*m_p11[m_currentParamIndex] + m_w12*m_p12[m_currentParamIndex] + m_w21*m_p21[m_currentParamIndex] + m_w22*m_p22[m_currentParamIndex])*m_pa;


    if (m_tempKu<0.0) {
        m_tempKu=0.0;
    }

}



CFreal SnbDiatomicSystem::getKappa(ThermoData &thermo, CFuint cellID, double sig)
{

    int b = int(sig / 1000.0);

    // Convert band to local indexing
    if (b < lowBand() || b > highBand()) {
        return 0.0;
    }
    b -= lowBand();

    if (thermo.usePrecomputedDiatomicParameters()) {
        m_tempKu=getLocalParameter(cellID, b, 0);
    }
    else {

        if (m_nparams < 4) { // Local parameters for thin systems
            getThinParamsSingleBand(thermo,cellID, b);
        }
        else {
            getThickParamsSingleBand(thermo,cellID,b);
        }
    }

    return m_tempKu;
}

double SnbDiatomicSystem::getPrecomputedBufferByteSize() const
{
    return m_paramCount*8;
}

NonUniformPath SnbDiatomicSystem::getMechanismType() const
{
    return m_nup;
}

CFreal SnbDiatomicSystem::getBand(CFreal sig)
{
    // Find the band corresponding to this wavelength
    int b = int(sig / 1000.0);

    // Convert band to local indexing
    if (b < lowBand() || b > highBand())
        return 0.0;
    b -= lowBand();

    return b;
}

bool SnbDiatomicSystem::bandEmissionComputed(int forCellID)
{
    return (m_currentBandEmissionCellIndex==forCellID);
}


void SnbDiatomicSystem::setupLocalParameters(ThermoData& thermo)
{
    if (thermo.usePrecomputedDiatomicParameters()) {
        localParamsSetup=true;

        double pa, tr, tv, fac = 1.0;
        double *p_params = new double [m_nparams*m_nbands];
        double* p_fac = new double[thermo.nCells()];
        fill(p_fac, p_fac+thermo.nCells(), 1.0);

        if (m_test_qss) {
            QssMolecules qss(m_system, m_species);
            for (int i=0; i<thermo.nCells(); i++) {
                thermo.setState(i);
                p_fac[i] = qss.levelCorrection(thermo);
                cout << "i = " << i << ", QSS level correct = " << p_fac[i] << endl;
            }
        }

        if (mp_locparams) {
            delete [] mp_locparams;
            mp_locparams = NULL;
        }


        if (m_nparams < 4) { // Local parameters for thin systems
            const size_t ndata = m_nbands*2;

            mp_locparams = new double [thermo.nCells()*ndata];
            m_paramCount=thermo.nCells()*ndata;
//            std::cout << ndata << std::endl;

            for (int i=0; i<thermo.nCells(); i++) {

                thermo.setState(i);

                tr  = thermo.Tr();
                tv  = thermo.Tv();
                pa  = thermo.N(m_species_index) * KB * tr;

                // std::cout << "m_species_index "<< m_species_index << " PARAMS tr:" << tr << " tv: " << tv << "pa: " << pa << " thermo.N: " << thermo.N(m_species_index) << " KB: " << KB  << std::endl;

                getParameters(tv,tr,p_params);

#pragma omp parallel for
                for (int b = 0; b < m_nbands; ++b) {
                    mp_locparams[i*ndata+b*2+0] = p_params[b*m_nparams+0] * pa;
                    mp_locparams[i*ndata+b*2+1] = p_params[b*m_nparams+1] * pa * p_fac[i];
                }
            }

        } else { // Local parameters for thick systems
            const size_t ndata = m_nbands*4;
            mp_locparams = new double [thermo.nCells()*ndata];
            m_paramCount=thermo.nCells()*ndata;
//            std::cout << ndata << std::endl;

            for (int i=0; i<thermo.nCells(); i++) {

                thermo.setState(i);

                tr  = thermo.Tr();
                tv  = thermo.Tv();
                pa  = thermo.N(m_species_index) * KB * tr;

                getParameters(tv,tr,p_params);

#pragma omp parallel for
                for (int b = 0; b < m_nbands; ++b) {
                    mp_locparams[i*ndata+b*4+0] = p_params[b*m_nparams+0] * pa;
                    mp_locparams[i*ndata+b*4+1] = p_params[b*m_nparams+2] * p_fac[i];
                    mp_locparams[i*ndata+b*4+2] = p_params[b*m_nparams+3] ;
                    if (m_nparams < 6)
                        mp_locparams[i*ndata+b*4+3] = p_params[b*m_nparams+4]*pa/100000.0;
                    else
                        mp_locparams[i*ndata+b*4+3] = p_params[b*m_nparams+4]*pa
                                + p_params[b*m_nparams+5];
                }
            }
        }

        delete [] p_params;
        delete [] p_fac;
    }
}

double SnbDiatomicSystem::getLocalParameter(const int& i, const int& j, const int& k) const
// Cell (i), Band (j) and Parameter (k) indices
{
    //Careful if mp_locparams is not intialised!

    if (m_nparams < 4) // Local parameters for thin systems
        return mp_locparams[i*m_nbands*2+j*2+k];
    else 
        return mp_locparams[i*m_nbands*4+j*4+k];
    
}

double SnbDiatomicSystem::emittedPower(int cellID, ThermoData& thermo)
{
    double *p_params = new double [m_nparams*m_nbands];

    //Not necessary: state is set in the HSNBRadiator setState function
    //thermo.setState(cellID);

    double tr  = thermo.Tr();
    double tv  = thermo.Tv();
    double pa  = thermo.N(m_species_index) * KB * tr;

//    CFLog(VERBOSE, "SnbDiatomicSystem::emittedPower => "<< this->speciesName() << ", " << this->systemName() <<  ", State: N("<<m_species_index << ")=" << thermo.N(m_species_index) << ", KB=" << KB << ", tr=" << tr << "\n");
//    CFLog(VERBOSE, "SnbDiatomicSystem::emittedPower => "<< this->systemName() <<  ", State: Tr=" << tr << ", Tv=" << tv << ", pa=" << pa << "\n");
    getParameters(tv, tr, p_params);

    // QSS level correction
    double qss_fac = 1.0;
    if (m_test_qss) {
        QssMolecules qss(m_system, m_species);
        //thermo.setState(cellID);
        qss_fac = qss.levelCorrection(thermo);
    }

    double sum = 0;
    for (int b = 0; b < m_nbands; ++b)
        sum += p_params[b*m_nparams+1];

//    CFLog(VERBOSE, "SnbDiatomicSystem::emittedPower => " << this->systemName() << " emission coeff=" <<  sum * pa * qss_fac * 1000.0 << "\n");

    delete[] p_params;
    return (sum * pa * qss_fac * 1000.0);
}

void SnbDiatomicSystem::bandEmission(const int cellID, ThermoData& thermo, double* const p_emis)
{
    if (bandEmissionComputed(cellID)==false) {
    double *p_params = new double [m_nparams*m_nbands];
    thermo.setState(cellID);

    m_tr  = thermo.Tr();
    m_tv  = thermo.Tv();
    m_pa  = thermo.N(m_species_index) * KB * m_tr;

    getParameters(m_tv, m_tr, p_params);

    // QSS level correction
    double qss_fac = 1.0;
    if (m_test_qss) {
        QssMolecules qss(m_system, m_species);
        thermo.setState(cellID);
        qss_fac = qss.levelCorrection(thermo);
    }

    double sum = 0;
    for (int b = 0; b < m_nbands; ++b) {
        p_emis[b] = p_params[b*m_nparams+1] * m_pa * qss_fac * 1000.0;
    }

    delete[] p_params;
    }
}

void SnbDiatomicSystem::determineBandRange()
{
    // Lower bound, start from 1000
    m_bandn = 1;
    m_band1 = 200;
    for (size_t band = 1; band < 200; ++band) {
        if (std::ifstream((m_directory+"/"+bandFilename(band)).c_str())) {
            m_band1 = band;
            break;
        }
    }
    
    // Upper bound start from 200,000
    for (size_t band = 200; band > m_band1; --band) {
        if (std::ifstream((m_directory+"/"+bandFilename(band)).c_str())) {
            m_bandn = band;
            break;
        }
    }
    
    m_nbands = m_bandn - m_band1 + 1;
}

bool SnbDiatomicSystem::waveNumberIsInBandRange(CFuint waveBand)
{
    return ((waveBand>=lowBand()) && (waveBand<=highBand()));
}

void SnbDiatomicSystem::loadTemperatureGrid()
{
    // Open the file corresponding to the first band
    std::ifstream file((m_directory+"/"+bandFilename(m_band1)).c_str());
    
    // Skip first two lines
    std::string line;
    std::getline(file, line);
    std::getline(file, line);
    
    // Now read the first temperatures from the file and count the number of
    // parameters (extra columns)
    std::vector<float>  tv(2);
    std::vector<float>  tr(2);
    std::vector<size_t> ntr(1, 1);
    
    std::getline(file, line);
    std::stringstream ss(line);
    ss >> tv[0];
    ss >> tr[0];
    
    m_nparams = 0;
    while (ss >> line)
        m_nparams++;
    


    // Read the rest of the temperatures from the file
    while (std::getline(file, line)) {
        sscanf(line.c_str(), "%f %f", &tv.back(), &tr.back());
        if (tr[tr.size()-1] < tr[tr.size()-2]) {
            ntr.push_back(0);
            tv.push_back(0.0f);
        }
        ntr.back()++;
        tr.push_back(0.0f);
    }
    
    tv.pop_back();
    tr.pop_back();
    file.close();
    
    // Compute sizes
    m_nv = tv.size();
    m_npoints = tr.size();
    
    // Allocate storage for temperatures
    mp_tv = new float [m_nv];
    mp_iv = new size_t [m_nv+1];
    mp_tr = new float [ntr[0]];
    
    // Copy the temperatures to their storage
    std::copy(tv.begin(), tv.end(), mp_tv);
    std::copy(tr.begin(), tr.begin()+ntr[0], mp_tr);
    
    // Store the indices to the first point in each Tv group
    mp_iv[0] = 0;
    for (int i = 0; i < m_nv; ++i)
        mp_iv[i+1] = mp_iv[i] + ntr[i];
}


void SnbDiatomicSystem::loadBandParameters(const size_t& iband)
{
    // Open the file corresponding to the band
    std::ifstream file((m_directory+"/"+bandFilename(iband+m_band1)).c_str());
    
    // Skip the first two lines
    std::string line;
    std::getline(file, line);
    std::getline(file, line);
    
    // Read the parameter information
    for (size_t ipoint = 0; ipoint < m_npoints; ++ipoint) {
        // Get pointer to storage location for this point/band
        double* const p_params =
            mp_params + m_nparams * (iband + m_nbands * ipoint);
        // Skip temperatures
        file >> line;
        file >> line;
        // Read parameters
        for (size_t k = 0; k < m_nparams; ++k) {
            file >> line;
            // Take care of truncated exponents in scientific notation
            size_t pos = line.find_first_of("eDd", 1);
            if (pos != string::npos)
                line = line.replace(pos, 1, "E");
            pos = line.find_first_of("+-", 1);
            if (pos != string::npos && line[pos-1] != 'E')
                line = line.insert(pos, "E");
            p_params[k] = atof(line.c_str());
        }
    }
    
    // Close the file
    file.close();
}




double SnbDiatomicSystem::tau(const HSNBNonThickParameterSet &pathParams)
{
    return std::exp(-opticalThickness(pathParams));
}

double SnbDiatomicSystem::tau(const HSNBThickParameterSet &pathParams)
{
    return std::exp(-opticalThickness(pathParams));
}

double SnbDiatomicSystem::tauShortened(const HSNBThickParameterSet &pathParams, CFreal kappa0, CFreal distance0)
{
    return std::exp(-opticalThicknessShortened(pathParams,kappa0,distance0));
}





double SnbDiatomicSystem::opticalThickness(const HSNBNonThickParameterSet &pathParams)
{
    CFLog(DEBUG_MAX, "SnbDiatomicSystem::opticalThickness::Thin => kappa=" <<  pathParams.kappa << " \n");
    //pathParams.print();
    return pathParams.kappa;
}

double SnbDiatomicSystem::opticalThicknessShortened(const HSNBThickParameterSet &pathParams, CFreal kappa0, CFreal distance0)
{
    CFLog(DEBUG_MAX, "SnbDiatomicSystem::opticalThicknessShortened::Thick => Start. pathParams.nbCells=" << pathParams.nbCells << "\n");
    //pathParams.print();

    double p_kb [5] = {0, 0, 0, 0, 0};

    double ku, bl, bd;


    switch (m_nup) {

    case CURTIS_GODSON:


        for (int i = 0;  i < pathParams.nbCells; i++) {

            if (i==0) {
                ku = kappa0*distance0;
            }
            else {
                ku = pathParams.kappa[i];
            }

            p_kb[0] += ku;
            p_kb[1] += pathParams.betaD[i]*ku;
            p_kb[2] += pathParams.betaL[i]*ku;

//            CFLog(VERBOSE, "SnbDiatomicSystem::opticalThicknessShortened::Thick, i=" << i << " , CG => kappa=" << pathParams.kappa[i]
//                  << ", betaD=" << pathParams.betaD[i] << ", betaL=" << pathParams.betaL[i] << " \n");

            cf_assert(pathParams.kappa[i] >= 0.0);
            cf_assert(pathParams.betaD[i] >= 0.0);
            cf_assert(pathParams.betaL[i] >= 0.0);
        }

        ku = p_kb[0];
        if (ku == 0.0)
            return 0.0;
        bd = p_kb[1]/ku; //wd = (bd > 1.0e-10 ? wdod(ku, bd) : 0.0);
        bl = p_kb[2]/ku; //wl = (bl > 1.0e-10 ? wlod(ku, bl) : 0.0);
        if(bd == 0.0) bd =1e-20;
        if(bl == 0.0) bl =1e-20;

        //std::cout << systemName() << " " << b << " " << ku << " " << bd << " " << bl << std::endl;
        return wvod(ku, wdod(ku, bd), wlod(ku, bl)); // Voigt profile

    case LINDQUIST_SIMMONS:
        for (int i = 0; i < pathParams.nbCells; i++) {
            // Get the cell number and length associated to the path index

//            CFLog(VERBOSE, "SnbDiatomicSystem::opticalThicknessShortened::Thick, i="<< i <<", LS => kappa=" << pathParams.kappa[i]
//                  << ", betaD=" << pathParams.betaD[i] << ", betaL=" << pathParams.betaL[i] << " \n");

            if (i==0) {
                ku = kappa0*distance0;
            }
            else {
                ku = pathParams.kappa[i];
            }

            assert(ku >= 0.0);
            if (i == 0) {
                bd = pathParams.betaD[i];
                bl = pathParams.betaL[i];
                if (bl == 0.0) bl = 1e-20;
                if (bd == 0.0) bd = 1e-20;
                p_kb[0] = ku;
                p_kb[1] = bd*ku;
                p_kb[2] = bl*ku;
                if (ku == 0.0) {
                    p_kb[3] = 0.0;
                    p_kb[4] = 0.0;
                } else {
                    p_kb[3] = wdod(ku, bd);
                    p_kb[4] = wlod(ku, bl);
                }

            } else if (ku > 0.0) {
                bd = pathParams.betaD[i];
                bl = pathParams.betaL[i];
                if (bd == 0.0) bd = 1e-20;
                if (bl == 0.0) bl = 1e-20;
                wquad(ku, bd, bl, p_kb);
            }

        }
        return wvod(p_kb[0], p_kb[3], p_kb[4]);


    default:
        CFLog(INFO, "SnbDiatomicSystem::opticalThicknessShortened:: Error! No NUP defined? \n");
        return 1.0;
    }
}

double SnbDiatomicSystem::opticalThickness(const HSNBThickParameterSet &pathParams)
{
//    CFLog(VERBOSE, "SnbDiatomicSystem::opticalThickness::Thick, "<<this->speciesName() <<" => Start. pathParams.nbCells=" << pathParams.nbCells<< "\n");
//    pathParams.print();


    double p_kb [5] = {0, 0, 0, 0, 0};

    double ku, bl, bd;


    switch (m_nup) {

    case CURTIS_GODSON:
        for (int i = 0;  i < pathParams.nbCells; i++) {
            ku = pathParams.kappa[i];
            p_kb[0] += ku;
            p_kb[1] += pathParams.betaD[i]*ku;
            p_kb[2] += pathParams.betaL[i]*ku;

//            CFLog(VERBOSE, "SnbDiatomicSystem::opticalThickness::Thick, CG => kappa=" << pathParams.kappa[i]
//                  << ", p_kb[1]=" << p_kb[1] << ", p_kb[2]=" << p_kb[2] << " \n");

//            std::cout << this->systemName() << " / " << this->speciesName() << "SnbDiatomicSystem::opticalThickness::Thick, CG => kappa=" << pathParams.kappa[i]
//                         << ", p_kb[1]=" << p_kb[1] << ", p_kb[2]=" << p_kb[2] << " i=" << i << std::endl;

            cf_assert(pathParams.kappa[i] >= 0.0);
            cf_assert(pathParams.betaD[i] >= 0.0);
            cf_assert(pathParams.betaL[i] >= 0.0);
        }

        ku = p_kb[0];
        if (ku == 0.0)
            return 0.0;
        bd = p_kb[1]/ku; //wd = (bd > 1.0e-10 ? wdod(ku, bd) : 0.0);
        bl = p_kb[2]/ku; //wl = (bl > 1.0e-10 ? wlod(ku, bl) : 0.0);
        if(bd == 0.0) bd =1e-20;
        if(bl == 0.0) bl =1e-20;

        //std::cout << systemName() << " " << b << " " << ku << " " << bd << " " << bl << std::endl;
        return wvod(ku, wdod(ku, bd), wlod(ku, bl)); // Voigt profile

    case LINDQUIST_SIMMONS:
        for (int i = 0; i < pathParams.nbCells; i++) {
            // Get the cell number and length associated to the path index

//            CFLog(VERBOSE, "SnbDiatomicSystem::opticalThickness::Thick, LS => kappa=" << pathParams.kappa[i]
//                  << ", betaD=" << pathParams.betaD[i] << ", betaL=" << pathParams.betaL[i] << " \n");

            ku = pathParams.kappa[i];
            assert(ku >= 0.0);
            if (i == 0) {
                bd = pathParams.betaD[i];
                bl = pathParams.betaL[i];
                if (bl == 0.0) bl = 1e-20;
                if (bd == 0.0) bd = 1e-20;
                p_kb[0] = ku;
                p_kb[1] = bd*ku;
                p_kb[2] = bl*ku;
                if (ku == 0.0) {
                    p_kb[3] = 0.0;
                    p_kb[4] = 0.0;
                } else {
                    p_kb[3] = wdod(ku, bd);
                    p_kb[4] = wlod(ku, bl);
                }

            } else if (ku > 0.0) {
                bd = pathParams.betaD[i];
                bl = pathParams.betaL[i];
                if (bd == 0.0) bd = 1e-20;
                if (bl == 0.0) bl = 1e-20;
                wquad(ku, bd, bl, p_kb);
            }

//            CFLog(VERBOSE, "SnbDiatomicSystem::opticalThickness::Thick, LS => kappa=" << pathParams.kappa[i]
 //                 << ", p_kb[1]=" << p_kb[1] << ", p_kb[2]=" << p_kb[2] << " \n");

        }
        return wvod(p_kb[0], p_kb[3], p_kb[4]);


    default:
        CFLog(INFO, "SnbDiatomicSystem::opticalThickness:: Error! No NUP defined? \n");

    }
}



void SnbDiatomicSystem::addStateParams(ThermoData& thermo, HSNBNonThickParameterSet &nonThickParams, CFreal cellDistance, CFuint localCellID, CFreal sig)
{
//    CFLog(VERBOSE, "SnbDiatomicSystem::addStateParams::NONTHICK, " << this->speciesName() << " / " << this->systemName() << "  \n");

    // Find the band corresponding to this wavelength
    int b = int(sig / 1000.0);


    // Convert band to local indexing
    if (b < lowBand() || b > highBand()) {
//        CFLog(VERBOSE, "SnbDiatomicSystem::addStateParams =>NONTHICK, " << this->speciesName()<<", b=" << b << " out of bounds. \n");
        nonThickParams.addState(0.0);
        return;
    }

    b -= lowBand();

    if (nonThickParams.isEmpty()) {
        if  (thermo.usePrecomputedDiatomicParameters()) {
            m_tempKu=mp_locparams[localCellID*2*m_nbands + b*2] * cellDistance;
        }
        else {
            getThinParamsSingleBand(thermo,localCellID,b);
            m_tempKu*=cellDistance;
        }

//        CFLog(VERBOSE, "SnbDiatomicSystem::addStateParams =>NONTHICK: Empty ParamSet. Add Ku=" << m_tempKu << " at localID " << localCellID  << "\n");
//        CFLog(VERBOSE, "SnbDiatomicSystem::addStateParams =>NONTHICK: mp_locparams["<< localCellID <<"*2*"<< m_nbands<<" + "<< b << "*2]=" << mp_locparams[localCellID*2*m_nbands + b*2] << "\n");
//        CFLog(VERBOSE, "SnbDiatomicSystem::addStateParams =>NONTHICK: cellDistance=" << cellDistance << "\n");

        nonThickParams.addState(m_tempKu);
    }
    else {
        //ku=sum_i<=j(ku)

        if  (thermo.usePrecomputedDiatomicParameters()) {
            m_tempKu=nonThickParams.kappa+mp_locparams[localCellID*2*m_nbands + b*2] * cellDistance;
        }
        else {
            getThinParamsSingleBand(thermo,localCellID,b);
            m_tempKu*=cellDistance;
            m_tempKu+=nonThickParams.kappa;
        }

//        CFLog(VERBOSE, "SnbDiatomicSystem::addStateParams =>NONTHICK: NonEmpty ParamSet. Add Ku=" << m_tempKu << " at localID " << localCellID  << "\n");
//        CFLog(VERBOSE, "SnbDiatomicSystem::addStateParams =>NONTHICK: mp_locparams["<< localCellID <<"*2*"<< m_nbands<<" + "<< b << "*2]=" << mp_locparams[localCellID*2*m_nbands + b*2] << "\n");
//        CFLog(VERBOSE, "SnbDiatomicSystem::addStateParams =>NONTHICK: cellDistance=" << cellDistance << "\n");

        nonThickParams.addState(m_tempKu);
    }

//    CFLog(INFO, "SnbDiatomicSystem::addStateParams =>THIN: cellDistance=" << cellDistance << " m_tempKu=" <<m_tempKu << "\n");

}

void SnbDiatomicSystem::addStateParams(ThermoData &thermo,HSNBThickParameterSet &thickParams, CFreal cellDistance, CFuint localCellID, CFreal sig)
{
//    CFLog(VERBOSE, "SnbDiatomicSystem::addStateParams::THICK, " << this->speciesName() << " / " << this->systemName()<< "  \n");
    // Find the band corresponding to this wavelength
    int b = int(sig / 1000.0);

    // Convert band to local indexing
    if (b < lowBand() || b > highBand()) {
//        CFLog(VERBOSE, "SnbDiatomicSystem::addStateParams =>THICK, " << this->speciesName()<<", b=" << b << " out of bounds. \n");
        thickParams.addState(0.0,0.0,0.0);
        return;
    }
    b -= lowBand();

   if (thermo.usePrecomputedDiatomicParameters()) {
       m_tempKu=mp_locparams[localCellID*4*m_nbands+b*4+0] * cellDistance;
       m_tempBetaD=mp_locparams[localCellID*4*m_nbands+b*4+2];
       m_tempBetaL=mp_locparams[localCellID*4*m_nbands+b*4+3];
   }
   else {
       //TODO use proper function header
       getThickParamsSingleBand(thermo,localCellID,b);
       m_tempKu*=cellDistance;
   }

   cf_assert(m_tempKu >= 0.0);
   cf_assert(m_tempBetaL >= 0.0);
   cf_assert(m_tempBetaD >= 0.0);

   thickParams.addState(m_tempKu,m_tempBetaD, m_tempBetaL);

//   CFLog(VERBOSE, "SnbDiatomicSystem::addStateParams =>THICK: mp_locparams[" << localCellID << "*4*" << m_nbands << " + " << b << "*4+0]=" << mp_locparams[localCellID*4*m_nbands+b*4+0] << "\n");
//   CFLog(VERBOSE, "SnbDiatomicSystem::addStateParams =>THICK: mp_locparams[" << localCellID << "*4*" << m_nbands << " + " << b << "*4+2]=" << mp_locparams[localCellID*4*m_nbands+b*4+2] << "\n");
//   CFLog(VERBOSE, "SnbDiatomicSystem::addStateParams =>THICK: mp_locparams[" << localCellID << "*4*" << m_nbands << " + " << b << "*4+3]=" << mp_locparams[localCellID*4*m_nbands+b*4+3] << "\n");
//   CFLog(VERBOSE, "SnbDiatomicSystem::addStateParams =>THICK: Add Ku=" << m_tempKu <<", betaD=" << m_tempBetaD << ", betaL=" << m_tempBetaL << " at localID " << localCellID  << "\n");
//   CFLog(VERBOSE, "SnbDiatomicSystem::addStateParams =>THICK: cellDistance=" << cellDistance << "\n");
//   CFLog(INFO, "SnbDiatomicSystem::addStateParams =>THICK: cellDistance=" << cellDistance << " m_tempKu=" <<m_tempKu << " m_tempBetaD="<<m_tempBetaD << " m_tempBetaL="<<m_tempBetaL << "\n");

}


double hquad(double x)
{
    double expxi2 [32] = {
        7.5943201E-49, 2.9875169E-43, 6.9915249E-39, 3.4952039E-35,
        6.1363102E-32, 4.9064041E-29, 2.0921338E-26, 5.2846712E-24,
        8.5189081E-22, 9.2610320E-20, 7.0832621E-18, 3.9408644E-16,
        1.6382810E-14, 5.2018952E-13, 1.2847697E-11, 2.5062323E-10,
        3.9117887E-09, 4.9395557E-08, 5.0942322E-07, 4.3261255E-06,
        3.0466224E-05, 1.7901699E-04, 8.8231820E-04, 3.6643531E-03,
        1.2874141E-02, 3.8392534E-02, 9.7457619E-02, 2.1108208E-01,
        3.9082976E-01, 6.1955334E-01, 8.4179756E-01, 9.8105426E-01,
    };
    
    double weights [32] = {
        4.1125314E-01, 3.1721873E-01, 2.7607288E-01, 2.5139122E-01,
        2.3442298E-01, 2.2182848E-01, 2.1201056E-01, 2.0409313E-01,
        1.9754880E-01, 1.9203791E-01, 1.8733019E-01, 1.8326318E-01,
        1.7971845E-01, 1.7660734E-01, 1.7386204E-01, 1.7142969E-01,
        1.6926843E-01, 1.6734470E-01, 1.6563126E-01, 1.6410584E-01,
        1.6275006E-01, 1.6154868E-01, 1.6048903E-01, 1.5956052E-01,
        1.5875434E-01, 1.5806318E-01, 1.5748099E-01, 1.5700284E-01,
        1.5662481E-01, 1.5634386E-01, 1.5615777E-01, 1.5606509E-01,
    };
    
    double ret = 0.0;
    for (size_t i = 0; i < 32; ++i)
        ret += std::log(1.0 + x*expxi2[i]) * weights[i];
    ret *= 2.0; // Symetrical function
    return ret;
}

double equad(double x)
{
    double expxi2 [32] = {
        7.5943201E-49, 2.9875169E-43, 6.9915249E-39, 3.4952039E-35,
        6.1363102E-32, 4.9064041E-29, 2.0921338E-26, 5.2846712E-24,
        8.5189081E-22, 9.2610320E-20, 7.0832621E-18, 3.9408644E-16,
        1.6382810E-14, 5.2018952E-13, 1.2847697E-11, 2.5062323E-10,
        3.9117887E-09, 4.9395557E-08, 5.0942322E-07, 4.3261255E-06,
        3.0466224E-05, 1.7901699E-04, 8.8231820E-04, 3.6643531E-03,
        1.2874141E-02, 3.8392534E-02, 9.7457619E-02, 2.1108208E-01,
        3.9082976E-01, 6.1955334E-01, 8.4179756E-01, 9.8105426E-01,
    };
    
    double weights [32] = {
        4.1125314E-01, 3.1721873E-01, 2.7607288E-01, 2.5139122E-01,
        2.3442298E-01, 2.2182848E-01, 2.1201056E-01, 2.0409313E-01,
        1.9754880E-01, 1.9203791E-01, 1.8733019E-01, 1.8326318E-01,
        1.7971845E-01, 1.7660734E-01, 1.7386204E-01, 1.7142969E-01,
        1.6926843E-01, 1.6734470E-01, 1.6563126E-01, 1.6410584E-01,
        1.6275006E-01, 1.6154868E-01, 1.6048903E-01, 1.5956052E-01,
        1.5875434E-01, 1.5806318E-01, 1.5748099E-01, 1.5700284E-01,
        1.5662481E-01, 1.5634386E-01, 1.5615777E-01, 1.5606509E-01,
    };
    
    double ret = 0.0;
    double a;
    for (size_t i = 0; i < 32; ++i) {
        a = x*expxi2[i];
        ret += a / (1.0 + a) * weights[i];
    }
    ret *= 2.0; // Symetrical function
    return ret;
}

double SnbDiatomicSystem::wlod(double ku, double bl)
{
    return (bl > 0 ? 2.0*bl*(sqrt(1.0+ku/bl)-1.0) : 0.0);
}

double SnbDiatomicSystem::wdod(double ku, double bd)
{
    if (bd == 0.0)
        return 0.0;

    switch (m_dpdf) {
        case EXPONENTIAL: return (bd*equad(ku/bd));
        case TAILED_INVERSE_EXP: return (bd*hquad(ku/bd));
        default: return 0.0;
    }
}

double SnbDiatomicSystem::wvod(double ku, double wd, double wl)
{
    double wdk, wlk, omega;

    // When ku is zero (avoid dividing by zero)
    if (ku > 0) {
        // Avoid dividing by zero when ku == wd or ku == wl
        if (ku == wd) return wd;
        if (ku == wl) return wl;

        wdk = wd/ku; wdk = 1.0-wdk*wdk;
        wlk = wl/ku; wlk = 1.0-wlk*wlk;
        omega = 1.0/(wdk*wdk)+1.0/(wlk*wlk)-1.0;
        return ku*sqrt(1.0-1.0/sqrt(omega));
    } else
        return 0.0;
}
//
//double SnbDiatomicSystem::wvod(double ku, double wd, double wl)
//{
//    double wdk, wlk, omega;
//    wdk = wd/ku; wdk = 1.0-wdk*wdk;
//    wlk = wl/ku; wlk = 1.0-wlk*wlk;
//    omega = 1.0/(wdk*wdk)+1.0/(wlk*wlk)-1.0;
//    return ku*sqrt(1.0-1.0/sqrt(omega));
//}

double dhquad(double x, double r)
{
    double expxi2 [32] = {
        7.5943201E-49, 2.9875169E-43, 6.9915249E-39, 3.4952039E-35,
        6.1363102E-32, 4.9064041E-29, 2.0921338E-26, 5.2846712E-24,
        8.5189081E-22, 9.2610320E-20, 7.0832621E-18, 3.9408644E-16,
        1.6382810E-14, 5.2018952E-13, 1.2847697E-11, 2.5062323E-10,
        3.9117887E-09, 4.9395557E-08, 5.0942322E-07, 4.3261255E-06,
        3.0466224E-05, 1.7901699E-04, 8.8231820E-04, 3.6643531E-03,
        1.2874141E-02, 3.8392534E-02, 9.7457619E-02, 2.1108208E-01,
        3.9082976E-01, 6.1955334E-01, 8.4179756E-01, 9.8105426E-01,
    };

    double weights [32] = {
        4.1125314E-01, 3.1721873E-01, 2.7607288E-01, 2.5139122E-01,
        2.3442298E-01, 2.2182848E-01, 2.1201056E-01, 2.0409313E-01,
        1.9754880E-01, 1.9203791E-01, 1.8733019E-01, 1.8326318E-01,
        1.7971845E-01, 1.7660734E-01, 1.7386204E-01, 1.7142969E-01,
        1.6926843E-01, 1.6734470E-01, 1.6563126E-01, 1.6410584E-01,
        1.6275006E-01, 1.6154868E-01, 1.6048903E-01, 1.5956052E-01,
        1.5875434E-01, 1.5806318E-01, 1.5748099E-01, 1.5700284E-01,
        1.5662481E-01, 1.5634386E-01, 1.5615777E-01, 1.5606509E-01,
    };

    double ret = 0.0, tmp, r2 = r*r;
    for (size_t i = 0; i < 32; ++i) {
        tmp = 1.0+x*pow(expxi2[i],r2);
        ret += expxi2[i] / tmp * weights[i];
    }
    ret *= 2.0; // Symetrical function
    return ret;
}

double dequad(double x, double r)
{
    double expxi2 [32] = {
        7.5943201E-49, 2.9875169E-43, 6.9915249E-39, 3.4952039E-35,
        6.1363102E-32, 4.9064041E-29, 2.0921338E-26, 5.2846712E-24,
        8.5189081E-22, 9.2610320E-20, 7.0832621E-18, 3.9408644E-16,
        1.6382810E-14, 5.2018952E-13, 1.2847697E-11, 2.5062323E-10,
        3.9117887E-09, 4.9395557E-08, 5.0942322E-07, 4.3261255E-06,
        3.0466224E-05, 1.7901699E-04, 8.8231820E-04, 3.6643531E-03,
        1.2874141E-02, 3.8392534E-02, 9.7457619E-02, 2.1108208E-01,
        3.9082976E-01, 6.1955334E-01, 8.4179756E-01, 9.8105426E-01,
    };

    double weights [32] = {
        4.1125314E-01, 3.1721873E-01, 2.7607288E-01, 2.5139122E-01,
        2.3442298E-01, 2.2182848E-01, 2.1201056E-01, 2.0409313E-01,
        1.9754880E-01, 1.9203791E-01, 1.8733019E-01, 1.8326318E-01,
        1.7971845E-01, 1.7660734E-01, 1.7386204E-01, 1.7142969E-01,
        1.6926843E-01, 1.6734470E-01, 1.6563126E-01, 1.6410584E-01,
        1.6275006E-01, 1.6154868E-01, 1.6048903E-01, 1.5956052E-01,
        1.5875434E-01, 1.5806318E-01, 1.5748099E-01, 1.5700284E-01,
        1.5662481E-01, 1.5634386E-01, 1.5615777E-01, 1.5606509E-01,
    };

    double ret = 0.0, tmp, r2 = r*r;
    for (size_t i = 0; i < 32; ++i) {
        tmp = 1.0+x*pow(expxi2[i],r2); tmp *= tmp;
        ret += expxi2[i] / tmp * weights[i];
    }
    ret *= 2.0; // Symetrical function
    return ret;
}

double SnbDiatomicSystem::dwlod(double x, double r)
{
    return  (2*x*r+(1-r*r)*sqrt(1+2*x))/((1-r*r+2*x)*sqrt(1+2*x));
    //return 0.0;
}

double SnbDiatomicSystem::dwdod(double x, double r)
{
    switch (m_dpdf) {
        case EXPONENTIAL: return dequad(x, r);
        case TAILED_INVERSE_EXP: return dhquad(x, r);
        default: return 0.0;
    }

}

void SnbDiatomicSystem::wquad(const double ku, const double bd, const double bl, double* const p_kbw)
{

   double points [4] = {
      6.94318442029737E-002, 0.330009478207572, 0.669990521792428, 0.930568155797026
   };

   double weights [4] = {
       0.173927422568727, 0.326072577431273, 0.326072577431273, 0.173927422568727
   };
/*
   double points [8] = {
       1.9855071751231856E-002, 0.10166676129318664, 0.23723379504183550, 0.40828267875217511,
       0.59171732124782495    , 0.76276620495816450, 0.89833323870681336, 0.98014492824876820
   };

   double weights [8] = {
       5.0614268145188129E-002, 0.11119051722668724, 0.15685332293894363, 0.18134189168918100    ,
       0.18134189168918100    , 0.15685332293894363, 0.11119051722668724, 5.0614268145188129E-002
   };

   double points [64] = {
       3.4747913211369275E-004, 1.8299416140222236E-003, 4.4933142616276900E-003, 8.3318730576869005E-003,
       1.3336586105044457E-002, 1.9495600173973004E-002, 2.6794312570798562E-002, 3.5215413934030104E-002,
       4.4738931460748532E-002, 5.5342277002442930E-002, 6.7000300922953560E-002, 7.9685351873709898E-002,
       9.3367342438601231E-002, 0.10801382052832936    , 0.12359004636973403    , 0.14005907491419456    ,
       0.15738184347288336    , 0.17551726437267123    , 0.19442232241380336    , 0.21405217689868300    ,
       0.23436026799005272    , 0.25529842714647361    , 0.27681699137326787    , 0.29886492101800410    ,
       0.32138992083116591    , 0.34433856400489471    , 0.36765641889561618    , 0.39128817812999650    ,
       0.41517778978800363    , 0.43926859035193977    , 0.46350343910610048    , 0.48782485366828776    ,
       0.51217514633171224    , 0.53649656089389963    , 0.56073140964806034    , 0.58482221021199643    ,
       0.60871182187000361    , 0.63234358110438360    , 0.65566143599510573    , 0.67861007916883398    ,
       0.70113507898199567    , 0.72318300862673213    , 0.74470157285352645    , 0.76563973200994728    ,
       0.78594782310131706    , 0.80557767758619681    , 0.82448273562732854    , 0.84261815652711669    ,
       0.85994092508580544    , 0.87640995363026586    , 0.89198617947167080    , 0.90663265756139866    ,
       0.92031464812629005    , 0.93299969907704638    , 0.94465772299755690    , 0.95526106853925141    ,
       0.96478458606596984    , 0.97320568742920166    , 0.98050439982602733    , 0.98666341389495571    ,
       0.99166812694231299    , 0.99550668573837220    , 0.99817005838597772    , 0.99965252086788614
    };

    double weights [64] = {
       8.9164036084755525E-004, 2.0735166302810961E-003, 3.2522289844890591E-003, 4.4233799131819093E-003,
       5.5840697300655181E-003, 6.7315239483592805E-003, 7.8630152380123313E-003, 8.9758578878486421E-003,
       1.0067411576765098E-002, 1.1135086904191616E-002, 1.2176351284355401E-002, 1.3188734857527340E-002,
       1.4169836307129757E-002, 1.5117328536201268E-002, 1.6028964177425789E-002, 1.6902580918570807E-002,
       1.7736106628441217E-002, 1.8527564270120020E-002, 1.9275076589307834E-002, 1.9976870566360144E-002,
       2.0631281621311788E-002, 2.1236757561826813E-002, 2.1791862264661715E-002, 2.2295279081878280E-002,
       2.2745813963709102E-002, 2.3142398290657205E-002, 2.3484091408105021E-002, 2.3770082857415123E-002,
       2.3999694298229155E-002, 2.4172381117401474E-002, 2.4287733720751756E-002, 2.4345478504569855E-002,
       2.4345478504569844E-002, 2.4287733720751728E-002, 2.4172381117401467E-002, 2.3999694298229141E-002,
       2.3770082857415137E-002, 2.3484091408105007E-002, 2.3142398290657194E-002, 2.2745813963709085E-002,
       2.2295279081878269E-002, 2.1791862264661704E-002, 2.1236757561826802E-002, 2.0631281621311774E-002,
       1.9976870566360133E-002, 1.9275076589307820E-002, 1.8527564270120037E-002, 1.7736106628441228E-002,
       1.6902580918570807E-002, 1.6028964177425813E-002, 1.5117328536201213E-002, 1.4169836307129773E-002,
       1.3188734857527340E-002, 1.2176351284355411E-002, 1.1135086904191630E-002, 1.0067411576765103E-002,
       8.9758578878486629E-003, 7.8630152380122723E-003, 6.7315239483591704E-003, 5.5840697300654375E-003,
       4.4233799131819683E-003, 3.2522289844891311E-003, 2.0735166302811459E-003, 8.9164036084813313E-004
    };
*/
    double dwd = 0.0, dwl = 0.0;
    double mku, mkubd, mkubl, xd, rd, xl, rl;
    for (int i = 0; i<4; i++) {
       mku = p_kbw[0] + ku * points[i];
       mkubd = p_kbw[1] + bd * ku * points[i];
       mkubl = p_kbw[2] + bl * ku * points[i];

       xd = mku * mku / mkubd;
       rd = bd * mku / mkubd;  
       dwd += weights[i] * dwdod(xd, rd);

       xl = mku * mku / mkubl / 2.0;
       rl = bl * mku / mkubl;
       dwl += weights[i] * dwlod(xl, rl);
    }

    p_kbw[0] += ku;
    p_kbw[1] += ku*bd;
    p_kbw[2] += ku*bl;
    p_kbw[3] += ku*dwd;
    p_kbw[4] += ku*dwl;
}

