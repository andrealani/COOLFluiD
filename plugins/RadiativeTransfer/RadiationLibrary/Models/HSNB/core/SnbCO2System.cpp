#include "SnbCO2System.h"

#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cstdlib>
#include <cmath>

#include "StringUtils.h"

using namespace std;
using namespace COOLFluiD;
using namespace COOLFluiD::RadiativeTransfer;

SnbCO2System::SnbCO2System(const std::string& directory, const std::string& nup, const ThermoData& thermo)

{

    m_species_index = thermo.speciesIndex("CO2");
    
    // Add the full path to the default data directory

    m_directory = directory+"/data/modele_msbe/CO2";

//    cout << "DIRECTORY SnbCO2System: " << m_directory << " \n";

    if (m_species_index < 0) {
        cout << "Error loading CO2 system ! "
             << "The species apparently isn't present in the atm. composition." << endl;
        cf_assert(m_species_index>=0);
    }

    m_nbLocalParams=nbLocalParameters;

    // Determine the min and max band
    determineBandRange();

    if (thermo.sigMin()!=-1) {
        CFuint lowWav=round(thermo.sigMin() / 25.0);
        CFuint hiWav=round(thermo.sigMax() / 25.0);

        if (lowWav<hiWav) {
            if (waveNumberIsInBandRange(lowWav)) {
                m_band1_down=lowWav;
            }
            if (waveNumberIsInBandRange(hiWav)) {
                m_bandn_down=hiWav;
            }

            // If the custom spectrum is not contained by the database information skip species
            if (m_band1_down>hiWav || m_bandn_down<lowWav) {
                m_band1_down=m_bandn_down+1;
            }
        }

        m_nbands_down = m_bandn_down - m_band1_down + 1;

    }
    
    // Load the temperature grid and determine number of parameters
    loadTemperatureGrid();
    m_nParamdata = m_nparams*m_nbands_down;

    if (nup == "LS")
        m_nup = LINDQUIST_SIMMONS;
    else
        m_nup = CURTIS_GODSON;

    m_lorentz = false;

    CFLog(INFO, "SnbCO2System::SnbCO2System => Loading CO2 system "
         << " (" << m_band1_down << " - " << m_bandn_down << ") "
         << ", Non uniform path approximation = " 
         << (m_nup == CURTIS_GODSON ? "CG" : "LS" ) << "\n");
    CFLog(INFO, "SnbCO2System::SnbCO2System => "<< m_nbands_down << "\n");
    
    // Allocate storage for the parameter information
    mp_params = new double [m_npoints * m_nbands_down * m_nparams];
    mp_locparams_up = NULL;
    mp_locparams_down = NULL;
    
    // Now load the parameter information

    for (size_t i = 0; i < m_nbands_down; ++i) {
        loadBandParameters(i);
    }


}

bool SnbCO2System::waveNumberIsInBandRange(CFuint waveNumber)
{
    return ((waveNumber>=lowBand()) && (waveNumber<=highBand()));
}

SnbCO2System::SnbCO2System(const SnbCO2System& system)
    : m_directory(system.m_directory),
      m_species_index(system.m_species_index),
      m_nup(system.m_nup),
      m_lorentz(system.m_lorentz),
      m_band1_up(system.m_band1_up),
      m_bandn_up(system.m_bandn_up),
      m_nbands_up(system.m_nbands_up),
      m_band1_down(system.m_band1_down),
      m_bandn_down(system.m_bandn_down),
      m_nbands_down(system.m_nbands_down),
      m_nparams(system.m_nparams),
      m_npoints(system.m_npoints),
      mp_t(system.mp_t == NULL ? NULL : new float [m_npoints]),
      mp_params(system.mp_params == NULL ? NULL :
          new double [m_npoints*m_nbands_down*m_nparams]),
      mp_locparams_up(NULL),
      mp_locparams_down(NULL)
{
    copy(system.mp_t, system.mp_t+m_npoints, mp_t);
    copy(
        system.mp_params, system.mp_params+m_npoints*m_nbands_down*m_nparams,
        mp_params);
}

void swap(SnbCO2System& s1, SnbCO2System& s2)
{
    std::swap(s1.m_directory, s2.m_directory);
    std::swap(s1.m_species_index, s2.m_species_index);
    std::swap(s1.m_nup, s2.m_nup);
    std::swap(s1.m_lorentz, s2.m_lorentz);
    std::swap(s1.m_band1_up, s2.m_band1_up);
    std::swap(s1.m_bandn_up, s2.m_bandn_up);
    std::swap(s1.m_nbands_up, s2.m_nbands_up);
    std::swap(s1.m_band1_down, s2.m_band1_down);
    std::swap(s1.m_bandn_down, s2.m_bandn_down);
    std::swap(s1.m_nbands_down, s2.m_nbands_down);
    std::swap(s1.m_nparams, s2.m_nparams);
    std::swap(s1.m_npoints, s2.m_npoints);
    std::swap(s1.mp_t, s2.mp_t);
    std::swap(s1.mp_params, s2.mp_params);
    std::swap(s1.mp_locparams_up, s2.mp_locparams_up);
    std::swap(s1.mp_locparams_down, s2.mp_locparams_down);
}

SnbCO2System::~SnbCO2System()
{
    delete [] mp_t;
    delete [] mp_params;
    if(mp_locparams_up)
        delete [] mp_locparams_up;
    if(mp_locparams_down)
        delete [] mp_locparams_down;
}

void SnbCO2System::addStateParams(ThermoData &thermo, HSNBCO2ParameterSet &co2Params, CFreal cellDistance, CFuint localCellID, CFreal sig)
{
//            std::cout << "NBCELLS" << co2Params.nbCells<< std::endl;
            int b = int(sig / 25.0);


            // Convert band to local indexing, check whether Band is in database
            if (b < lowBand() || b > highBand()) {
        //        CFLog(VERBOSE, "SnbDiatomicSystem::addStateParams =>NONTHICK, " << this->speciesName()<<", b=" << b << " out of bounds. \n");
                co2Params.addState(0.0);
                return;
            }

            b -= lowBand();


            if (thermo.usePrecomputedDiatomicParameters()) {


                switch (m_nup) {
                case CURTIS_GODSON:
                {

                    m_tempKD = mp_locparams_down[localCellID*m_ndata+b*m_nbLocalParams+0] * cellDistance;
                    m_tempKL = mp_locparams_down[localCellID*m_ndata+b*m_nbLocalParams+1] * cellDistance;

                    // Doppler broadening, formal Curtis-Godson approximation
                    m_tempBetaD = m_tempKD/(mp_locparams_down[localCellID*m_ndata+b*m_nbLocalParams+2]+1e-20);
                    // Lorentz broadening, classical Curtis-Godson approximation
                    m_tempBetaL = m_tempKL*mp_locparams_down[localCellID*m_ndata+b*m_nbLocalParams+3];


                    if (m_tempKD == 0.0 || m_tempBetaD == 0.0)
                        m_tempBetaD = 1.0e-20;
                    else
                        m_tempBetaD = m_tempKD / m_tempBetaD;

                    if (m_tempKL == 0.0 || m_tempBetaL == 0.0)
                        m_tempBetaL = 1.0e-20;
                    else
                        m_tempBetaL = m_tempBetaL / m_tempKL;


                    //New optical thickness
                    if(m_tempBetaD == 0.0)
                        //wlod(kl,bl);
                        m_tempKappa= wlod(m_tempKD,m_tempBetaL);
                    else if (m_tempKL == 0.0)
                        //wdod(kd,bd)
                        m_tempKappa = wdod(m_tempKD,m_tempBetaD);
                    else
                        //wvod(kd, kl, wdod(kd, bd), wlod(kl, bl))
                        m_tempKappa = wvod(m_tempKD, m_tempKL, wdod(m_tempKD, m_tempBetaD), wlod(m_tempKL, m_tempBetaL)); // Voigt profile


                    break;
                }
                case LINDQUIST_SIMMONS:
                {

                    m_tempKD = mp_locparams_down[localCellID*m_ndata+b*m_nbLocalParams+0] * cellDistance;
                    m_tempKL = mp_locparams_down[localCellID*m_ndata+b*m_nbLocalParams+1] * cellDistance;
                    m_tempBetaD = mp_locparams_down[localCellID*m_ndata+b*m_nbLocalParams+2];
                    m_tempBetaL = mp_locparams_down[localCellID*m_ndata+b*m_nbLocalParams+3];
                    if (m_tempBetaD == 0.0) m_tempBetaD = 1e-20;
                    if (m_tempBetaL == 0.0) m_tempBetaL = 1e-20;


                    double p_kb [6] = {0, 0, 0, 0, 0, 0};

                    if (co2Params.isEmpty()) {


                        p_kb[0] = m_tempKD;
                        p_kb[1] = m_tempKL;
                        p_kb[2] = m_tempKD/m_tempBetaD;
                        p_kb[3] = m_tempBetaL*m_tempKL;
                        p_kb[4] = wdod(m_tempKD, m_tempBetaD);
                        p_kb[5] = wlod(m_tempKL, m_tempBetaL); // Voigt profile

                        if(m_tempKD == 0.0)
                            m_tempKappa = wlod(m_tempKL,m_tempBetaL);
                        else if (m_tempKL == 0.0)
                            m_tempKappa = wdod(m_tempKD,m_tempBetaD);
                        else
                            m_tempKappa = wvod(p_kb[0], p_kb[1], p_kb[4], p_kb[5]);

                    } else {

                        wquad(m_tempKD, m_tempKL, m_tempBetaD, m_tempBetaL, p_kb);

                        if (p_kb[0] == 0.0)
                            m_tempKappa = p_kb[5];
                        else if (p_kb[1] == 0.0)
                            m_tempKappa = p_kb[4];
                        else
                            m_tempKappa = wvod(p_kb[0], p_kb[1], p_kb[4], p_kb[5]);
                    }

                    break;

                }
                default:
                    std::cout << "This integration model is not supported yet!" << std::endl;

                }

//                std::cout << "SnbCO2System::addStateParams => Precomputed mp_locparams_down[" <<localCellID<<"*"<<m_ndata<<"+"<<b<<"*"<<m_nbLocalParams<<"+0]=" << mp_locparams_down[localCellID*m_ndata+b*m_nbLocalParams+0] << std::endl;
//                std::cout << "SnbCO2System::addStateParams => Precomputed mp_locparams_down[" <<localCellID<<"*"<<m_ndata<<"+"<<b<<"*"<<m_nbLocalParams<<"+1]=" << mp_locparams_down[localCellID*m_ndata+b*m_nbLocalParams+1] << std::endl;
//                std::cout << "SnbCO2System::addStateParams => Precomputed mp_locparams_down[" <<localCellID<<"*"<<m_ndata<<"+"<<b<<"*"<<m_nbLocalParams<<"+2]=" << mp_locparams_down[localCellID*m_ndata+b*m_nbLocalParams+2] << std::endl;
//                std::cout << "SnbCO2System::addStateParams => Precomputed mp_locparams_down[" <<localCellID<<"*"<<m_ndata<<"+"<<b<<"*"<<m_nbLocalParams<<"+3]=" << mp_locparams_down[localCellID*m_ndata+b*m_nbLocalParams+3] << std::endl;
            }
            else {
                //TODO: USE PROPER FUNCTION HEADER
                double kd, kl, betad, betal;

                getAbsorptionSingleBand(thermo,localCellID,b,kd,kl,betad,betal);

//                std::cout << "SnbCO2System::addStateParams => Single Band kappaD="<< kd << std::endl;
//                std::cout << "SnbCO2System::addStateParams => Single Band kappaL=" << kl << std::endl;
//                std::cout << "SnbCO2System::addStateParams => Single Band betaD=" << betad << std::endl;
//                std::cout << "SnbCO2System::addStateParams => Single Band betaL=" << betal << std::endl;


                switch (m_nup) {
                case CURTIS_GODSON:
                {

                    m_tempKD = kd * cellDistance;
                    m_tempKL = kl * cellDistance;

                    // Doppler broadening, formal Curtis-Godson approximation
                    m_tempBetaD = m_tempKD/(betad+1e-20);
                    // Lorentz broadening, classical Curtis-Godson approximation
                    m_tempBetaL = m_tempKL*betal;


                    if (m_tempKD == 0.0 || m_tempBetaD == 0.0)
                        m_tempBetaD = 1.0e-20;
                    else
                        m_tempBetaD = m_tempKD / m_tempBetaD;

                    if (m_tempKL == 0.0 || m_tempBetaL == 0.0)
                        m_tempBetaL = 1.0e-20;
                    else
                        m_tempBetaL = m_tempBetaL / m_tempKL;


                    //New optical thickness
                    if(m_tempBetaD == 0.0)
                        //wlod(kl,bl);
                        m_tempKappa= wlod(m_tempKD,m_tempBetaL);
                    else if (m_tempKL == 0.0)
                        //wdod(kd,bd)
                        m_tempKappa = wdod(m_tempKD,m_tempBetaD);
                    else
                        //wvod(kd, kl, wdod(kd, bd), wlod(kl, bl))
                        m_tempKappa = wvod(m_tempKD, m_tempKL, wdod(m_tempKD, m_tempBetaD), wlod(m_tempKL, m_tempBetaL)); // Voigt profile


                    break;
                }
                case LINDQUIST_SIMMONS:
                {

                    m_tempKD = kd * cellDistance;
                    m_tempKL = kl * cellDistance;
                    m_tempBetaD = betad;
                    m_tempBetaL = betal;
                    if (m_tempBetaD == 0.0) m_tempBetaD = 1e-20;
                    if (m_tempBetaL == 0.0) m_tempBetaL = 1e-20;

                    double p_kb [6] = {0, 0, 0, 0, 0, 0};

                    if (co2Params.isEmpty()) {


                        p_kb[0] = m_tempKD;
                        p_kb[1] = m_tempKL;
                        p_kb[2] = m_tempKD/m_tempBetaD;
                        p_kb[3] = m_tempBetaL*m_tempKL;
                        p_kb[4] = wdod(m_tempKD, m_tempBetaD);
                        p_kb[5] = wlod(m_tempKL, m_tempBetaL); // Voigt profile

                        if(m_tempKD == 0.0)
                            m_tempKappa = wlod(m_tempKL,m_tempBetaL);
                        else if (m_tempKL == 0.0)
                            m_tempKappa = wdod(m_tempKD,m_tempBetaD);
                        else
                            m_tempKappa = wvod(p_kb[0], p_kb[1], p_kb[4], p_kb[5]);

                    } else {

                        wquad(m_tempKD, m_tempKL, m_tempBetaD, m_tempBetaL, p_kb);

                        if (p_kb[0] == 0.0)
                            m_tempKappa = p_kb[5];
                        else if (p_kb[1] == 0.0)
                            m_tempKappa = p_kb[4];
                        else
                            m_tempKappa = wvod(p_kb[0], p_kb[1], p_kb[4], p_kb[5]);
                    }


                    break;

                }
                default:
                    std::cout << "This integration model is not supported yet!" << std::endl;

                }
            }

//            std::cout << "co2Params.kappa=" << co2Params.kappa << "+" << "m_tempKappa=" << m_tempKappa << std::endl;
            m_tempKappa+=co2Params.kappa;
//            std::cout << "m_tempKappa=" << m_tempKappa << " co2Params.nbCells=" << co2Params.nbCells << std::endl;

            co2Params.addState(m_tempKappa);

//            co2Params.print();
}

void SnbCO2System::getAbsorptionSingleBand(ThermoData &thermo, CFuint cellID, double band, double &KD, double &KL, double &betaD, double &betaL)
{
    double ptot, T, pco2, xco2;

        thermo.setState(cellID);

        ptot = thermo.P();
        T    = thermo.Tr();
        pco2 = thermo.N(m_species_index) * KB * T;
        xco2 = thermo.X(m_species_index);


        float* p_lower = std::lower_bound(mp_t, mp_t+m_npoints, (float)T);

        size_t it;
        if (p_lower == mp_t+m_npoints)
            it = m_npoints-2;
        else
            it = (mp_t == p_lower ? 0 : std::distance(mp_t, p_lower)-1);

        const double t1 = mp_t[it];
        const double t2 = mp_t[it+1];

        const double w1 = (t2-T)/(t2-t1);
        const double w2 = (T-t1)/(t2-t1);

        const size_t nParamdata = m_nparams*m_nbands_down;
        const double* const p1 = mp_params + it*nParamdata;
        const double* const p2 = p1 + nParamdata;


        int i = m_nparams*band;

//        std::cout << "cellID= " << cellID << ", it=" << it << std::endl;
//        std::cout << "i=" << i << " m_nparams=" << m_nparams << ", band=" << band << std::endl;

        //Interpolate values
        KD = (w1*p1[i] + w2*p2[i]) * pco2; // kappa doppler
        KL = (w1*p1[i+3] + w2*p2[i+3]) * pco2; // kappa lorentz
        betaD = w1*p1[i+2] + w2*p2[i+2] ; // beta doppler
        betaL = ( (w1*p1[i+5] + w2*p2[i+5]) * xco2
                + (w1*p1[i+6] + w2*p2[i+6]) ) * ptot ; // beta lorentz


        if (KD<0.0) {
            KD=0.0;
        }
        if (KL<0.0) {
            KL=0.0;
        }

        if (betaD<0.0) {
            betaD=0.0;
        }
        if (betaL<0.0) {
            betaL=0.0;
        }





}

double SnbCO2System::getEmissionSingleBand(ThermoData &thermo, CFuint cellID, double band)
{

    double eta;
    thermo.setState(cellID);
    double T=thermo.Tr();

    float* p_lower = std::lower_bound(mp_t, mp_t+m_npoints, (float)T);
    size_t it;
    if (p_lower == mp_t+m_npoints)
        it = m_npoints-2;
    else
        it = (mp_t == p_lower ? 0 : std::distance(mp_t, p_lower)-1);

    // Step 2: Compute the linear interpolation
    const double t1 = mp_t[it];
    const double t2 = mp_t[it+1];

    const double w1 = (t2-T)/(t2-t1);
    const double w2 = (T-t1)/(t2-t1);

    const size_t nParamdata = m_nparams*m_nbands_down;
    const double* const p1 = mp_params + it*nParamdata;
    const double* const p2 = p1 + nParamdata;

    int i = m_nparams*band;


    if (m_lorentz) {
        eta = w1*p1[i+4] + w2*p2[i+4]; // eta lorentz
    }
    else {
        eta = w1*p1[i+1] + w2*p2[i+1]; // eta doppler
    }

    if (eta<0.0) {
        return 0.0;
    }
    else {
        return eta;
    }
}

double SnbCO2System::getEmissionSingleBandFixedCell(CFuint band)
{

    int i = m_nparams*band;


    if (m_lorentz) {
        m_eta = m_w1*m_p1[i+4] + m_w2*m_p2[i+4]; // eta lorentz
    }
    else {
        m_eta = m_w1*m_p1[i+1] + m_w2*m_p2[i+1]; // eta doppler
    }

    if (m_eta<0.0) {
        return 0.0;
    }
    else {
        return m_eta;
    }
}

double SnbCO2System::setupParamInterpolation(ThermoData& thermo, CFuint cellID)
{

//    std::cout << "SnbCO2System::setupParamInterpolation => start" << std::endl;
    thermo.setState(cellID);
    double T=thermo.Tr();
    float* p_lower = std::lower_bound(mp_t, mp_t+m_npoints, (float)T);

    if (p_lower == mp_t+m_npoints)
        m_it = m_npoints-2;
    else
        m_it = (mp_t == p_lower ? 0 : std::distance(mp_t, p_lower)-1);

    // Step 2: Compute the linear interpolation
    m_t1 = mp_t[m_it];
    m_t2 = mp_t[m_it+1];

    m_w1 = (m_t2-T)/(m_t2-m_t1);
    m_w2 = (T-m_t1)/(m_t2-m_t1);

    m_p1 = mp_params + m_it*m_nParamdata;
    m_p2 = m_p1 + m_nParamdata;
//    std::cout << "SnbCO2System::setupParamInterpolation => end" << std::endl;
}

void SnbCO2System::getParameters(
    double T, double* const p_params) const
{
    // Step 1: Determine the indices for T which either bound the
    // interval in which T falls or the correct interval which will
    // be used to extrapolate if it falls out of the grid boundaries

    // Index for T:  0 <= it < npoints-2
    float* p_lower = std::lower_bound(mp_t, mp_t+m_npoints, (float)T);
    size_t it;
    if (p_lower == mp_t+m_npoints)
        it = m_npoints-2;
    else
        it = (mp_t == p_lower ? 0 : std::distance(mp_t, p_lower)-1);
    
    // Step 2: Compute the linear interpolation
    const double t1 = mp_t[it];
    const double t2 = mp_t[it+1];

    const double w1 = (t2-T)/(t2-t1);
    const double w2 = (T-t1)/(t2-t1);

    const size_t nParamdata = m_nparams*m_nbands_down;
    const double* const p1 = mp_params + it*nParamdata;
    const double* const p2 = p1 + nParamdata;

//    std::cout<< "it= " << it << std::endl;
//    std::cout << "nParamdata=" << nParamdata << " m_nparams=" << m_nparams << " m_nbands_down=" << m_nbands_down << std::endl;

    for (int i = 0; i < nParamdata; ++i) {
        p_params[i] = w1*p1[i] + w2*p2[i];

//        std::cout << "SnbCO2System::getParameters =>"<< " w1=" << w1 << ", p1[" << i << "]=" << p1[i] << ", w2=" << w2
//                  << ", p2[" << i << "]=" << p2[i]  << std::endl;

        if (p_params[i]<0.0)
            p_params[i]=0.0;
    }
}

double eqInt(double nu,double t)
// nu (cm-1), t (K), eqInt (W.cm-2.sr-1.cm) 
{
    double const c1=2*HP*pow(C0*100,2.0);
    double const c2= HP*C0*100/KB;
    double res;

    if (t == 0.) { return res=0.;}
    else { return res=c1*pow(nu,3.)/(exp(c2*nu/t)-1.);}
}

void SnbCO2System::setupLocalParameters(ThermoData& thermo)
{
    if (thermo.usePrecomputedDiatomicParameters()) {
        double ptot, t, tv, pco2, xco2, pmoy;
        double *p_params = new double [m_nparams*m_nbands_down];


        if (mp_locparams_down) {
            delete [] mp_locparams_down;
        }

        //    cout << "SnbCO2System::setupLocalParameters => thermo.nCells()*m_ndata=" << thermo.nCells()*m_ndata << endl;


        // Down local parameters (transmissivity calculations)
        mp_locparams_down = new double [thermo.nCells()*m_ndata];


        for (int i=0; i<thermo.nCells(); i++) {
            thermo.setState(i);

            ptot = thermo.P();
            t    = thermo.Tr();
            pco2 = thermo.N(m_species_index) * KB * t;



            xco2 = thermo.X(m_species_index);

//            cout << m_species_index <<" ptot=" << ptot << " t=" << t << " pco2=" << pco2 <<" xco2=" << xco2 << endl;

//            std::cout << "SnbCO2System::setupLocalParameters => State=" << i << std::endl;
            getParameters(t, p_params);

            for (int b = 0; b < m_nbands_down; ++b) {
                mp_locparams_down[i*m_ndata+b*m_nbLocalParams+0] = p_params[b*m_nparams+0] * pco2; // kappa doppler
                mp_locparams_down[i*m_ndata+b*m_nbLocalParams+1] = p_params[b*m_nparams+3] * pco2; // kappa lorentz
                mp_locparams_down[i*m_ndata+b*m_nbLocalParams+2] = p_params[b*m_nparams+2] ; // beta doppler
                mp_locparams_down[i*m_ndata+b*m_nbLocalParams+3] = ( p_params[b*m_nparams+5] * xco2
                        + p_params[b*m_nparams+6] ) * ptot ; // beta lorentz
                if (m_lorentz) {
                    mp_locparams_down[i*m_ndata+b*m_nbLocalParams+4] = p_params[b*m_nparams+4]; // eta lorentz
                }
                else {
                    mp_locparams_down[i*m_ndata+b*m_nbLocalParams+4] = p_params[b*m_nparams+1]; // eta doppler
                }

            }

            // Avg pressure is precomputed in setupProfileType

        }


        delete [] p_params;
    }
}

double SnbCO2System::getLocalParameter(const int& i, const int& j, const int& k) const
// Cell (i), Band (j) and Parameter (k) indices
{
    return mp_locparams_down[i*m_nbands_down+j];
}

void SnbCO2System::determineBandRange()
{

    // Down SNB for CO2 : bands of 25 cm-1
    m_band1_down = 10;
    m_bandn_down = 332;
    m_nbands_down = m_bandn_down - m_band1_down + 1;

    m_ndata=m_nbLocalParams*m_nbands_down;

}

void SnbCO2System::loadTemperatureGrid()
{
    // Open the file corresponding to the first band
    std::ifstream file((m_directory+"/"+bandFilename(m_band1_down)).c_str());
    
    // Read the first temperature from the file and count the number of
    // parameters (extra columns)
    std::string line;
    std::vector<float>  t(2);
    
    std::getline(file, line);
    std::stringstream ss(line);
    ss >> t[0];
    
    m_nparams = 0;
    while (ss >> line)
        m_nparams++;
    
    // Read the rest of the temperatures from the file
    while (std::getline(file, line)) {
        sscanf(line.c_str(), "%f", &t.back());
        t.push_back(0.0f);
    }
    
    t.pop_back();
    file.close();
    
    // Compute sizes
    m_npoints = t.size();
    
    // Allocate storage for temperatures
    mp_t = new float [m_npoints];

    // Copy the temperatures to their storage
    std::copy(t.begin(), t.end(), mp_t);
}


void SnbCO2System::loadBandParameters(const size_t& iband)
{
    // Open the file corresponding to the band
    std::ifstream file((m_directory+"/"+bandFilename(iband+m_band1_down)).c_str());
    
    // Read the parameter information
    std::string line;
    for (size_t ipoint = 0; ipoint < m_npoints; ++ipoint) {
        // Get pointer to storage location for this point/band
        double* const p_params =
            mp_params + m_nparams * (iband + m_nbands_down * ipoint);
        // Skip temperature
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


double hquad_alpha(double x)
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

    double alpha = 0.3;
    double ret = 0.0;
    for (size_t i = 0; i < 32; ++i)
        ret += (std::pow(1.0 + x*expxi2[i], alpha)-1.0) * weights[i];
    ret *= 2.0 / alpha ; // Symetrical function
    return ret;
}

double SnbCO2System::wlod(double kl, double bl)
{
    return (bl*(sqrt(1.0+2.0*kl/bl)-1.0));
}

double SnbCO2System::wdod(double kd, double bd)
{
    return (bd*hquad_alpha(kd/bd)); // H-alpha dstribution
}

double SnbCO2System::wvod(
         double kd, double kl, double wd, double wl)
{
    double wdk, wlk, omega;
    wdk = wd/kd; wdk = 1.0-wdk*wdk;
    wlk = wl/kl; wlk = 1.0-wlk*wlk;
    omega = 1.0/(wdk*wdk)+1.0/(wlk*wlk)-1.0;
    if(m_lorentz)
        return kl*sqrt(1.0-1.0/sqrt(omega));
    else
        return kd*sqrt(1.0-1.0/sqrt(omega));
}


double dhquad_alpha(double x, double r)
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

    double alpha = 0.3;
    double ret = 0.0, tmp, r2 = r*r;
    for (size_t i = 0; i < 32; ++i) {
        tmp = pow(1.0+x*pow(expxi2[i],r2) , 1.0-alpha) ;
        ret += expxi2[i] / tmp * weights[i];
    }
    ret *= 2.0; // Symetrical function
    return ret;
}

double SnbCO2System::dwlod(double x, double r)
{
    return  (2*x*r+(1-r*r)*sqrt(1+2*x))/((1-r*r+2*x)*sqrt(1+2*x));
}

double SnbCO2System::dwdod(double x, double r)
{
    return dhquad_alpha(x, r);
}

void SnbCO2System::wquad(const double kd, const double kl,
                                const double bd, const double bl, double* const p_kbw)
{
    double points [4] = {
       6.94318442029737E-002, 0.330009478207572, 0.669990521792428, 0.930568155797026
    };

    double weights [4] = {
        0.173927422568727, 0.326072577431273, 0.326072577431273, 0.173927422568727
    };
    double dwd = 0.0, dwl = 0.0;
    double mkd, mkl, mkdbd, mklbl, xd, rd, xl, rl;
    for (int i = 0; i<4; i++) {
       mkd = p_kbw[0] + kd * points[i];
       mkl = p_kbw[1] + kl * points[i];
       mkdbd = p_kbw[2] + kd / bd * points[i];
       mklbl = p_kbw[3] + kl * bl * points[i];

       xd = mkdbd;
       rd = bd * mkdbd / mkd;  
       dwd += weights[i] * dwdod(xd, rd);

       xl = mkl * mkl / mklbl ;
       rl = bl * mkl / mklbl ;
       dwl += weights[i] * dwlod(xl, rl);
    }
    if (kd == 0.0) dwd = 0.0;
    if (kl == 0.0) dwl = 0.0;

    p_kbw[0] += kd;
    p_kbw[1] += kl;
    p_kbw[2] += kd/bd;
    p_kbw[3] += kl*bl;
    p_kbw[4] += kd*dwd;
    p_kbw[5] += kl*dwl;
}


void SnbCO2System::setupProfileType(const std::vector<CFreal>& cellVolumes, ThermoData& thermo)
{

    //Compute average pressure over the total volume by weighting cell pressure with the cell's volume
    std::vector<CFreal>::const_iterator it = cellVolumes.begin();
    CFreal ptot;
    CFreal pavg=0;
    CFreal totalVolume=0;
    for (int i=0; i<thermo.nCells(); i++) {
        thermo.setState(i);
        ptot = thermo.P();
        pavg += ptot * (*it);
        totalVolume+=(*it);

//        std::cout <<"SnbCO2System::setupProfileType => Current Volume: " << (*it)  << ",  ptot= " << ptot << std::endl;
        it++;
    }



    //cm^-3 or m^-3 does not matter (divided out)
    pavg = pavg / (totalVolume);

    std::cout << "SnbCO2System::setupProfileType => Total volume: " << totalVolume << std::endl;
//    std::cout << "pavg: " << pavg << std::endl;

    if (pavg >= 1000.0)
    {
        m_lorentz = true;
    }
    else {
        m_lorentz = false;
    }
}

double SnbCO2System::emittedPower(int cellID, ThermoData &thermo)
{
    thermo.setState(cellID);
//    double ptot  = thermo.P();

    //p = N/V*KV*T "ideal gas"
    double pa  = thermo.N(m_species_index) * KB * thermo.Tr();

//    According to JB: Use total pressure instead?
//    double pa = thermo.P();


    double sum = 0;
    if (thermo.usePrecomputedDiatomicParameters()) {

        //use eta
        for (int b = 0; b < m_nbands_down; ++b) {
            sum += mp_locparams_down[cellID*m_ndata+b*m_nbLocalParams+4];
            //            CFLog(INFO, "SnbCO2System::emittedPower => band " << b <<" emission coeff=" <<  mp_locparams_down[cellID*m_ndata+b*m_nbLocalParams+4] << "\n");
        }

    }
    else {
//        for (int b = 0; b < m_nbands_down; ++b) {
//            sum += getEmissionSingleBand(thermo,cellID,b);
//        }

        setupParamInterpolation(thermo,cellID);
        for (int b = 0; b < m_nbands_down; ++b) {
            sum += getEmissionSingleBandFixedCell(b);
//            CFLog(INFO, "SnbCO2System::emittedPower => band " << b <<" emission coeff=" <<  getEmissionSingleBandFixedCell(b) << " INPLACE" << "\n");

        }

    }

//    CFLog(INFO, "SnbCO2System::emittedPower => " << " emission coeff=" <<  sum * pa  * 1000.0 << " pa=" << pa << " ptot  =" << thermo.P()<< "\n");


    return (sum * pa * 1000.0);
}

void SnbCO2System::bandEmission(const int cellID, ThermoData &thermo, double * const p_emis)
{
    thermo.setState(cellID);

//    double ptot  = thermo.P();
    double pa  = thermo.N(m_species_index) * KB * thermo.Tr();


    if (thermo.usePrecomputedDiatomicParameters()) {        
        // eta
        for (int b = 0; b < m_nbands_down; ++b) {
            p_emis[b] = mp_locparams_down[cellID*m_ndata+b*m_nbLocalParams+4] * pa * 1000.0;
        }
    }
    else {
//        for (int b = 0; b < m_nbands_down; ++b) {
//            p_emis[b] = getEmissionSingleBand(thermo,cellID,b) * pa * 1000.0;
//        }

        //More efficient to setup the interpolation just once!
        setupParamInterpolation(thermo, cellID);
        for (int b = 0; b < m_nbands_down; ++b) {
            p_emis[b] = getEmissionSingleBandFixedCell(b) * pa * 1000.0;
        }
    }


}

double SnbCO2System::opticalThickness(const HSNBCO2ParameterSet &pathParams)
{

   CFLog(DEBUG_MAX, "SnbCO2::opticalThickness => kappa=" <<  pathParams.kappa << " \n");

   return pathParams.kappa;

}

double SnbCO2System::tau(const HSNBCO2ParameterSet &pathParams)
{
   return std::exp(-opticalThickness(pathParams));
}

