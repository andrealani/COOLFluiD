#include "SnbContinuumSystem.h"
#include "PhotonPath.h"
#include "Environment/SingleBehaviorFactory.hh"

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

SnbContinuumSystem::SnbContinuumSystem
(SpeciesLoadData loadData, const ThermoData& thermo, ContinuumSort sort)
  : m_directory(loadData.baseDirectory), m_sort(sort)
{
    m_directory += "/data/modele_msbe/";
    m_datadir = loadData.baseDirectory + "/data";
    m_nbParams = 0;
    m_species = loadData.speciesName;
    m_system  = loadData.systemName;

    m_directory+=m_species+"/"+m_system;

    mp_chi = new Chineq(m_system, m_species,m_datadir);

    m_species_index = thermo.speciesIndex(m_species);

    m_usePrecomputedParams=true;
    
    if (m_species_index < 0) {
        cout << "Error loading SNB system '" << m_system << "'! "
             << "The species isn't loaded for this system." << endl;
        exit(1);
    }

    setupProducts();
    
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
       }

       m_nbands = m_bandn - m_band1 + 1;

//       cout << "system: " << this->systemName()<<", lowWav " << lowWav << ", hiWav " << hiWav <<", m_band1 " << m_band1
//            << ", m_bandn " << m_bandn << "\n";
    }
    
    // Load the temperature grid and determine number of parameters
    loadTemperatureGrid();

    m_nParamdata = m_nparams*m_nbands;
    
    CFLog(VERBOSE, "SnbContinuumSystem::SnbContinuumSystem => Loading continuum system " << m_system
          << " (" << m_band1 << " - " << m_bandn << ") "
          << m_nparams << ", Continuum sort = "
          << (m_sort == BOUNDFREE ? "BF" : "FF") <<"\n");

    
    // Allocate storage for the parameter information
    mp_params = new double [m_npoints * m_nbands * m_nparams];
    mp_locparams = NULL;
    
    // Now load the parameter information
    #pragma omp parallel for
    for (size_t i = 0; i < m_nbands; ++i)
        loadBandParameters(i);
}

SnbContinuumSystem::SnbContinuumSystem(const SnbContinuumSystem& system)
    : m_directory(system.m_directory),
      m_species(system.m_species),
      m_products(system.m_products),
      m_system(system.m_system),
      m_species_index(system.m_species_index),
      m_sort(system.m_sort),
      m_band1(system.m_band1),
      m_bandn(system.m_bandn),
      m_nbands(system.m_nbands),
      m_nparams(system.m_nparams),
      m_npoints(system.m_npoints),
      m_datadir(system.m_datadir),
      mp_t(system.mp_t == NULL ? NULL : new float [m_npoints]),
      mp_params(system.mp_params == NULL ? NULL :
          new double [m_npoints*m_nbands*m_nparams]),
      mp_locparams(NULL)
{
    mp_chi = new Chineq(m_system, m_species,m_datadir);
    copy(system.mp_t, system.mp_t+m_npoints, mp_t);
    copy(
        system.mp_params, system.mp_params+m_npoints*m_nbands*m_nparams,
        mp_params);
}

void swap(SnbContinuumSystem& s1, SnbContinuumSystem& s2)
{
    std::swap(s1.m_directory, s2.m_directory);
    std::swap(s1.m_datadir, s2.m_datadir);
    std::swap(s1.m_species, s2.m_species);
    std::swap(s1.m_products, s2.m_products);
    std::swap(s1.m_system, s2.m_system);
    std::swap(s1.m_species_index, s2.m_species_index);
    std::swap(s1.m_sort, s2.m_sort);
    std::swap(s1.m_band1, s2.m_band1);
    std::swap(s1.m_bandn, s2.m_bandn);
    std::swap(s1.m_nbands, s2.m_nbands);
    std::swap(s1.m_nparams, s2.m_nparams);
    std::swap(s1.m_npoints, s2.m_npoints);
    std::swap(s1.mp_t, s2.mp_t);
    std::swap(s1.mp_params, s2.mp_params);
    std::swap(s1.mp_locparams, s2.mp_locparams);
}

SnbContinuumSystem::~SnbContinuumSystem()
{
    delete mp_chi;
    delete [] mp_t;
    delete [] mp_params;
    if (mp_locparams)
        delete [] mp_locparams;
}

void SnbContinuumSystem::getParameters(
    double T, double* const p_params)
{
    // Step 1: Determine the indices for T which either bound the
    // interval in which T falls or the correct interval which will
    // be used to extrapolate if it falls out of the grid boundaries
    
    // Index for T:  0 <= it < npoints-2
    m_pLower = std::lower_bound(mp_t, mp_t+m_npoints, (float)T);

    if (m_pLower == mp_t+m_npoints)
        m_it = m_npoints-2;
    else
        m_it = (mp_t == m_pLower ? 0 : std::distance(mp_t, m_pLower)-1);

    // Step 2: Compute the linear interpolation
    m_t1 = mp_t[m_it];
    m_t2 = mp_t[m_it+1];
    
    m_w1 = (m_t2-T)/(m_t2-m_t1);
    m_w2 = (T-m_t1)/(m_t2-m_t1);
    
    m_ndata = m_nparams*m_nbands;
    m_p1 = mp_params + m_it*m_ndata;
    m_p2 = m_p1 + m_ndata;
   
    // Log-linear interpolation if non zero values
    for (int i = 0; i < m_ndata; ++i) {
        if (m_p1[i] != 0.0 & m_p2[i] != 0.0) {
            p_params[i] = exp(m_w1*log(m_p1[i]) + m_w2*log(m_p2[i]));
        } else {
            p_params[i] = m_w1*m_p1[i] + m_w2*m_p2[i];
        }
    }

    // Enforces nul emission coefficients at low temperature
    // Extrapolation is not sufficiently accurate when the chi factor is taken into account
    if (T<m_t1 & m_nparams>2) {
      for (int i = 0; i < m_ndata; i=i+3) {
          p_params[i+1] = 0.0;
          p_params[i+2] = 0.0;
      } 
    }

}

void SnbContinuumSystem::getParametersSingleBand(double T, CFuint band,double* p_params)
{
    // Step 1: Determine the indices for T which either bound the
    // interval in which T falls or the correct interval which will
    // be used to extrapolate if it falls out of the grid boundaries

    // Index for T:  0 <= it < npoints-2
    m_pLower = std::lower_bound(mp_t, mp_t+m_npoints, (float)T);

    if (m_pLower == mp_t+m_npoints)
        m_it = m_npoints-2;
    else
        m_it = (mp_t == m_pLower ? 0 : std::distance(mp_t, m_pLower)-1);

    // Step 2: Compute the linear interpolation
    m_t1 = mp_t[m_it];
    m_t2 = mp_t[m_it+1];

    m_w1 = (m_t2-T)/(m_t2-m_t1);
    m_w2 = (T-m_t1)/(m_t2-m_t1);

    m_ndata = m_nparams*m_nbands;
    m_p1 = mp_params + m_it*m_ndata;
    m_p2 = m_p1 + m_ndata;


    // Log-linear interpolation if non zero values
    for (int i = m_nparams*band; i < m_nparams*(band+1); ++i) {
        if (m_p1[i] != 0.0 & m_p2[i] != 0.0) {
            p_params[i-m_nparams*band] = exp(m_w1*log(m_p1[i]) + m_w2*log(m_p2[i]));
        } else {
            p_params[i-m_nparams*band] = m_w1*m_p1[i] + m_w2*m_p2[i];
        }
    }

    // Enforces nul emission coefficients at low temperature
    // Extrapolation is not sufficiently accurate when the chi factor is taken into account
    if (T<m_t1 & m_nparams>2) {
          p_params[1] = 0.0;
          p_params[2] = 0.0;
    }
}

double SnbContinuumSystem::getThinAbsorptionSingleBand(ThermoData& thermo, CFuint cellID, double band)
{
    m_ie = thermo.speciesIndex("e-");


    m_pParamsCurrentCell = new double[m_nparams];

        thermo.setState(cellID);

        m_tr  = thermo.Tr();
        m_tv  = thermo.Tv();
        m_nd  = thermo.N();

//         std::cout << " I " << i << " tr " << tr << " tv " << tv << " nd " << nd << std::endl;

        // Computation of the negative ion concentrations if not specified
        if (m_system == "N-_bf" || m_system == "O-_bf") {
            if (thermo.N()[m_species_index] == 0.0) {
                double N = min( mp_chi->concNegIon(thermo) , thermo.N(thermo.speciesIndex(m_species.substr(0,1))));
                m_pa = N * KB * m_tr;
            } else {
                m_pa  = thermo.N(m_species_index) * KB * m_tr;
            }
        } else {
            m_pa  = thermo.N(m_species_index) * KB * m_tr;
        }
        m_pe  = thermo.N(m_ie) * KB * m_tv;

        getParametersSingleBand(m_tv, band, m_pParamsCurrentCell);

        switch (m_sort) {

        case BOUNDFREE :
            if (m_nparams > 2) {
                // Computation of the chi factor
                m_chi_factor = mp_chi->computeChi(thermo);

                // BF systems which can take into account chemical disequilibrium
                if (m_system == "N_bf" || m_system == "O_bf" || m_system == "C_bf"){
                    AtomicPartFunc qat(m_species, m_datadir);
                    double q = qat.Q(thermo);
                    double eion = mp_chi->ionizationEnergy(); // J/mol
                    double lowT_fac = std::exp(eion/(2.0*KB*m_tv));


                    return (m_pParamsCurrentCell[0] - m_chi_factor * m_pParamsCurrentCell[1] * lowT_fac) * m_pa / q * m_tv/ m_tr;

                } else if (m_system == "O2_bf_SR") {
                    // SR continuum: Tv for absorption, Tr for spontaneous and induced emission

                    double kappa = m_pParamsCurrentCell[0] *m_pa* m_tv / m_tr;
                    getParametersSingleBand(m_tr,band, m_pParamsCurrentCell);
                    kappa -=  m_chi_factor * m_pParamsCurrentCell[1] * m_pa;
                    return kappa;


                } else {
                    double eion = mp_chi->ionizationEnergy(); // J/mol
                    double lowT_fac = std::exp(eion/(2.0*KB*m_tv));


                    return (m_pParamsCurrentCell[0] - m_chi_factor * m_pParamsCurrentCell[1] * lowT_fac) * m_pa * m_tv / m_tr;

                }
            } else {
                // BF systems for which chemical equilibrium is always assumed

                return m_pParamsCurrentCell[0] *m_pa * m_tv / m_tr;

            }
            break;

        case FREEFREE :
            if (m_system != "C3_ff") {
                return m_pParamsCurrentCell[0] * m_pe * m_pa * m_tv / m_tr;
            } else {
                // C3_ff is a fake free-free system, it is actually treated
                // as thin but in thermal equilibrium so we are using this
                // trick for now
                return m_pParamsCurrentCell[0] * m_pa * m_tv / m_tr * 1.0e-4;

            }
            break;
        }


    delete [] m_pParamsCurrentCell;

}

double SnbContinuumSystem::getThinEmissionSingleBand(ThermoData& thermo,CFuint cellID, double band)
{

    const size_t ie = thermo.speciesIndex("e-");
    const double* nd;
    double pa, pe, tr, tv, chi_factor;
    double *p_params = new double [m_nparams];


    thermo.setState(cellID);

    tr  = thermo.Tr();
    tv  = thermo.Tv();
    nd  = thermo.N();

    //         std::cout << " I " << i << " tr " << tr << " tv " << tv << " nd " << nd << std::endl;

    // Computation of the negative ion concentrations if not specified
    if (m_system == "N-_bf" || m_system == "O-_bf") {
        if (thermo.N()[m_species_index] == 0.0) {
            double N = min( mp_chi->concNegIon(thermo) , thermo.N(thermo.speciesIndex(m_species.substr(0,1))));
            pa = N * KB * tr;
        } else {
            pa  = thermo.N(m_species_index) * KB * tr;
        }
    } else {
        pa  = thermo.N(m_species_index) * KB * tr;
    }
    pe  = thermo.N(ie) * KB * tv;

    getParametersSingleBand(tv, band, p_params);

    switch (m_sort) {

    case BOUNDFREE :
        if (m_nparams > 2) {
            // Computation of the chi factor
            chi_factor = mp_chi->computeChi(thermo);

            // BF systems which can take into account chemical disequilibrium
            if (m_system == "N_bf" || m_system == "O_bf" || m_system == "C_bf"){
                AtomicPartFunc qat(m_species, m_datadir);
                double q = qat.Q(thermo);
                double eion = mp_chi->ionizationEnergy(); // J/mol
                double lowT_fac = std::exp(eion/(2.0*KB*tv));


                return p_params[2] * chi_factor * lowT_fac * pa / q * tv/tr;

            } else if (m_system == "O2_bf_SR") {
                // SR continuum: Tv for absorption, Tr for spontaneous and induced emission

                return  p_params[2] *chi_factor * pa;

            } else {
                double eion = mp_chi->ionizationEnergy(); // J/mol
                double lowT_fac = std::exp(eion/(2.0*KB*tv));


                return p_params[2] * chi_factor * lowT_fac * pa * tv / tr;

            }
        } else {
            // BF systems for which chemical equilibrium is always assumed

            return p_params[1] *pa * tv / tr;


        }
        break;

    case FREEFREE :
        if (m_system != "C3_ff") {
            return p_params[1] * pe * pa * tv / tr;

        } else {
            // C3_ff is a fake free-free system, it is actually treated
            // as thin but in thermal equilibrium so we are using this
            // trick for now

            return p_params[1] * pa * tv / tr * 1.0e-4;

        }
        break;
    }


    delete [] p_params;
}

double SnbContinuumSystem::getEmissionSingleBandFixedCell(ThermoData& thermo, CFuint band)
{

    const double* nd;
    double tr, tv;

//    thermo.setState(cellID);

    tr  = thermo.Tr();
    tv  = thermo.Tv();
    nd  = thermo.N();
//    std::cout << "SnbContinuumSystem::getEmissionSingleBandFixedCell => Interpolate..." << std::endl;


    // Log-linear interpolation if non zero values
    for (int i = m_nparams*band; i < m_nparams*(band+1); ++i) {
        if (m_p1[i] != 0.0 & m_p2[i] != 0.0) {
            m_pBandParams[i-m_nparams*band] = exp(m_w1*log(m_p1[i]) + m_w2*log(m_p2[i]));
        } else {
            m_pBandParams[i-m_nparams*band] = m_w1*m_p1[i] + m_w2*m_p2[i];
        }
//        std::cout << "SnbContinuumSystem::getEmissionSingleBandFixedCell => i=" << i<< std::endl;
//        std::cout << "m_p1 " << m_p1[i] << std::endl;
//        std::cout << "m_p2 " << m_p2[i] << std::endl;
//        std::cout << "SnbContinuumSystem::getEmissionSingleBandFixedCell => crash?"<< std::endl;


    }
//    std::cout << "SnbContinuumSystem::getEmissionSingleBandFixedCell => Interpolated" << std::endl;

    // Enforces nul emission coefficients at low temperature
    // Extrapolation is not sufficiently accurate when the chi factor is taken into account
    if (tv <m_t1 & m_nparams>2) {
          m_pBandParams[1] = 0.0;
          m_pBandParams[2] = 0.0;
    }

    switch (m_sort) {

    case BOUNDFREE :
        if (m_nparams > 2) {
            // Computation of the chi factor
            // BF systems which can take into account chemical disequilibrium
            if (m_system == "N_bf" || m_system == "O_bf" || m_system == "C_bf"){
                return m_pBandParams[2] * m_chi_factor * m_lowT_fac * m_pa / m_q * tv/tr;
            } else if (m_system == "O2_bf_SR") {
                // SR continuum: Tv for absorption, Tr for spontaneous and induced emission
                return  m_pBandParams[2] * m_chi_factor * m_pa;
            } else {
                return m_pBandParams[2] * m_chi_factor * m_lowT_fac * m_pa * tv / tr;
            }
        } else {
            // BF systems for which chemical equilibrium is always assumed
            return m_pBandParams[1] *m_pa * tv / tr;
        }
        break;

    case FREEFREE :
        if (m_system != "C3_ff") {
            return m_pBandParams[1] * m_pe * m_pa * tv / tr;
        } else {
            // C3_ff is a fake free-free system, it is actually treated
            // as thin but in thermal equilibrium so we are using this
            // trick for now
            return m_pBandParams[1] * m_pa * tv / tr * 1.0e-4;
        }
        break;
    }


}

void SnbContinuumSystem::setupEmissionCoefficientsFixedCell(ThermoData &thermo, CFuint cellID)
{
    const size_t ie = thermo.speciesIndex("e-");
    const double* nd;
    double tr, tv;
//    double *p_params = new double [m_nparams];


    thermo.setState(cellID);

    tr  = thermo.Tr();
    tv  = thermo.Tv();
    nd  = thermo.N();

    // Computation of the negative ion concentrations if not specified
    if (m_system == "N-_bf" || m_system == "O-_bf") {
        if (thermo.N()[m_species_index] == 0.0) {
            double N = min( mp_chi->concNegIon(thermo) , thermo.N(thermo.speciesIndex(m_species.substr(0,1))));
            m_pa = N * KB * tr;
        } else {
            m_pa  = thermo.N(m_species_index) * KB * tr;
        }
    } else {
        m_pa  = thermo.N(m_species_index) * KB * tr;
    }
    m_pe  = thermo.N(ie) * KB * tv;



    if (m_sort==BOUNDFREE) {
        if (m_nparams > 2) {
            // Computation of the chi factor
            m_chi_factor = mp_chi->computeChi(thermo);

            // BF systems which can take into account chemical disequilibrium
            if (m_system == "N_bf" || m_system == "O_bf" || m_system == "C_bf"){
                AtomicPartFunc qat(m_species, m_datadir);
                m_q = qat.Q(thermo);
                m_eion = mp_chi->ionizationEnergy(); // J/mol
                m_lowT_fac = std::exp(m_eion/(2.0*KB*tv));
            } else if (m_system == "O2_bf_SR") {
                // Nothing to be precomputed here
            } else {
                m_eion = mp_chi->ionizationEnergy(); // J/mol
                m_lowT_fac = std::exp(m_eion/(2.0*KB*tv));
            }
        }
    }

}

double SnbContinuumSystem::setupParamInterpolation(ThermoData &thermo, CFuint cellID)
{
//    std::cout << "SnbCO2System::setupParamInterpolation => start" << std::endl;
    thermo.setState(cellID);
//    double T=thermo.Tr();
    //TODO double check
    double T=thermo.Tv();
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

    m_nParamdata = m_nparams*m_nbands;
//    std::cout << "SnbCO2System::setupParamInterpolation => m_it=" << m_it << " m_nParamdata=" << m_nParamdata << " m_nparams=" << m_nparams << " m_nbands=" << m_nbands << std::endl;
    m_p1 = mp_params + m_it*m_nParamdata;
    m_p2 = m_p1 + m_nParamdata;
//    std::cout << "SnbCO2System::setupParamInterpolation => end" << std::endl;


}

double SnbContinuumSystem::getParameterAt(double T, size_t band, size_t point) const
{
//    // Step 1: Determine the indices for T which either bound the
//    // interval in which T falls or the correct interval which will
//    // be used to extrapolate if it falls out of the grid boundaries

//    // Index for T:  0 <= it < npoints-2
//    float* p_lower = std::lower_bound(mp_t, mp_t+m_npoints, (float)T);
//    size_t it;
//    if (p_lower == mp_t+m_npoints)
//        it = m_npoints-2;
//    else
//        it = (mp_t == p_lower ? 0 : std::distance(mp_t, p_lower)-1);

//    // Step 2: Compute the linear interpolation
//    const double t1 = mp_t[it];
//    const double t2 = mp_t[it+1];

//    const double w1 = (t2-T)/(t2-t1);
//    const double w2 = (T-t1)/(t2-t1);

//    const size_t ndata = m_nparams*m_nbands;
//    const double* const p1 = mp_params + it*ndata;
//    const double* const p2 = p1 + ndata;



//    // Log-linear interpolation if non zero values
//    for (int i = 0; i < ndata; ++i) {
//        if (p1[i] != 0.0 & p2[i] != 0.0) {
//            p_params[i] = exp(w1*log(p1[i]) + w2*log(p2[i]));
//        } else {
//            p_params[i] = w1*p1[i] + w2*p2[i];
//        }

//    }

//    // Enforces nul emission coefficients at low temperature
//    // Extrapolation is not sufficiently accurate when the chi factor is taken into account
//    if (T<t1 & m_nparams>2) {
//      for (int i = 0; i < ndata; i=i+3) {
//          p_params[i+1] = 0.0;
//          p_params[i+2] = 0.0;
//      }
//    }

}

void SnbContinuumSystem::exportLocalParameters(boost::filesystem::path exportPath, ThermoData& thermo)
{
    if (thermo.usePrecomputedContinuumParameters()) {
    const size_t ndata = m_nbands*2;
    
    m_outFileHandle  = COOLFluiD::Environment::SingleBehaviorFactory<COOLFluiD::Environment::FileHandlerOutput>::getInstance().create();

    //Only for testing, in reality mp_locparams is setup
    //mp_locparams = new double [thermo.nCells()*ndata];

    boost::filesystem::path loc("loc"+m_species+m_system);
    boost::filesystem::path outfile = exportPath / loc;
    ofstream& fout = m_outFileHandle->open(outfile);

    fout << thermo.nCells() << std::endl;
    fout << m_nbands << std::endl;
    fout << 2 << std::endl;

    for (int i=0; i<thermo.nCells(); i++) {

        for (int b = 0; b < m_nbands; ++b) {
           fout << mp_locparams[i*ndata+b*2+0] <<" "<< mp_locparams[i*ndata+b*2+1] <<" ";
        }
        fout << "\n";
    }

    fout.close();
    }
}

void SnbContinuumSystem::readLocalParameters(boost::filesystem::path importPath)
{
    boost::filesystem::path loc("loc"+m_species+m_system);
    boost::filesystem::path infile = importPath / loc;
    std::vector<std::string> tokens;
    uint localIndex=0;

    m_inFileHandle  = COOLFluiD::Environment::SingleBehaviorFactory<COOLFluiD::Environment::FileHandlerInput>::getInstance().create();

    // Open the file
    std::ifstream& fin = m_inFileHandle->open(infile);

    uint numberCells;
    uint numberParameters;


    if (fin.is_open()) {

        std::string line;

        //Load number of cells
        std::getline(fin, line);
        numberCells=atoi(line.c_str());

        //Load number of bands
        std::getline(fin, line);
        m_nbands=atoi(line.c_str());


        //Load number of parameters
        std::getline(fin, line);
        numberParameters=atoi(line.c_str());

        mp_locparams = new double[numberCells*m_nbands*numberParameters];

        std::getline(fin, line);

//        std::cout << "Cells: " <<  numberCells << "\n" << "Bands: " << m_nbands << "\n" <<
//                "Parameters: " << numberParameters << "\n";


        while (!line.empty()) {

            String::tokenize(line, tokens, " \t\r\n");

            for (int i = 0; i < tokens.size(); ++i) {
                mp_locparams[localIndex]=atof(tokens[i].c_str());
//                std::cout << mp_locparams[localIndex] << " ";
                localIndex++;
            }

//            std::cout << " \n";

            std::getline(fin, line);

        }
    }

    else {std::cout << "locFile not open";}

    m_inFileHandle->close();
}

void SnbContinuumSystem::setupLocalParameters(ThermoData& thermo)
{
    if (thermo.usePrecomputedContinuumParameters()) {
//        std::cout <<"CONTINUUM USE PRECOMPUTED" << std::endl;

        const size_t ndata = m_nbands*2;
        const size_t ie = thermo.speciesIndex("e-");
        const double* nd;
        double pa, pe, tr, tv, chi_factor;
        double *p_params = new double [m_nparams*m_nbands];

        // Allocation of the local parameter matrix
        if (mp_locparams) {
            delete [] mp_locparams;
            mp_locparams = NULL;
        }
        mp_locparams = new double [thermo.nCells()*ndata];
        m_nbParams= thermo.nCells()*ndata;

        for (int i = 0; i<thermo.nCells(); i++) {

            thermo.setState(i);

            tr  = thermo.Tr();
            tv  = thermo.Tv();
            nd  = thermo.N();

            //         std::cout << " I " << i << " tr " << tr << " tv " << tv << " nd " << nd << std::endl;

            // Computation of the negative ion concentrations if not specified
            if (m_system == "N-_bf" || m_system == "O-_bf") {
                if (thermo.N()[m_species_index] == 0.0) {
                    double N = min( mp_chi->concNegIon(thermo) , thermo.N(thermo.speciesIndex(m_species.substr(0,1))));
                    pa = N * KB * tr;
                } else {
                    pa  = thermo.N(m_species_index) * KB * tr;
                }
            } else {
                pa  = thermo.N(m_species_index) * KB * tr;
            }
            pe  = thermo.N(ie) * KB * tv;

            getParameters(tv, p_params);

            switch (m_sort) {

            case BOUNDFREE :
                if (m_nparams > 2) {
                    // Computation of the chi factor
                    chi_factor = mp_chi->computeChi(thermo);

                    // BF systems which can take into account chemical disequilibrium
                    if (m_system == "N_bf" || m_system == "O_bf" || m_system == "C_bf"){
                        AtomicPartFunc qat(m_species, m_datadir);
                        double q = qat.Q(thermo);
                        double eion = mp_chi->ionizationEnergy(); // J/mol
                        double lowT_fac = std::exp(eion/(2.0*KB*tv));

#pragma omp parallel for
                        for (int b = 0; b < m_nbands; ++b) {
                            mp_locparams[i*ndata+b*2+0] = (p_params[b*m_nparams+0]
                                    - chi_factor * p_params[b*m_nparams+1] * lowT_fac) * pa / q * tv/ tr;
                            mp_locparams[i*ndata+b*2+1] = p_params[b*m_nparams+2] * chi_factor * lowT_fac * pa / q * tv/tr;
                        }
                    } else if (m_system == "O2_bf_SR") {
                        // SR continuum: Tv for absorption, Tr for spontaneous and induced emission
#pragma omp parallel for
                        for (int b = 0; b < m_nbands; ++b)
                            mp_locparams[i*ndata+b*2+0] = p_params[b*m_nparams+0] *pa* tv / tr;
                        getParameters(tr, p_params);
#pragma omp parallel for
                        for (int b = 0; b < m_nbands; ++b) {
                            mp_locparams[i*ndata+b*2+0] -=  chi_factor * p_params[b*m_nparams+1] * pa;
                            mp_locparams[i*ndata+b*2+1] = p_params[b*m_nparams+2] *chi_factor * pa;
                        }
                    } else {
                        double eion = mp_chi->ionizationEnergy(); // J/mol
                        double lowT_fac = std::exp(eion/(2.0*KB*tv));

#pragma omp parallel for
                        for (int b = 0; b < m_nbands; ++b) {
                            mp_locparams[i*ndata+b*2+0] = (p_params[b*m_nparams+0]
                                    - chi_factor * p_params[b*m_nparams+1] * lowT_fac) * pa * tv / tr;
                            mp_locparams[i*ndata+b*2+1] = p_params[b*m_nparams+2] * chi_factor * lowT_fac * pa * tv / tr;
                        }
                    }
                } else {
                    // BF systems for which chemical equilibrium is always assumed
#pragma omp parallel for
                    for (int b = 0; b < m_nbands; ++b) {
                        mp_locparams[i*ndata+b*2+0] = p_params[b*m_nparams+0] *pa * tv / tr;
                        mp_locparams[i*ndata+b*2+1] = p_params[b*m_nparams+1] *pa * tv / tr;
                    }

                }
                break;

            case FREEFREE :
                if (m_system != "C3_ff") {
#pragma omp parallel for
                    for (int b = 0; b < m_nbands; ++b) {
                        mp_locparams[i*ndata+b*2+0]
                                = p_params[b*m_nparams+0] * pe * pa * tv / tr;
                        mp_locparams[i*ndata+b*2+1]
                                = p_params[b*m_nparams+1] * pe * pa * tv / tr;
                    }
                } else {
                    // C3_ff is a fake free-free system, it is actually treated
                    // as thin but in thermal equilibrium so we are using this
                    // trick for now
#pragma omp parallel for
                    for (int b = 0; b < m_nbands; ++b) {
                        mp_locparams[i*ndata+b*2+0]
                                = p_params[b*m_nparams+0] * pa * tv / tr * 1.0e-4;
                        mp_locparams[i*ndata+b*2+1]
                                = p_params[b*m_nparams+1] * pa * tv / tr * 1.0e-4;
                    }
                }
                break;
            }
        }

        delete [] p_params;
    }
}

double SnbContinuumSystem::getLocalParameter(const int& i, const int& j, const int& k) const
// Cell (i), Band (j) and Parameter (k) indices
{
    return mp_locparams[i*m_nbands*2+j*2+k];
}

double SnbContinuumSystem::emittedPower(ThermoData& thermo,int i)
{
    double sum = 0.0;

    if (thermo.usePrecomputedContinuumParameters()) {
        for (int b = 0; b < m_nbands; ++b) {
            sum += mp_locparams[i*m_nbands*2 + b*2 + 1];
        }
    }
    else {
//        for (int b = 0; b < m_nbands; ++b) {
//            sum += getThinEmissionSingleBand(thermo,i,b);
//        }


        setupEmissionCoefficientsFixedCell(thermo,i);
        setupParamInterpolation(thermo,i);
        m_pBandParams = new double [m_nparams];

//        std::cout << "SnbContinuumSystem::emittedPower => m_nparams=" << m_nparams << std::endl;

        for (int b = 0; b < m_nbands; ++b) {
//            std::cout << "SnbContinuumSystem::emittedPower => b=" << b << std::endl;
            sum += getEmissionSingleBandFixedCell(thermo, b);
        }

        delete [] m_pBandParams;
    }


//    CFLog(VERBOSE, "SnbContinuumSystem::emittedPower( " << i << "), " << this->systemName() << "=" << sum*1000.0<< "\n");
    return sum*1000.0;
}

double SnbContinuumSystem::emittedPowerSaveMemory(int i)
{
    double sum = 0.0;
    for (int b = 0; b < m_nbands; ++b) {
        sum += mp_locparams[i*m_nbands*2 + b*2 + 1];
    }
    return sum*1000.0;
}

double SnbContinuumSystem::localParamsInPlace(double &localEmissivity, int cellID, int band, ThermoData &thermo)
{
//    const size_t ndata = m_nbands*2;
//    const size_t ie = thermo.speciesIndex("e-");
//    const double* nd;
//    double pa, pe, tr, tv, chi_factor;
//    double *p_params = new double [m_nparams*m_nbands];

//    // Allocation of the local parameter matrix
//    if (mp_locparams) {
//        delete [] mp_locparams;
//        mp_locparams = NULL;
//    }
//    mp_locparams = new double [thermo.nCells()*ndata];


//    thermo.setState(cellID);

//    tr  = thermo.Tr();
//    tv  = thermo.Tv();
//    nd  = thermo.N();

//    std::cout << " I " << cellID << " tr " << tr << " tv " << tv << " nd " << nd << std::endl;

//    // Computation of the negative ion concentrations if not specified
//    if (m_system == "N-_bf" || m_system == "O-_bf") {
//        if (thermo.X()[m_species_index] == 0.0) {
//            double N = min( mp_chi->concNegIon(thermo) , thermo.N(thermo.speciesIndex(m_species.substr(0,1))));
//            pa = N * KB * tr;
//        } else {
//            pa  = thermo.N(m_species_index) * KB * tr;
//        }
//    } else {
//        pa  = thermo.N(m_species_index) * KB * tr;
//    }
//    pe  = thermo.N(ie) * KB * tv;

//    getParameters(tv, p_params);

//    switch (m_sort) {

//        case BOUNDFREE :
//            if (m_nparams > 2) {
//                // Computation of the chi factor
//                chi_factor = mp_chi->computeChi(thermo);

//                // BF systems which can take into account chemical disequilibrium
//                if (m_system == "N_bf" || m_system == "O_bf" || m_system == "C_bf"){
//                    AtomicPartFunc qat(m_species);
//                    double q = qat.Q(thermo);
//                    double eion = mp_chi->ionizationEnergy(); // J/mol
//                    double lowT_fac = std::exp(eion/(2.0*KB*tv));


//                    localEmissivity = p_params[b*m_nparams+2] * chi_factor * lowT_fac * pa / q * tv/tr;

//                } else if (m_system == "O2_bf_SR") {
//                    // SR continuum: Tv for absorption, Tr for spontaneous and induced emission

//                    getParameters(tr, p_params);
//                    localEmissivity = p_params[b*m_nparams+2] *chi_factor * pa;

//                } else {
//                    double eion = mp_chi->ionizationEnergy(); // J/mol
//                    double lowT_fac = std::exp(eion/(2.0*KB*tv));

//                    localEmissivity = p_params[b*m_nparams+2] * chi_factor * lowT_fac * pa * tv / tr;

//                }
//            } else {
//                // BF systems for which chemical equilibrium is always assumed
//                localEmissivity = p_params[b*m_nparams+1] *pa * tv / tr;


//            }
//            break;

//        case FREEFREE :
//            if (m_system != "C3_ff") {


//                    localEmissivity
//                      = p_params[b*m_nparams+1] * pe * pa * tv / tr;

//            } else {
//                // C3_ff is a fake free-free system, it is actually treated
//                // as thin but in thermal equilibrium so we are using this
//                // trick for now

//                    localEmissivity
//                      = p_params[b*m_nparams+1] * pa * tv / tr * 1.0e-4;

//            }
//            break;
//    }
}

/// Computes the emission for each band
void SnbContinuumSystem::bandEmission(ThermoData& thermo, int i, double* const p_emis)
{
    if (thermo.usePrecomputedContinuumParameters()) {
        for (int b = 0; b < m_nbands; ++b) {
            p_emis[b] = mp_locparams[i*m_nbands*2 + b*2 + 1] * 1000.0;
        }
    }
    else {
//        for (int b = 0; b < m_nbands; ++b) {
//            p_emis[b] = getThinEmissionSingleBand(thermo,i,b) * 1000.0;
//        }

        setupEmissionCoefficientsFixedCell(thermo,i);
        setupParamInterpolation(thermo,i);
        m_pBandParams = new double [m_nparams];

        for (int b = 0; b < m_nbands; ++b) {
            p_emis[b] = getEmissionSingleBandFixedCell(thermo, b) * 1000.0;
        }

        delete [] m_pBandParams;

    }
}




double SnbContinuumSystem::opticalThickness(
    const PhotonPath& path, int ic, double sig)
{
    // Find the band corresponding to this wavelength
    int b = int(sig / 1000.0);

    // Convert band to local indexing
    if (b < lowBand() || b > highBand())
        return 0.0;
    b -= lowBand();

    double ku = 0.0;
    for (int i = 0; i < ic; i++)
        ku += mp_locparams[path.cellID(i)*m_nbands*2 + b*2 + 0] * path.cellDistance(i);

    return ku;
}

double SnbContinuumSystem::opticalThickness(const HSNBNonThickParameterSet &pathParams)
{
    CFLog(DEBUG_MED, "SnbContinuumSystem::opticalThickness => opticalThickness="<< pathParams.kappa << "\n");
    return pathParams.kappa;
}

void SnbContinuumSystem::addStateParams(ThermoData& thermo, HSNBNonThickParameterSet &nonThickParams, CFreal cellDistance, CFuint localCellID, CFreal sig)
{
//    CFLog(VERBOSE, "SnbContinuumSystem::addStateParams, " << this->speciesName() << " / " << this->systemName() << "  \n");
    int b = int(sig / 1000.0);

    // Convert band to local indexing
    if (b < lowBand() || b > highBand()) {
//        CFLog(VERBOSE, "SnbContinuumSystem::addStateParams => b=" << b << " out of bounds. \n");
        nonThickParams.addState(0.0);
        return;
    }
    b -= lowBand();

    if (nonThickParams.isEmpty()) {

        if (thermo.usePrecomputedContinuumParameters()) {
            m_tempKu=mp_locparams[localCellID*2*m_nbands + b*2] * cellDistance;
        }
        else {
            m_tempKu=getThinAbsorptionSingleBand(thermo,localCellID,b)*cellDistance;
        }

//        CFLog(VERBOSE, "SnbContinuumSystem::addStateParams =>NONTHICK: mp_locparams["<< localCellID <<"*2*"<< m_nbands<<" + "<< b << "*2]=" << mp_locparams[localCellID*2*m_nbands + b*2] << "\n");
//        CFLog(VERBOSE, "SnbContinuumSystem::addStateParams =>Empty ParamSet: Add Ku=" << m_tempKu <<  " at localID " << localCellID  << "\n");

        nonThickParams.addState(m_tempKu);
    }
    else {
        //ku=sum_i<=j(ku)
        if (thermo.usePrecomputedContinuumParameters()) {
            m_tempKu=nonThickParams.kappa+mp_locparams[localCellID*2*m_nbands + b*2] * cellDistance;
        }
        else {
            m_tempKu=nonThickParams.kappa+getThinAbsorptionSingleBand(thermo,localCellID,b)*cellDistance;
        }

//        CFLog(VERBOSE, "SnbContinuumSystem::addStateParams =>NONTHICK: mp_locparams["<< localCellID <<"*2*"<< m_nbands<<" + "<< b << "*2]=" << mp_locparams[localCellID*2*m_nbands + b*2] << "\n");
//        CFLog(VERBOSE, "SnbContinuumSystem::addStateParams =>NonEmpty ParamSet: Add Ku=" << m_tempKu <<  " at localID " << localCellID  << "\n");

        nonThickParams.addState(m_tempKu);
    }


}

void SnbContinuumSystem::determineBandRange()
{
    // Lower bound, start from 1000
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
//    cout << "Original band: " << m_band1 << " " << m_bandn <<", mnbands=" << m_nbands << endl;
}

bool SnbContinuumSystem::waveNumberIsInBandRange(CFuint waveNumber)
{
    return ((waveNumber>=lowBand()) && (waveNumber<=highBand()));
}

void SnbContinuumSystem::setupProducts()
{

    if(m_sort == BOUNDFREE) {

        if (m_system == "N-_bf") { m_products = "N"; }
        else if (m_system == "O-_bf") { m_products = "O"; }
        else if (m_system == "N_bf") { m_products = "N+"; }
        else if (m_system == "O_bf") { m_products = "O+"; }
        else if (m_system == "N2_bf") { m_products = "N2+"; }
        else if (m_system == "O2_bf") { m_products = "O2+"; }
        else if (m_system == "NO_bf") { m_products = "NO+"; }
        else if (m_system == "O2_bf_SR") { m_products = "O"; }
        else if (m_system == "C-_bf") { m_products = "C"; }
        else if (m_system == "C_bf") { m_products = "C+"; }
        else if (m_system == "C2_bf") { m_products = "C2+"; }
        else if (m_system == "CN_bf") { m_products = "CN+"; }
        else if (m_system == "CH_bf") { m_products = "CH+"; }
        else if (m_system == "CO_bf") { m_products = "CO+"; }
        else if (m_system == "H_bf") { m_products = "H+"; }
        else if (m_system == "H2_bf") { m_products = "H2+"; }
        else { 
           m_products = ""; 
           CFLog(ERROR, "SnbContinuumSystem::setupProducts => WARNING: The system " << m_system << " is not supported ! \n");
        }

//        int index = ThermoData::getInstance().speciesIndex(m_products);
//        if(index < 0 || index >= ThermoData::getInstance().nSpecies())
//           cout << "WARNING: the species produced by the system "
//                << m_system << " is not in the species list !" << endl;

    } else {
        m_products = "";
    }

}

void SnbContinuumSystem::loadTemperatureGrid()
{
    // Open the file corresponding to the first band
    std::ifstream file((m_directory+"/"+bandFilename(m_band1)).c_str());
    
    // Skip first two lines
    std::string line;
    std::getline(file, line);
    std::getline(file, line);
    
    // Now read the first temperature from the file and count the number of
    // parameters (extra columns)
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


void SnbContinuumSystem::loadBandParameters(const size_t& iband)
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

