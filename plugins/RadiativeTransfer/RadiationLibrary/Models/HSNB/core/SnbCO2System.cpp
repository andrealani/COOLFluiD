#include "SnbCO2System.h"

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

SnbCO2System::SnbCO2System(const std::string& directory, const std::string& nup, const ThermoData& thermo)
    : m_directory(directory)
{
    m_species_index = thermo.speciesIndex("CO2");
    
    // Add the full path to the default data directory
    if (m_directory == "CO2") {
        m_directory = std::getenv("HTGR_DATA_DIRECTORY");
        if (m_directory[m_directory.size()-1] != '/')
            m_directory += '/';
        m_directory += "modele_msbe/CO2";
    }

    if (m_species_index < 0) {
        cout << "Error loading CO2 system ! "
             << "The species isn't loaded." << endl;
        exit(1);
    }

    // Determine the min and max band
    determineBandRange();
    
    // Load the temperature grid and determine number of parameters
    loadTemperatureGrid();

    if (nup == "LS")
        m_nup = LINDQUIST_SIMMONS;
    else
        m_nup = CURTIS_GODSON;

    m_lorentz = false;

    cout << "Loading CO2 system "
         << " (" << m_band1_up << " - " << m_bandn_up << ") "
         << ", Non uniform path approximation = " 
         << (m_nup == CURTIS_GODSON ? "CG" : "LS" ) << endl;
    
    // Allocate storage for the parameter information
    mp_params = new double [m_npoints * m_nbands_down * m_nparams];
    mp_locparams_up = NULL;
    mp_locparams_down = NULL;
    
    // Now load the parameter information
    #pragma omp parallel for
    for (size_t i = 0; i < m_nbands_down; ++i)
        loadBandParameters(i);
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

    const size_t ndata = m_nparams*m_nbands_down;
    const double* const p1 = mp_params + it*ndata;
    const double* const p2 = p1 + ndata;

    for (int i = 0; i < ndata; ++i) {
        p_params[i] = w1*p1[i] + w2*p2[i];
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
    double ptot, t, tv, pco2, xco2, pmoy;
    double *p_params = new double [m_nparams*m_nbands_down];
    const size_t ndata = m_nbands_down*4;

    if (mp_locparams_up) {
        delete [] mp_locparams_up;
    }
    if (mp_locparams_down) {
        delete [] mp_locparams_down;
    }

    // Up local parameters (Planck function over large bands)
    mp_locparams_up = new double [thermo.nCells()*m_nbands_up];
    for (int i=0; i<thermo.nCells(); i++) {
        thermo.setState(i);

        t    = thermo.Tr();
        tv   = thermo.Tv();
        if (tv != t)
           cout << "WARNING: Thermal non equilibrium is not supported for CO2 !" << endl;

        for (int b = 0; b < m_nbands_up; ++b)
            mp_locparams_up[i*m_nbands_up+b] = eqInt(b*1000.0+500.0, t);
    }

    // Down local parameters (transmissivity calculations)
    mp_locparams_down = new double [thermo.nCells()*ndata];
    for (int i=0; i<thermo.nCells(); i++) {
        thermo.setState(i);

        ptot = thermo.P();
        t    = thermo.Tr();
        pco2 = thermo.N(m_species_index) * KB * t;
        xco2 = thermo.X(m_species_index);
        
        getParameters(t, p_params);

        #pragma omp parallel for
        for (int b = 0; b < m_nbands_down; ++b) {
            mp_locparams_down[i*ndata+b*4+0] = p_params[b*m_nparams+0] * pco2; // kappa doppler
            mp_locparams_down[i*ndata+b*4+1] = p_params[b*m_nparams+3] * pco2; // kappa lorentz
            mp_locparams_down[i*ndata+b*4+2] = p_params[b*m_nparams+2] ; // beta doppler
            mp_locparams_down[i*ndata+b*4+3] = ( p_params[b*m_nparams+5] * xco2
                                               + p_params[b*m_nparams+6] ) * ptot ; // beta lorentz
        }

//        pmoy += ptot * (thermo.loc(i+1)-field.loc(i));
    }

//    //NOTEJB: ??
//    pmoy = pmoy / (field.loc(field.nPoints()-1)-field.loc(0));
//    if (pmoy >= 1000.0) m_lorentz = true;

    delete [] p_params;
}

double SnbCO2System::getLocalParameter(const int& i, const int& j, const int& k) const
// Cell (i), Band (j) and Parameter (k) indices
{
    return mp_locparams_up[i*m_nbands_up+j];
}

void SnbCO2System::determineBandRange()
{
    // Up SNB common to all systems : bands of 1000 cm-1
    m_band1_up = 0;
    m_bandn_up = 8;
    m_nbands_up = m_bandn_up - m_band1_up + 1;

    // Down SNB for CO2 : bands of 25 cm-1
    m_band1_down = 10;
    m_bandn_down = 332;
    m_nbands_down = m_bandn_down - m_band1_down + 1;
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

void SnbCO2System::downToUpSnb(const double* const p_down, double* const p_up)
{
    int j_up = 0;
    double sigjp1_down, sigjp1_up = 1000.0;
    double delta_nu = 0.025;

    for (int j_down=0; j_down<m_nbands_down; j_down++) {
        sigjp1_down = (j_down+10)*25+12.5;

        if (sigjp1_down > sigjp1_up) { 
            p_up[j_up] += p_down[j_down]*delta_nu/2.0;
            ++j_up;
            sigjp1_up +=1000.0;
            p_up[j_up] += p_down[j_down]*delta_nu/2.0;
        } else { 
            p_up[j_up] += p_down[j_down]*delta_nu;
        }
    }
    // Rescale first and last up band intervals
    p_up[0] /= 0.7625;
    p_up[m_nbands_up-1] /= 0.3125;

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

