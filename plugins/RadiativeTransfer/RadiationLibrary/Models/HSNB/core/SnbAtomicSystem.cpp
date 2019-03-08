#include "SnbAtomicSystem.h"
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cstdlib>
#include <cmath>

#include "StringUtils.h"

using namespace std;
using namespace COOLFluiD::RadiativeTransfer;

void swap(SnbAtomicSystem& s1, SnbAtomicSystem& s2)
{
    std::swap(s1.m_systems, s2.m_systems);
    std::swap(s1.m_specdir, s2.m_specdir);
    std::swap(s1.mp_grid, s2.mp_grid);
    std::swap(s1.mp_bounds, s2.mp_bounds);
    std::swap(s1.m_nloclblparams, s2.m_nloclblparams);
    std::swap(s1.mp_loclblparams, s2.mp_loclblparams);
    std::swap(s1.m_nlocparams, s2.m_nlocparams);
    std::swap(s1.mp_locparams, s2.mp_locparams);
    std::swap(s1.m_species_index, s2.m_species_index);
    std::swap(s1.m_nsystems, s2.m_nsystems);
    std::swap(s1.m_band1, s2.m_band1);
    std::swap(s1.m_bandn, s2.m_bandn);
    std::swap(s1.m_nbands, s2.m_nbands);
    std::swap(s1.m_computelbl, s2.m_computelbl);
    std::swap(s1.m_spectra_time, s2.m_spectra_time);
}

SnbAtomicSystem::SnbAtomicSystem(const std::string& spectradir, 
                                 const std::string& computespectra, 
                                 const std::string& grid_type,
                                 const vector<std::string>& systemlist, ThermoData& thermo)
  : m_nsystems(systemlist.size()),
    m_specdir(spectradir),
    mp_grid(NULL),
    mp_bounds(NULL),
    m_nloclblparams(0),
    mp_loclblparams(NULL),
    m_nlocparams(0),
    mp_locparams(NULL),
    m_spectra_time(0)
{
    // Load the atomic systems
    bool test_qss;
    string system_name;
    for (int i = 0; i < systemlist.size(); i++) {
      test_qss = false;
      system_name = systemlist[i];
      if(system_name[system_name.size()-1] == '*') {
          test_qss = true;
          system_name.erase(system_name.size()-1);
      }
      m_systems.push_back(LblAtomicSystem(thermo, system_name, m_datadir, test_qss));
      m_species_index.push_back(thermo.speciesIndex(system_name));
      if (m_species_index[i] < 0) {
          cout << "Error loading SNB system '" << m_systems[i].name() << "'! "
               << "The species isn't loaded for this system." << endl;
          exit(1);
      }
    }

    // Setup the lbl spectral grid
    if(computespectra == "false") { 
        m_computelbl = false; 
        cout << "Read precomputed lbl spectra for atoms" << endl;

        cout << "Reading lbl spectral grid..." << endl;
        mp_grid = new LblSpectralGrid(m_specdir);
    }   
    else { 
        m_computelbl = true; 
        cout << "Compute lbl spectra for atoms" << endl;

        cout << "Generating lbl spectral grid..." << endl;
        if (grid_type == "HiRes") {
            cout << "  using high resolution grid" << endl;
            mp_grid = new LblSpectralGrid(1.0e3,2.0e5);
        } else {
            std::vector<AtomicLineData> lines;
            for (int i = 0; i < m_systems.size(); ++i)
                m_systems[i].getLineData(thermo, lines);
            if (grid_type == "daSilva") {
                cout << "  using da Silva grid" << endl;
                mp_grid = new LblSpectralGrid(lines, 1.0e3, 2.0e5, DA_SILVA);
            } else {
                cout << "  using adaptive grid" << endl;
                mp_grid = new LblSpectralGrid(lines, 1.0e3, 2.0e5, SMART);
            }
        }

        cout << "  points: " << mp_grid->size() << "\n";
        cout << "  min [1/cm]: " << mp_grid->min() << "\n";
        cout << "  max [1/cm]: " << mp_grid->max() << "\n" << endl;

        ofstream gf((m_specdir + "/spectral-grid.dat").c_str(), ios::binary);
        size_t np = mp_grid->size();
        gf.write( (char*) &np, sizeof(size_t));
        gf.write( (char*) mp_grid->ptr(), mp_grid->size()*sizeof(double));
        gf.close();
    }

    determineBandRange();

    lblToSnbInit();

    cout << "Loading atomic systems ";
    for (int i=0; i < m_systems.size(); i++) {
      cout << m_systems[i].name() << ", ";
    }
    cout << " (" << m_band1 << " - " << m_bandn << ") " << endl;
    
}

SnbAtomicSystem::~SnbAtomicSystem()
{
    if (mp_grid != NULL)         delete mp_grid;
    if (mp_bounds != NULL)       delete [] mp_bounds;
    if (mp_loclblparams != NULL) delete [] mp_loclblparams;
    if (mp_locparams != NULL)    delete [] mp_locparams;
}

double SnbAtomicSystem::getLocalParameter(const int& i, const int& j, const int& k) const
// Cell (i), Band (j) and Parameter (k) indices
{
    return mp_locparams[i*m_nbands+j];
}

void SnbAtomicSystem::setupLocalParameters(ThermoData &thermo)
{
    // Spectra allocation
    int size = thermo.nCells()*mp_grid->size()*2;
    if (size > m_nloclblparams) {
        if (mp_loclblparams != NULL) delete [] mp_loclblparams;
        mp_loclblparams = new double [size];
    }
    m_nloclblparams = size;

    size = thermo.nCells()*m_nbands;
    if (size > m_nlocparams) {
        if (mp_locparams != NULL) delete [] mp_locparams;
        mp_locparams = new double [size];
    }
    m_nlocparams = size;

    // Initialize to zero
    fill(mp_loclblparams, mp_loclblparams+m_nloclblparams, 0.0);
    fill(mp_locparams, mp_locparams+m_nlocparams, 0.0);

    double* p_spectra = new double [mp_grid->size()];

    for (int i=0; i<thermo.nCells(); i++) {

        // Read or compute the lbl spectra according to the input requirements
        if (!m_computelbl) {
            readLblSpectra(i, mp_loclblparams + mp_grid->size()*2*i);
        } else {
            //computeLblSpectra(i, field, mp_loclblparams + mp_grid->size()*2*i);
            computeLblSpectraEq(i, thermo, mp_loclblparams + mp_grid->size()*2*i);
        }

        #pragma omp parallel for
        for (int j=0; j < mp_grid->size(); j++) {
            if (mp_loclblparams[mp_grid->size()*2*i+2*j+1] != 0.0) {
                p_spectra[j] = mp_loclblparams[mp_grid->size()*2*i+2*j+0]
                             / mp_loclblparams[mp_grid->size()*2*i+2*j+1];
            } else { p_spectra[j] = 0.0; }
            if (p_spectra[j] > 1.0e100) cout << "ETAKAPPA " << mp_loclblparams[mp_grid->size()*2*i+2*j+0]
                                        << " " << mp_loclblparams[mp_grid->size()*2*i+2*j+1]
                                        << " " << i << " " << j << endl;
        }

        lblToSnb(p_spectra, mp_locparams + m_nbands*i);

    }

    for (int i=0; i< thermo.nCells()*m_nbands-1; i++) {
        thermo.setState(i);
        if (mp_locparams[i] > 1.0e10) cout << "PROBLEME" << i << endl;
    }

    delete [] p_spectra;

}

void SnbAtomicSystem::readLblSpectra(int index, double* const p_spectra)
{
    // Define spectrum file name 
    ostringstream filename;
    filename << "spectra";
    filename << setfill('0') << setw(3) << index;
    filename << ".dat";
    cout << "Reading lbl spectrum from " + filename.str() + "..." << endl;

    string path = m_specdir + "/" + filename.str();
//    ifstream f( path.c_str() );
//    for (int k = 0; k < mp_grid->size(); ++k) {
//        f >> p_spectra[k*2] >> p_spectra[k*2+1];
//    }
    ifstream f( path.c_str() , ios::binary);
    f.read( (char*)p_spectra, 2*mp_grid->size()*sizeof(double));
    f.close();

}

void SnbAtomicSystem::writeLblSpectra(int index, const double* const p_spectra)
{
    // Define spectrum file name 
    ostringstream filename;
    filename << "spectra";
    filename << setfill('0') << setw(3) << index;
    filename << ".dat";
    cout << "Writing lbl spectrum to " + filename.str() + "..." << endl;

    string path = m_specdir + "/" + filename.str();
    ofstream f( path.c_str() );
    f.precision(8);
    f << scientific;
    for (int k = 0; k < mp_grid->size(); ++k) {
        f << (*mp_grid)[k] << " " << p_spectra[k*2] << " " << p_spectra[k*2+1] << "\n";
    }
//    ofstream f( path.c_str() , ios::binary);
//    f.write( (char*)p_spectra, 2*mp_grid->size()*sizeof(double) );
//    f.close();

}

void SnbAtomicSystem::computeLblSpectra(int index, ThermoData& thermo,
     double* const p_spectra)
{

    cout << "Computing lbl spectrum at point " << index << endl;

    // Set thermo properties according to Field Data at cell index index
    
    thermo.setState(index);

    // Call spectra function
    fill(p_spectra, p_spectra+mp_grid->size()*2, 0.0);
    for (int k = 0; k < m_systems.size(); ++k)
        m_systems[k].spectra(thermo, (*mp_grid), p_spectra, ETA_KAPPA);

    // Write computed spectra
    writeLblSpectra(index, p_spectra);

}

// Computation of the spectra satisfying equilibirum relations
void SnbAtomicSystem::computeLblSpectraEq(int index, ThermoData& thermo,
     double* const p_spectra)
{
    cout << "Computing lbl spectrum at point " << index << endl;

    // Set thermo properties according to Field Data at cell index index  
    thermo.setState(index);

    clock_t t1 = clock();
    // Call spectraEq function (only eta/sigma is computed)
    fill(p_spectra, p_spectra+mp_grid->size()*2, 0.0);
    for (int k = 0; k < m_systems.size(); ++k)
        m_systems[k].spectraEq(thermo, *mp_grid, p_spectra);
    std::cout << "Hello" << endl;

    // Compute eta and kappa such that equilibirum is retrieved
    double sig;
    const double fac2 = HP*C0*100.0/(KB * thermo.Tv() ); // [cm]
    const double fac3 = 2.0*HP*C0*C0*10000.0; // [W cm^2 sr-1]
    double f;
    #pragma omp parallel for private(sig)
    for (int j=0; j < mp_grid->size(); j++) {
        sig = (*mp_grid)[j]; // [cm-1]
        // Compute kappa
        f = std::exp(0.5*fac2*sig);
        p_spectra[j*2+1] = p_spectra[j*2]*(f-1.0/f)/(fac3*sig*sig)*f;
        //p_spectra[j*2+1] = p_spectra[j*2] *(std::exp(std::min(fac2*sig, MAXLOG))-1.0) /(fac3*sig*sig);
        // Warning: huge value of fac2*sig (high wavenumer, low temperature) => infinite exponential factor 
        //if (isnan(p_spectra[j*2+1]) || isinf(p_spectra[j*2+1])) p_spectra[j*2+1] =0.0;
        // Compute eta
        p_spectra[j*2] *= sig;
    }
    m_spectra_time += clock() - t1;

    // Write computed spectra
    writeLblSpectra(index, p_spectra);

}

void SnbAtomicSystem::lblToSnb(double* const p_lbl, double* const p_snb)
{

    int j;
    double delta_nu = 1000.0;
    double sigj, sigjp1, sig_tmp1, sig_tmp2, lbl_tmp1, lbl_tmp2, fac;

    // Contribution of lbl spectral intervals included in a spectral band
    #pragma omp parallel for private(j, sigj, sigjp1)
    for (int b=0; b < m_nbands; b++) {
       for (j = mp_bounds[b]; j < mp_bounds[b+1]-1; j++) {
          sigj = (*mp_grid)[j];
          sigjp1 = (*mp_grid)[j+1];
          p_snb[b]
             += 0.5*(p_lbl[j+1]+p_lbl[j])*(sigjp1-sigj);
       }
       p_snb[b] /= delta_nu;
    }

    // Contribution of lbl spectral intervals that cross two or more spectral bands
    int b=0;
    while(b < m_nbands-1) {
      j = mp_bounds[b+1]-1;
      sigj = (*mp_grid)[j];
      sigjp1 = (*mp_grid)[j+1];

      sig_tmp1 = (m_band1+b+1)*1000.0;
      lbl_tmp1 = (p_lbl[j+1]-p_lbl[j])/(sigjp1-sigj)
                    *(sig_tmp1-sigj) + p_lbl[j];

      p_snb[b]
          += 0.5*(lbl_tmp1+p_lbl[j])*(sig_tmp1-sigj)/delta_nu;

      // Deal with bands that might be contained within [sigj-sigjp1]
      while(sigjp1 > sig_tmp1+1000.0 ) {
         b++;
         sig_tmp2 = sig_tmp1+1000.0;
         lbl_tmp2 = (p_lbl[j+1]-p_lbl[j])/(sigjp1-sigj)
                       *(sig_tmp2-sigj) + p_lbl[j];

         p_snb[b] = 0.5*(lbl_tmp2+lbl_tmp1);

         sig_tmp1 = sig_tmp2;
         lbl_tmp1 = lbl_tmp2;
      }

      b++;
      p_snb[b]
          += 0.5*(p_lbl[j+1]+lbl_tmp1)*(sigjp1-sig_tmp1)/delta_nu;

    }

    // Rescale first and last band integration
    fac = (mp_grid->min() - m_band1*1000.0)/delta_nu;
    if (fac != 0.0) p_snb[0] *= 1.0/(1.0-fac);
    fac = ((m_bandn+1)*1000.0-mp_grid->max())/delta_nu;
    if (fac != 0.0) p_snb[m_nbands-1] *= 1.0/(1.0-fac);

}
/*
void SnbAtomicSystem::lblToSnb(double* const p_lbl, double* const p_snb)
{

  int jsnb=0;
  double delta_nu=1000.0;
  double grid_temp, lbl_temp;
  double sigj, sigjp1;

  for (int j=0; j<mp_grid->size()-1; j++) {
    sigj = (*mp_grid)[j];
    sigjp1 = (*mp_grid)[j+1];

    if (sigjp1 > (jsnb+2)*1000.0) { // The lbl grid crosses a band bound

      // Linear interpolation at the band bound
      grid_temp=(jsnb+2)*1000.0;
      lbl_temp=(p_lbl[j+1]-p_lbl[j])/(sigjp1-sigj)
                   *(grid_temp-sigj) + p_lbl[j];

      // Contribution to the left band
      p_snb[jsnb]
          += 0.5*(lbl_temp+p_lbl[j])*(grid_temp-sigj)/delta_nu;
      jsnb = jsnb +1;

      // Warning: the lbl grid does not cover the whole last snb band
      if (jsnb == m_nbands-1) { 
        delta_nu += (*mp_grid)[mp_grid->size()-1]-200000.0; }

      // Contribution to the right band
      p_snb[jsnb]
          += 0.5*(p_lbl[j+1]+lbl_temp)*(sigjp1-grid_temp)/delta_nu;

    } else { // sigj and sigjp1 within the interval [(jsnb+1)*1000 ; (jsnb+2)*1000] 
      p_snb[jsnb]
          += 0.5*(p_lbl[j+1]+p_lbl[j])*(sigjp1-sigj)/delta_nu;
    }
  }

}
*/
void SnbAtomicSystem::determineBandRange()
{

    double sigmin = mp_grid->min();
    double sigmax = mp_grid->max();

    if (sigmin < 1.0e3 || sigmax > 2.0e5)
       cout << "Error in the minimum and maximum LBL wavenumber" << endl;

    m_band1 = sigmin / 1000;
    m_bandn = sigmax / 1000;

    if (m_bandn * 1000 == sigmax)
       m_bandn--;

    m_nbands = m_bandn - m_band1 + 1;

}

void SnbAtomicSystem::lblToSnbInit()
{

    mp_bounds = new int [m_nbands+1];
    mp_bounds[0] = 0;
    for (int b = 1; b < m_nbands; ++b)
        mp_bounds[b] = mp_grid->index((b+1)*1000.0);
    mp_bounds[m_nbands] = mp_grid->size();

//    for (int b = 0; b < m_nbands; ++b)
//        cout << "(" << b << " - " << b+1 << ") " << (*mp_grid)[mp_bounds[b]] << ", " << (*mp_grid)[mp_bounds[b+1]-1] << endl;
//
//    exit(1);
//
//
//    double sigjp1, sigbp1;
//    mp_start = new int [m_nbands];
//    mp_end   = new int [m_nbands];
//    int b=0;
//    sigbp1 = (m_band1+1)*1000.0;
//    mp_start[b] = 0;
//
//    for (int j=0; j<mp_grid->size()-1; j++) {
//       sigjp1 = (*mp_grid)[j+1];
//       if (sigjp1 > sigbp1) {
//          mp_end[b] = j;
//          b++;
//          sigbp1 = (m_band1+b+1)*1000.0;
//          while (sigjp1 > sigbp1) {
//             mp_start[b] = j;
//             mp_end[b] = j;
//             b++;
//             sigbp1;
//          }
//          mp_start[b] = j+1;
//       }
//    }
//    mp_end[b] = mp_grid->size()-1;
//
//    for (int b = 0; b < m_nbands; ++b)
//        cout << "(" << b << " - " << b+1 << ") " << (*mp_grid)[mp_start[b]] << ", " << (*mp_grid)[mp_end[b]] << endl;

}
