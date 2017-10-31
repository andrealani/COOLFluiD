#if __linux__
#include "omp.h"
#endif
#include "RteSolution.h"
#include "Constants.h"

using namespace std;

RteSolution::RteSolution(const FieldData& field,
                         SnbAtomicSystem* atoms,
                         vector<SnbDiatomicSystem>& diatomics,
                         vector<SnbContinuumSystem>& continua,
                         SnbCO2System* co2)
    : m_field(field), m_atoms(atoms), m_diatomics(diatomics), m_continua(continua),
      m_co2(co2), m_spectra_time(0), m_spec_flux_time(0)
{
    determineBandRange();

    // Setup local radiative properties for diatomics systems and continua
    clock_t t1 = clock();
    for (int i = 0; i < m_diatomics.size(); ++i)
        m_diatomics[i].setupLocalParameters(m_field);
    for (int i = 0; i < m_continua.size(); ++i)
        m_continua[i].setupLocalParameters(m_field);
    if (m_co2)
        m_co2->setupLocalParameters(m_field);
    if (m_atoms) {
        m_atoms->setupLocalParameters(m_field);
    }
    m_spectra_time = clock() - t1;
}

RteSolution::~RteSolution()
{

}

double RteSolution::eqInt(double nu, double T)
// nu (cm-1), t (K), eqInt (W.cm-2.sr-1.cm) 
{
    double const c1 = 2.0*HP*C0*C0*10000.0;
    double const c2 = HP*C0*100.0/KB;
    double const max_power = std::log(std::numeric_limits<double>::max());

    if (T <= 0.0)
        return 0.0;
    else {
        // Split exponent to avoid overflow with small T and large nu
        double e = 0.5*c2*nu/T;
        if (e < max_power) {
            e = std::exp(e);
            return c1/e*nu*nu*nu/(e-1.0/e);
        } else {
            // In the limit of very large exponents, switch to Wien's law
            e = std::exp(-e);
            return c1*e*nu*nu*nu*e;
        }
    }
}

void RteSolution::setupBC(double* const p_wallInt, double twall, bool by_system)
{
    const int n_systems = (by_system ? 3+m_diatomics.size()+m_continua.size() : 1);
    for (int b=0; b< m_nbands; b++ ) {
        p_wallInt[b*n_systems] = eqInt((m_bmin+b)*1000+500, twall);
    }

}

void RteSolution::computeChemicalSourceTerms(double* const p_chemSource,
                                             const double* const p_specUnu)
{

    GlobalThermoData& thermo = GlobalThermoData::getInstance();
    int ncells = m_field.nCells();
    int index_elec = thermo.speciesIndex("e-");
    int bmin, bmax, len, index_species, index_products;
    double eta, kappa, unu, sigma, production_rate;
    string type;

    fill(p_chemSource, p_chemSource+ncells*thermo.nSpecies(), 0.0);

    for (int i=0; i < m_continua.size(); i++) {

        bmin = m_continua[i].lowBand()-m_bmin;
        bmax = m_continua[i].highBand()+1-m_bmin;
        len = m_continua[i].systemName().length();
        type = m_continua[i].systemName().substr(len-2,2);
        index_species = thermo.speciesIndex(m_continua[i].speciesName());
        index_products = thermo.speciesIndex(m_continua[i].productsName());

        if ( type == "bf" ) { // photoionization processes

            for (int icell =0; icell < ncells; icell++) {

                production_rate = 0.0;

                for (int b=bmin; b<bmax; ++b) {
                    eta = m_continua[i].getLocalParameter(icell, b-bmin, 1);
                    kappa = m_continua[i].getLocalParameter(icell, b-bmin, 0);
                    unu = (p_specUnu[icell*m_nbands+b] + p_specUnu[(icell+1)*m_nbands+b]) /2.0;
                    sigma = (b+m_bmin)*1000.0+500.0;
                    // production_rate in mole per second per cubic centimeter
                    production_rate += (unu * kappa - FOURPI * eta ) *1000.0 / (sigma*HP*C0*100.0);
                }

                p_chemSource[index_species *ncells+icell] -= thermo[index_species ].mw * production_rate / NA; 
                p_chemSource[index_elec    *ncells+icell] += thermo[index_elec    ].mw * production_rate / NA; 
                p_chemSource[index_products*ncells+icell] += thermo[index_products].mw * production_rate / NA; 
            }
        }

        if ( type == "SR") { // oxygen photodissociation

            for (int icell =0; icell < ncells; icell++) {

                production_rate = 0.0;

                for (int b=bmin; b<bmax; ++b) {
                    eta = m_continua[i].getLocalParameter(icell, b-bmin, 1);
                    kappa = m_continua[i].getLocalParameter(icell, b-bmin, 0);
                    unu = (p_specUnu[icell*m_nbands+b] + p_specUnu[(icell+1)*m_nbands+b]) /2.0;
                    sigma = (b+m_bmin)*1000.0+500.0;
                    // production_rate in mole per second per cubic centimeter
                    production_rate += (unu * kappa - FOURPI * eta ) *1000.0 / (sigma*HP*C0*100.0);
                }

                p_chemSource[index_species *ncells+icell] -= thermo[index_species ].mw * production_rate / NA; 
                p_chemSource[index_products*ncells+icell] += 2.0 * thermo[index_products].mw * production_rate / NA; 
            }

        }

    }

}

void RteSolution::computeOmegaSourceTerm(double* const p_omegaSource,
                                         const double* const p_specUnu)
{

    GlobalThermoData& thermo = GlobalThermoData::getInstance();
    int ncells = m_field.nCells();
    int bmin, bmax, len;
    double eta, kappa, unu, sigma, production_rate;
    string type;

    fill(p_omegaSource, p_omegaSource+m_field.nCells(), 0.0);

    for (int i=0; i < m_continua.size(); i++) {

        bmin = m_continua[i].lowBand()-m_bmin;
        bmax = m_continua[i].highBand()+1-m_bmin;
        len = m_continua[i].systemName().length();
        type = m_continua[i].systemName().substr(len-2,2);

        if ( type == "bf" ) { // photoionization processes

            double Eion = m_continua[i].chi()->ionizationEnergy(); // J/mol
            for (int icell =0; icell < ncells; icell++) {

                production_rate = 0.0;
                for (int b=bmin; b<bmax; ++b) {
                    eta = m_continua[i].getLocalParameter( icell, b-bmin, 1);
                    kappa = m_continua[i].getLocalParameter( icell, b-bmin, 0);
                    unu = (p_specUnu[icell*m_nbands+b] + p_specUnu[(icell+1)*m_nbands+b]) /2.0;
                    sigma = (b+m_bmin)*1000.0+500.0;
                    // production_rate in mol per second per cubic centimeter
                    production_rate += (unu * kappa - FOURPI * eta ) *1000.0 / (sigma*HP*C0*100.0);
                }

                // Omega source term (W/cm3)
                p_omegaSource[icell] -= production_rate /NA * Eion ;
            }
        }

    }

}

// Integrates the spectral flux over wavenumber and compute the divergence
void RteSolution::computeEnergySourceTerm(double* const p_energySource, const double* const p_specFlux)
{

    double* p_flux = new double [m_field.nPoints()];
    fill(p_flux , p_flux + m_field.nPoints(), 0.0);

    for (int i=0; i<m_field.nPoints(); i++) {
        for (int b=0; b < m_nbands ; b++ ) {
            p_flux[i] += p_specFlux[i*m_nbands + b] * 1000.0;
        }
    }

    string filename = "flux.dat";
    ofstream results;
    results.open(filename.c_str());

    for (int i = 0; i < m_field.nPoints(); ++i) {
        results << setw(10) << m_field.loc(i)
                << setw(14) << p_flux[i] << endl;
    }
    results.close();


    double xi, xip1;
    switch (m_field.geom()) {
        case CARTESIAN:
            xip1 = m_field.loc(0);
            for (int icell=0; icell<m_field.nCells(); icell++) {
                xi   = xip1;
                xip1 = m_field.loc(icell+1);
                p_energySource[icell] = (p_flux[icell] - p_flux[icell+1]) / (xip1 - xi);
            }
            break;
        case SPHERICAL:
            xip1 = m_field.loc(0);
            for (int icell=0; icell<m_field.nCells(); icell++) {
                xi   = xip1;
                xip1 = m_field.loc(icell+1);
                p_energySource[icell] = 3.0 * (p_flux[icell]*xi*xi - p_flux[icell+1]*xip1*xip1)
                                      / (pow(xip1,3) - pow(xi,3));
            }
            break;
    }

    delete[] p_flux;

}

void RteSolution::computeFluxField(double* const p_specFlux, double* const p_specUnu)
{
    clock_t t1 = clock();
    // Define angular discretization
    const int nmu = 20; // nmu should be even otherwise mu=0 is chosen
    double dmu [nmu], mu[nmu];
    // Regular mu grid
    for (int imu =0; imu<nmu; imu++) {
        dmu[imu] = (double)2.0 /(double)nmu ;
        mu[imu] = 1.0 - ((double)imu+0.5) * dmu[imu] ;
    }

    // Regular mu*mu grid
//    for (int imu =0; imu<nmu/2; imu++) {
//        dmu[imu] = ((imu+1)*(imu+1)-imu*imu) * pow(2.0 /(double)nmu, 2.0);
//        mu[imu] = 1.0 - pow(((double)imu+0.5)*2.0/(double) nmu, 2.0);
//        dmu[nmu-1-imu] = dmu[imu];
//        mu[nmu-1-imu] = - mu[imu];
//    }

//    for (int imu =0; imu<nmu; imu++) {
//        cout << "DIRECTION " << mu[imu] << " " << dmu[imu] << endl;
//    }

    // Allocation of variables
    double* p_specInt = new double [m_nbands];
    double* p_specWallInt = new double [2*m_nbands];
    double* p_specWallTemp = new double [2*m_nbands];
    int iwall, pos, wpos, it=0;
    bool convergence ;
    fill(p_specFlux, p_specFlux+m_field.nPoints()*m_nbands, 0.0);
    fill(p_specUnu , p_specUnu +m_field.nPoints()*m_nbands, 0.0);

    // Initialization of wall intensities
    setupBC(p_specWallInt, m_field.twall(0));
    setupBC(p_specWallInt+m_nbands, m_field.twall(1));
    if (m_field.epsilon(0) == 1.0 & m_field.epsilon(1) == 1.0 ) {
        convergence = true; 
    } else {
        convergence = false;
    }

    // Computation of wall leaving intensities until convergence if the two walls are not black
    // p_specWallInt: leaving intensities
    // p_specWallTemp: leaving intensities at the previous iteration
    while(!convergence) {
        it++;
        cout << "LEAVING INTENSITY CALCULATION - ITERATION " << it << endl;

        copy(p_specWallInt, p_specWallInt+2*m_nbands, p_specWallTemp );
        fill(p_specWallInt, p_specWallInt+2*m_nbands, 0.0 ); 

        for (int imu=0; imu<nmu; imu++) {
            if (mu[imu] > 0.) { pos = m_field.nPoints()-1; wpos = 1; }
            if (mu[imu] < 0.) { pos = 0; wpos  = 0; }

            LineOfSight los(m_field, pos, mu[imu]);
            iwall = los.wallNum();
            copy(p_specWallTemp+iwall*m_nbands, p_specWallTemp+(iwall+1)*m_nbands, p_specInt );
            computePath(los, p_specInt);

            // Computation of reflected intensity
            for (int b=0; b<m_nbands; b++) {
                p_specWallInt[wpos*m_nbands + b] += 2. * std::abs(mu[imu]) * dmu[imu] * p_specInt[b] * (1-m_field.epsilon(wpos)) ;
            }
        }

        convergence = true;
        // Add emitted intensity and check convergence
        for (int b=0; b< 2*m_nbands; b++) {
            if (b >= m_nbands ) { wpos = 1; } else { wpos = 0; }
            p_specWallInt[b] += m_field.epsilon(wpos) * eqInt((m_bmin+b-wpos*m_nbands)*1000.+500.,m_field.twall(wpos)) ;
            if (abs(p_specWallInt[b]-p_specWallTemp[b]) > 1.e-6 ) {
                convergence = false;
            }    
        }
    }

    // Computation of the flux

    for (int imu=0; imu<nmu; imu++) {
        cout << "DIRECTION " << imu << endl;
        for (int pos=0; pos < m_field.nPoints(); pos++) {
            // Step 1: draw the path from point pos in the direction mu
            LineOfSight los(m_field, pos, mu[imu]);
            iwall = los.wallNum();
            // Step 2: get the leaving intensity from the wall crossed by the path
            copy(p_specWallInt+iwall*m_nbands, p_specWallInt+(iwall+1)*m_nbands, p_specInt );
            // Step 3: compute intensity at point pos in the direction mu
            computePath(los, p_specInt);
            // Step 4: add the corresponding contribution to the flux
            #pragma omp parallel for
            for (int b=0; b<m_nbands; b++) {
                p_specFlux[pos*m_nbands + b] += TWOPI * mu[imu] * dmu[imu] * p_specInt[b] ;
                p_specUnu[pos*m_nbands + b]  += TWOPI * dmu[imu] * p_specInt[b] ;
            }
        }
    }
    m_spec_flux_time = clock() - t1;

    writeFieldResults("flux", p_specFlux);

    delete [] p_specWallInt;
    delete [] p_specWallTemp;
    delete [] p_specInt;

}

void RteSolution::computeIntensity(bool by_system)
{

    LineOfSight los(m_field, m_field.nPoints()-1,1.0);
    //LineOfSight los(m_field, 0,-0.70711);
    int iwall = los.wallNum();

    // Step2: setup BC
    const int size =
        m_nbands*(by_system ? 3+m_diatomics.size()+m_continua.size() : 1);
    double* p_specInt = new double [size];
    std::fill(p_specInt, p_specInt+size, 0.0);
    setupBC(p_specInt, m_field.twall(iwall), by_system);

    computePath(los, p_specInt, by_system);

    writeResults("intensity", p_specInt, by_system);

    double I = 0.0;
    for (int k = 0; k < m_nbands; ++k)
        I += p_specInt[k] * 1000.0;
    cout << "total intensity = " << I << endl;

    delete [] p_specInt;

}

void RteSolution::intensityCheck()
{
    double* p_specInt = new double [m_nbands];

    // METHOD 1
    // Initialize the boundary values
    std::fill(p_specInt, p_specInt+m_nbands, 0.0);
    setupBC(p_specInt, m_field.twall(0));

    // Loop over each point of the flow field
    for (int i = 1; i < m_field.nPoints(); ++i) {
        // Compute intensity at point
        LineOfSight los(m_field, i-1, i);
        computePath(los, p_specInt);

        // Compute the total intensity
        double I = 0.0;
        for (int k = 0; k < m_nbands; ++k)
            I += p_specInt[k] * 1000.0;

        cout << setw(15) << I;
    }
    cout << endl;

    // METHOD 2
    // Loop over each point of the flow field
    for (int i = 1; i < m_field.nPoints(); ++i) {
        // Initialize the boundary values
        std::fill(p_specInt, p_specInt+m_nbands, 0.0);
        setupBC(p_specInt, m_field.twall(0));

        // Compute intensity at point
        LineOfSight los(m_field, i, 1.0);
        computePath(los, p_specInt);

        // Compute the total intensity
        double I = 0.0;
        for (int k = 0; k < m_nbands; ++k)
            I += p_specInt[k] * 1000.0;

        cout << setw(15) << I;
    }
    cout << endl;
}

void RteSolution::computeTau()
{
    LineOfSight los(m_field, 0, -1.0);

    double* p_tau = new double [(los.nPath()+1)*m_nbands];
    fill(p_tau, p_tau+(los.nPath()+1)*m_nbands, 1.0);

    // Computation of the total transmissivities along the path
    if (m_atoms) {
       static_cast<RadiativeSystem<SnbAtomicSystem>&>(m_atoms[0]).tau(
           los, p_tau, m_atoms->lowBand()-m_bmin, m_nbands, TimesEq());
    }
    for (int i = 0; i < m_diatomics.size(); ++i) {
       static_cast<RadiativeSystem<SnbDiatomicSystem>&>(m_diatomics[i]).tau(
           los, p_tau, m_diatomics[i].lowBand()-m_bmin, m_nbands, TimesEq());
    }
    for (int i = 0; i < m_continua.size(); ++i) {
       static_cast<RadiativeSystem<SnbContinuumSystem>&>(m_continua[i]).tau(
           los, p_tau, m_continua[i].lowBand()-m_bmin, m_nbands, TimesEq());
    }
    if (m_co2) {
       static_cast<RadiativeSystem<SnbCO2System>&>(m_co2[0]).tau(
           los, p_tau, m_co2->lowBand()-m_bmin, m_nbands, TimesEq());
    }

    ofstream f("tau.dat");
    for (int b = 0; b < m_nbands; ++b)
        f << setw(10) << (b+m_bmin)*1000+500
          << setw(15) << p_tau[los.nPath()*m_nbands + b] << endl;
    f.close();
}

void RteSolution::computePath(const LineOfSight& los, double* const p_specInt,
    bool by_system)
{

    double etadx, etakappa, total_tau, delta_tau, len, prod;
    int icell, bmin, bmax;
    double* p_tauSyst = new double [(los.nPath()+1)*m_nbands];
    double* p_tau = new double [(los.nPath()+1)*m_nbands];
    fill(p_tau, p_tau+(los.nPath()+1)*m_nbands, 1.0);

    // Total number of "systems" to consider
    const int n_systems = 3 + m_diatomics.size() + m_continua.size();

    // Computation of the total transmissivities along the path
    if (m_atoms) {
       static_cast<RadiativeSystem<SnbAtomicSystem>&>(m_atoms[0]).tau(
           los, p_tau, m_atoms->lowBand()-m_bmin, m_nbands, TimesEq());
       // Saving of atomic transmissivities (expensive computation)
       copy(p_tau, p_tau + (los.nPath()+1)*m_nbands, p_tauSyst);
    }
    for (int i = 0; i < m_diatomics.size(); ++i) {
       static_cast<RadiativeSystem<SnbDiatomicSystem>&>(m_diatomics[i]).tau(
           los, p_tau, m_diatomics[i].lowBand()-m_bmin, m_nbands, TimesEq());
    }
    for (int i = 0; i < m_continua.size(); ++i) {
       static_cast<RadiativeSystem<SnbContinuumSystem>&>(m_continua[i]).tau(
           los, p_tau, m_continua[i].lowBand()-m_bmin, m_nbands, TimesEq());
    }
    if (m_co2) {
       static_cast<RadiativeSystem<SnbCO2System>&>(m_co2[0]).tau(
           los, p_tau, m_co2->lowBand()-m_bmin, m_nbands, TimesEq());
    }


    int index;
    // Contributions from the wall
    #pragma omp parallel for
    for (int b=0; b<m_nbands; ++b) {
        index = (by_system ? b*n_systems + 0 : b);
        p_specInt[index] *= p_tau[los.nPath()*m_nbands+b];
    }

    // Contributions from atoms
    if (m_atoms) {
        bmin = m_atoms->lowBand()-m_bmin;
        bmax = m_atoms->highBand()+1-m_bmin;
        // (p_tauSyst precomputed)
        for (int k=los.nPath()-1; k>=0; k--) { // Loop over the path
            icell = los.cellNum(k);
            len = los.length(k);
            for (int b=bmin; b<bmax; ++b) {
                index = (by_system ? b*n_systems + 1 : b);
                etakappa = m_atoms->getLocalParameter( icell, b-bmin, 0);
                total_tau = sqrt(p_tau[(k+1)*m_nbands+b] * p_tau[k*m_nbands+b]);

                //cout << b << ", " << i << ", " << k << endl;
                //cout << p_tauSyst[k*m_nbands+b] << ", " << p_tauSyst[(k+1)*m_nbands+b] << endl;

                delta_tau = p_tauSyst[k*m_nbands+b] * p_tauSyst[(k+1)*m_nbands+b];
                if (delta_tau > 0.0) {
                	delta_tau = (p_tauSyst[k*m_nbands+b] - p_tauSyst[(k+1)*m_nbands+b])
                       / sqrt(delta_tau);
                	if (!(delta_tau != delta_tau || delta_tau > 1.0e100))
                		p_specInt[index] += etakappa * delta_tau * total_tau;
                }
                if (p_specInt[index]> 1.0e100) {
                    cout << "WARN " << b << " " 
                                    << icell << " "
                                    << etakappa << " "
                                    << delta_tau << " "
                                    << total_tau << " " << endl;
                }
            }
        } 
    }

    // Contributions from molecular systems
    for (int i = 0; i < m_diatomics.size(); ++i) {
        bmin = m_diatomics[i].lowBand()-m_bmin;
        bmax = m_diatomics[i].highBand()+1-m_bmin;
        if (m_diatomics[i].nParams() < 3) { // Thin Systems
            for (int k=los.nPath()-1; k>=0; k--) { // Loop over the path
                icell = los.cellNum(k);
                len = los.length(k);
                #pragma omp parallel for private(etadx, total_tau, index)
                for (int b=bmin; b<bmax; ++b) {
                    index = (by_system ? b*n_systems+2+i : b);
                    etadx = m_diatomics[i].getLocalParameter( icell, b-bmin, 1) * len;
                    total_tau = sqrt(p_tau[(k+1)*m_nbands+b] * p_tau[k*m_nbands+b]);
                    p_specInt[index] += etadx * total_tau;
                }
            }
        } else { // Thick Systems
            fill(p_tauSyst, p_tauSyst+(los.nPath()+1)*m_nbands, 1.0);
            static_cast<RadiativeSystem<SnbDiatomicSystem>&>(m_diatomics[i]).tau(
                los, p_tauSyst, bmin, m_nbands, Eq());
            for (int k=los.nPath()-1; k>=0; k--) { // Loop over the path
                icell = los.cellNum(k);
                len = los.length(k);
                #pragma omp parallel for private(etakappa, total_tau, delta_tau, index)
                for (int b=bmin; b<bmax; ++b) {
                    index = (by_system ? b*n_systems+2+i : b);
                    etakappa = m_diatomics[i].getLocalParameter( icell, b-bmin, 1);
                    total_tau = sqrt(p_tau[(k+1)*m_nbands+b] * p_tau[k*m_nbands+b]);
                    prod = p_tauSyst[k*m_nbands+b]*p_tauSyst[(k+1)*m_nbands+b];

                    if (prod > 0.0) {
                        delta_tau = (p_tauSyst[k*m_nbands+b] - p_tauSyst[(k+1)*m_nbands+b])
                               / sqrt(prod);
                        p_specInt[index] += etakappa * delta_tau * total_tau;
                    }
                }
            } 
        }
    }

    // Contributions from continua
    for (int i = 0; i < m_continua.size(); ++i) {
        bmin = m_continua[i].lowBand()-m_bmin;
        bmax = m_continua[i].highBand()+1-m_bmin;
        for (int k=los.nPath()-1; k>=0; k--) { // Loop over the path
            icell = los.cellNum(k);
            len = los.length(k);
            #pragma omp parallel for private(etadx, total_tau, index)
            for (int b=bmin; b<bmax; ++b) {
                index = (by_system ? b*n_systems+2+m_diatomics.size()+i : b);
                etadx = m_continua[i].getLocalParameter( icell, b-bmin, 1) * len;
                total_tau = sqrt(p_tau[(k+1)*m_nbands+b] * p_tau[k*m_nbands+b]);
                p_specInt[index] += etadx * total_tau;
            }
        }
    }

    // Contribution from CO2
    if (m_co2) {
        bmin = m_co2->lowBand()-m_bmin;
        bmax = m_co2->highBand()+1-m_bmin;

        fill(p_tauSyst, p_tauSyst+(los.nPath()+1)*m_nbands, 1.0);
        static_cast<RadiativeSystem<SnbCO2System>&>(m_co2[0]).tau(
            los, p_tauSyst, bmin, m_nbands, Eq());
        for (int k=los.nPath()-1; k>=0; k--) { // Loop over the path
            icell = los.cellNum(k);
            len = los.length(k);
            #pragma omp parallel for private(etakappa, total_tau, delta_tau, index)
            for (int b=bmin; b<bmax; ++b) {
                index = (by_system ? (b+1)*n_systems-1 : b);
                etakappa = m_co2->getLocalParameter(icell, b-bmin, 1);
                total_tau = sqrt(p_tau[(k+1)*m_nbands+b] * p_tau[k*m_nbands+b]);
                delta_tau = (p_tauSyst[k*m_nbands+b] - p_tauSyst[(k+1)*m_nbands+b])
                       / sqrt(p_tauSyst[k*m_nbands+b] * p_tauSyst[(k+1)*m_nbands+b]);
                if (!(std::isnan(delta_tau) || std::isinf(delta_tau)))
                    p_specInt[index] += etakappa * delta_tau * total_tau;
            }
        } 
    }

    delete [] p_tau;
    delete [] p_tauSyst;

}

void RteSolution::writeFieldResults(const string& field_name, const double* const p_field)
{

    string filename;
    ofstream results;

    filename = "spectral_" + field_name + "_field.dat";
    results.open(filename.c_str());

    for (int b = 0; b < m_nbands; ++b) {
        results << setw(10) << (m_bmin+b)*1000 + 500;
        for (int j = 0; j < m_field.nPoints(); ++j)
            results << setw(14) << p_field[j*m_nbands+b];
        results << endl;
    }
    results.close();


    double* p_cumul = new double [m_field.nPoints()];
    fill(p_cumul, p_cumul+m_field.nPoints(), 0.0);
    filename = "cumulated_" + field_name + "_field.dat";
    results.open(filename.c_str());

    for (int b = 0; b < m_nbands; ++b) {
        results << setw(10) << (m_bmin+b)*1000 + 500;
        for (int j = 0; j < m_field.nPoints(); ++j) {
            p_cumul[j] += p_field[j*m_nbands+b] * 1000.0;
            results << setw(14) << p_cumul[j];
        }
        results << endl;
    }
    results.close();

}

void RteSolution::writeResults(const string& field_name, const double* const p_field, bool by_system)
{

    string filename;
    ofstream results;

    filename = "spectral_" + field_name + ".dat";
    results.open(filename.c_str());

    // Write a header
    if (by_system) {
        results << setw(10) << "sigma";
        results << setw(14) << "Wall";
        results << setw(14) << "Atoms";

        for (int i = 0; i < m_diatomics.size(); ++i)
            results << setw(14) << m_diatomics[i].speciesName()+"/"+m_diatomics[i].systemName();

        for (int i = 0; i < m_continua.size(); ++i)
            results << setw(14) << m_continua[i].speciesName()+"/"+m_continua[i].systemName();

        results << setw(14) << "CO2";
        results << endl;
    } else {
        results << setw(10) << "sigma";
        results << setw(14) << "I";
        results << endl;
    }

    // Write the spectral quantities
    int n_contributions = (by_system ? 3+m_diatomics.size()+m_continua.size() : 1);
    for (int b = 0; b < m_nbands; ++b) {
        results << setw(10) << (m_bmin+b)*1000 + 500;
        for (int i = 0; i < n_contributions; ++i)
            results << setw(14) << p_field[b*n_contributions+i];
        results << endl;
    }
    results.close();


    double p_cumul = 0.;
    filename = "cumulated_" + field_name + ".dat";
    results.open(filename.c_str());

    for (int b = 0; b < m_nbands; ++b) {
        for (int i = 0; i < n_contributions; ++i)
            p_cumul += p_field[b*n_contributions+i] * 1000.0;
        results << setw(10) << (m_bmin+b)*1000 + 500
                << setw(14) << p_cumul << endl;
    }
    results.close();

}

void RteSolution::determineBandRange()
{

    if (m_diatomics.size() != 0) {
        m_bmin = m_diatomics[0].lowBand();
        m_bmax = m_diatomics[0].highBand();
    } else if (m_continua.size() !=0) {
        m_bmin = m_continua[0].lowBand();
        m_bmax = m_continua[0].highBand();
    } else if (m_co2) {
        m_bmin = m_co2->lowBand();
        m_bmax = m_co2->highBand();
    } else {
        m_bmin = m_atoms->lowBand();
        m_bmax = m_atoms->highBand();
    }

    for (int i = 0; i < m_diatomics.size(); ++i) {
        m_bmin = min(m_bmin, m_diatomics[i].lowBand());
        m_bmax = max(m_bmax, m_diatomics[i].highBand());
    }

    for (int i = 0; i < m_continua.size(); ++i) {
        m_bmin = min(m_bmin, m_continua[i].lowBand());
        m_bmax = max(m_bmax, m_continua[i].highBand());
    }

    if (m_co2) {
        m_bmin = min(m_bmin, m_co2->lowBand());
        m_bmax = max(m_bmax, m_co2->highBand());
    }

    if (m_atoms) {
        m_bmin = min(m_bmin, m_atoms->lowBand());
        m_bmax = max(m_bmax, m_atoms->highBand());
    }

    m_nbands = m_bmax-m_bmin+1;
}
