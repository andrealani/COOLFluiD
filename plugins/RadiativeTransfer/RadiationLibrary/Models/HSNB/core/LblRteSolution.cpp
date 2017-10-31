#if __linux__
#include "omp.h"
#endif
#include "LblRteSolution.h"
#include "Constants.h"

using namespace std;

LblRteSolution::LblRteSolution(const FieldData& field,
                         const LblSpectralGrid& lblgrid,
                         double* const spectrum)
    : m_field(field), m_grid(lblgrid), mp_spectrum(spectrum)
{

}

LblRteSolution::~LblRteSolution()
{

}

double LblRteSolution::eqInt(double nu,double t)
// nu (cm-1), t (K), eqInt (W.cm-2.sr-1.cm) 
{
    double const c1=2*HP*pow(C0*100,2.0);
    double const c2= HP*C0*100/KB;
    double res;

    if (t == 0.) { return res=0.;}
    else { return res=c1*pow(nu,3.)/(exp(c2*nu/t)-1.);}
}

void LblRteSolution::setupBC(double* const p_wallInt, double twall)
{

    for (int i=0; i< m_grid.size(); i++ ) {
        p_wallInt[i] = eqInt(m_grid[i], twall);
    }

}

void LblRteSolution::computeEnergySourceTerm(double* const p_energySource, const double* const p_specFlux)
{

    double* p_flux = new double [m_field.nPoints()];
    fill(p_flux , p_flux + m_field.nPoints(), 0.0);

    for (int i=0; i<m_field.nPoints(); i++) {
        for (int b=1; b < m_grid.size() ; b++ ) {
            p_flux[i] += 0.5*(p_specFlux[i*m_grid.size()+b]+p_specFlux[i*m_grid.size()+b-1]) * (m_grid[b]-m_grid[b-1]);
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
    xip1 = m_field.loc(0);
    for (int icell=0; icell<m_field.nCells(); icell++) {
       xi   = xip1;
       xip1 = m_field.loc(icell+1);
       p_energySource[icell] = (p_flux[icell] - p_flux[icell+1]) / (xip1 - xi);
    }

    delete[] p_flux;

}


void LblRteSolution::computeFluxField(double* const p_specFlux)
{

    int Np = m_field.nPoints(), Nb = m_grid.size();
    double* p_slabTau = new double [Np*Np];
    double opt_thick, delta_tau;

    // Allocation of variables
    double* p_specWallInt = new double [2*Nb];
    bool convergence ;
    fill(p_specFlux, p_specFlux+Np*Nb, 0.0);

    // Initialization of wall intensities
    setupBC(p_specWallInt, m_field.twall(0));
    setupBC(p_specWallInt+Nb, m_field.twall(1));
    if (m_field.epsilon(0) == 1.0 & m_field.epsilon(1) == 1.0 ) {
        convergence = true; 
    } else {
        convergence = false;
        cout << "FATAL ERROR, Non-black wall case not implemented yet !!!" << endl;
        return;
    }

    for (int b=0; b < Nb; b++) {

       // Setup the tau slab matrix
       fill(p_slabTau, p_slabTau + Np*Np, 0.0);

       for (int i=0; i<Np; i++) { // Loop over the points i 
          //case i=j
          opt_thick = 0.0;
          p_slabTau[i*Np+i] = E3(opt_thick);

          // Loop over the points j
          for (int j=i+1; j<Np; j++) {
            opt_thick += mp_spectrum[(j-1)*Nb*2+b*2+1] * (m_field.loc(j)-m_field.loc(j-1));
            p_slabTau[i*Np+j] = E3(opt_thick);
            // Reverse path
            p_slabTau[j*Np+i] = p_slabTau[i*Np+j] ;
          }

       }

       for (int i=0; i<Np; i++) { // Loop over the points i
  
          // Contributions from the walls
          p_specFlux [i*Nb + b] += TWOPI * p_specWallInt[   b] * p_slabTau[    0 *Np+i];
          p_specFlux [i*Nb + b] -= TWOPI * p_specWallInt[Nb+b] * p_slabTau[(Np-1)*Np+i];
  
          // Contributions from the plasma 
          for (int j=0; j<Np-1; j++) { // Loop over the cells j
             delta_tau = p_slabTau[(j+1)*Np+i] - p_slabTau[ j   *Np+i];
             p_specFlux [i*Nb + b] += TWOPI * mp_spectrum[j*Nb*2+b*2] * delta_tau;
          }

       }

    }

    //writeFieldResults("flux", p_specFlux);
    writeResults("flux", p_specFlux);

    delete [] p_specWallInt;
    delete [] p_slabTau;

}

/*
void LblRteSolution::computeFluxField(double* const p_specFlux, double* const p_specUnu)
{
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

    writeFieldResults("flux", p_specFlux);

    delete [] p_specWallInt;
    delete [] p_specWallTemp;
    delete [] p_specInt;

}
*/

void LblRteSolution::writeFieldResults(const string& field_name, const double* const p_field)
{

    string filename;
    ofstream results;

    filename = "spectral_" + field_name + "_field.dat";
    results.open(filename.c_str());

    for (int b = 0; b < m_grid.size(); ++b) {
        results << setw(10) << m_grid[b];
        for (int j = 0; j < m_field.nPoints(); ++j)
            results << setw(14) << p_field[j*m_grid.size()+b];
        results << endl;
    }
    results.close();


    double* p_cumul = new double [m_field.nPoints()];
    fill(p_cumul, p_cumul+m_field.nPoints(), 0.0);
    filename = "cumulated_" + field_name + "_field.dat";
    results.open(filename.c_str());

    for (int b = 1; b < m_grid.size(); ++b) {
        results << setw(10) << (m_grid[b]+m_grid[b-1])/2.0;
        for (int j = 0; j < m_field.nPoints(); ++j) {
            p_cumul[j] += 0.5*(p_field[j*m_grid.size()+b]+p_field[j*m_grid.size()+b-1]) * (m_grid[b]-m_grid[b-1]);
            results << setw(14) << p_cumul[j];
        }
        results << endl;
    }
    results.close();

}

void LblRteSolution::writeResults(const string& field_name, const double* const p_field)
{

    string filename;
    ofstream results;

    filename = "spectral_" + field_name + ".dat";
    results.open(filename.c_str());

    for (int b = 0; b < m_grid.size(); ++b) {
        results << setw(10) << m_grid[b]
                << setw(14) << p_field[b] << endl;
    }
    results.close();


    double p_cumul = 0.;
    filename = "cumulated_" + field_name + ".dat";
    results.open(filename.c_str());

    for (int b = 1; b < m_grid.size(); ++b) {
        p_cumul += 0.5*(p_field[b]+p_field[b-1]) * (m_grid[b]-m_grid[b-1]);
        results << setw(10) << (m_grid[b]+m_grid[b-1])/2.0
                << setw(14) << p_cumul << endl;
    }
    results.close();

}

double LblRteSolution::E3(double x)
{
    return (x <= 1.0e-80 ? 0.5 : 0.5*(exp(-x)*(1-x)+E1(x)*x*x));
}

double LblRteSolution::E2(double x)
{
        double res;

        if (x <=1.e-80)
        {
        res=1.;
        }
        else
        {
        res=exp(-x)-x*E1(x);
        }

        return res;
}

double LblRteSolution::E1(double x)
{
        double const xhi(7.038e2);
        double res, t, y;

        if (x <= 0.)
        {
        res=0.;
        }
        else if (x <= 1.)
        {
        t = 2.0e0*x - 1.0e0;
        y = ((((((((((((
        -2.68991715725239147e-14)*t +7.06801512848037888e-13)*t
        -1.70965489077552146e-11)*t +3.81717274972110040e-10)*t
        -7.77461368291180797e-09)*t +1.43197088445269424e-07)*t
        -2.36082288702093691e-06)*t +3.44231259904916577e-05)*t
        -4.37905639072711938e-04)*t +4.79589265565698012e-03)*t
        -4.51020052155249319e-02)*t +3.93469340287366573e-01)*t
        -1.33373585783784498e-01;
        res = y - log(x);
        }
        else if (x <= 2.)
        {
        t = 2.0e0*x - 3.0e0;
        res = (((((((((((((((((((((
        -8.37205032055466094e-12)*t +2.56178397184819264e-11)*t
        -3.46318473395090883e-11)*t +1.13647587973449109e-10)*t
        -4.71527752610644265e-10)*t +1.49398063525941572e-09)*t
        -4.63129693426107937e-09)*t +1.48977580091910592e-08)*t
        -4.82569495402871907e-08)*t +1.56826141006595685e-07)*t
        -5.13180763112114830e-07)*t +1.69349493659850299e-06)*t
        -5.64487355301463818e-06)*t +1.90487447194908340e-05)*t
        -6.52605660793704808e-05)*t +2.27604942478038604e-04)*t
        -8.07756431103733780e-04)*t +2.88381958526388800e-03)*t
        -9.98576333997560656e-03)*t +3.09903000206149474e-02)*t
        -7.43767200494766144e-02)*t +1.00019582406632653e-01;
        }
        else if (x <= 4.)
        {
        t = x - 3.0e0;
        res = (((((((((((((((((((((
        -8.37205032050043424e-12)*t +2.56178397172768657e-11)*t
        -3.46318473142044651e-11)*t +1.13647587461332047e-10)*t
        -4.71527742750325071e-10)*t +1.49398045511718119e-09)*t
        -4.63129381995412142e-09)*t +1.48977072175917445e-08)*t
        -4.82561707556553537e-08)*t +1.56814957829348347e-07)*t
        -5.13031023477588644e-07)*t +1.69163480242326026e-06)*t
        -5.62356212102635108e-06)*t +1.88251709859233953e-05)*t
        -6.31322400361725865e-05)*t +2.09438056045469880e-04)*t
        -6.70998555174443855e-04)*t +1.99762928637710204e-03)*t
        -5.22456890280011749e-03)*t +1.10637929706361252e-02)*t
        -1.65956894559546525e-02)*t +1.30483810941970387e-02;
        }
        else if (x <= xhi)
        {
        t = 14.5e0/(x+3.25e0) - 1.0e0;
        y = ((((((((((((((((-1.181683914463340e-11
        *t-1.048399621792147e-11)*t-6.804206766157201e-11)
        *t-2.211936213417643e-10)*t-1.119940934924999e-9)
        *t-3.663650270142246e-9)*t-1.668705862724518e-8)
        *t-6.118194012663569e-8)*t-2.703104833900864e-7)
        *t-1.055648511722045e-6)*t-4.720903675164850e-6)
        *t-1.950763781198026e-5)*t-9.164504841773118e-5)
        *t-4.058921304226172e-4)*t-2.142130549996267e-3)
        *t-1.063748751165807e-2)*t-8.506991549845731e-2)*t +
         9.237553078077841e-1;
        res = exp(-x)*y/x;
        }
        else{
        res=0.;
        }

        return res;
}
