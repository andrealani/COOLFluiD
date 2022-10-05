#include "Framework/MethodStrategyProvider.hh"
#include "Framework/MeshData.hh"
#include "Framework/DataHandle.hh"
#include "Framework/GeometricEntity.hh"
#include "MHD/MHDTerm.hh"
#include "FiniteVolumeMHD/FiniteVolumeMHD.hh"
#include "FiniteVolumeMHD/MHDConsACASourceTerm.hh"
#include "FiniteVolume/CellCenterFVMData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::MHD;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<MHDConsACASourceTerm,
		       CellCenterFVMData,
		       ComputeSourceTerm<CellCenterFVMData>,
		       FiniteVolumeMHDModule>
mHDConsACASTFVMCCProvider("MHDConsACAST");

//////////////////////////////////////////////////////////////////////////////

MHDConsACASourceTerm::MHDConsACASourceTerm(const std::string& name) :
  ComputeSourceTermFVMCC(name),
  socket_gravity("gravity")
{
  addConfigOptionsTo(this);

  // _gravity = 0; 
  setParameter("gravity",&_gravity);
  setParameter("PevtsovHeating",&_PevtsovHeating);
  setParameter("PevtsovHeatingFactor",&_PevtsovHeatingFactor);

  //_Manchester = 0;
  setParameter("Manchester",&_Manchester);
  setParameter("ManchesterHeatingAmplitude",&_ManchesterHeatingAmplitude);
  setParameter("ManchesterSigma",&_ManchesterSigma);
  setParameter("divQ",&_divQ);
  setParameter("divQConductivity",&_divQConductivity);
  setParameter("divQalphaCollisionless",&_divQalphaCollisionless);
  setParameter("ViscosityAndResistivity",&_ViscosityAndResistivity);
  setParameter("Viscosity",&_Viscosity);
  setParameter("Resistivity",&_Resistivity);
  
  //_RadiativeLossTerm = 0;
  setParameter("RadiativeLossTerm",&_RadiativeLossTerm);

  
  //_projectionIDs = vector<CFuint>();
  //setParameter("ProjectionIDs",&_projectionIDs);

}

//////////////////////////////////////////////////////////////////////////////


void MHDConsACASourceTerm::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFint >("gravity","Switch on gravity.");
  options.addConfigOption< CFint >("PevtsovHeating","Switch on PevtsovHeating term.");
  options.addConfigOption< CFreal >("PevtsovHeatingFactor","Scale PevtsovHeating term.");
  options.addConfigOption< CFint >("Manchester","Switch on Manchester heating term.");
  options.addConfigOption< CFreal >("ManchesterHeatingAmplitude","Scale Manchester volumetric heating amplitude.");
  options.addConfigOption< CFreal >("ManchesterSigma","Sigma value of Manchester heating term.");
  options.addConfigOption< CFint >("divQ","Switch on heat conduction.");
  options.addConfigOption< CFreal >("divQConductivity","Conductivity of heat conduction.");
  options.addConfigOption< CFreal >("divQalphaCollisionless","alpha value of heat conductivity in collisionless regime.");
  options.addConfigOption< CFint >("ViscosityAndResistivity","Switch on viscosity and resisitivity.");
  options.addConfigOption< CFreal >("Viscosity","nu value.");
  options.addConfigOption< CFreal >("Resistivity","eta value.");
  options.addConfigOption< CFint >("RadiativeLossTerm","Switch on optically thin approximation for radiation losses.");
}

//////////////////////////////////////////////////////////////////////////////

MHDConsACASourceTerm::~MHDConsACASourceTerm()
{
}

//////////////////////////////////////////////////////////////////////////////

void MHDConsACASourceTerm::setup()
{
  ComputeSourceTermFVMCC::setup();

  socket_gravity.getDataHandle().resize(socket_volumes.getDataHandle().size()*4);
}

//////////////////////////////////////////////////////////////////////////////

void MHDConsACASourceTerm::computeSource(Framework::GeometricEntity *const element,
					 RealVector& source,
					 RealMatrix& jacobian)
{
  CFLog(DEBUG_MAX, "MHDConsACASourceTerm::computeSource() => START\n");
 
  DataHandle<CFreal> volumes = socket_volumes.getDataHandle();

  SafePtr<MHDTerm> model = PhysicalModelStack::getActive()->getImplementor()->
	                  getConvectiveTerm().d_castTo<MHDTerm>();

  // this source term is for MHD flows
  const vector<State*>* const states = element->getStates();
  const CFuint elementID = element->getID();

  // all elements in FVM should have only one state
  cf_assert(states->size() == 1);

  State *const currState = (*states)[0];
  
  const CFreal refSpeed = model->getRefSpeed();
  const CFreal refSpeedSq = refSpeed*refSpeed;

  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFuint dim = PhysicalModelStack::getActive()->getDim(); 
  for (CFuint i = 0; i < (nbEqs-1); ++i) {
    source[i] = 0.0;
  }

  const std::string correctionType = model->getCorrectionType();

  if (correctionType == "Mixed") {
    // mixed (hyperbolic and parabolic) correction
    const CFreal dissipCoeff = model->getDissipationCoefficient();
    const CFreal dissipCoeffSq = dissipCoeff*dissipCoeff;

    source[nbEqs-1] = -(refSpeedSq/dissipCoeffSq)*
      (*currState)[nbEqs-1]*
      volumes[elementID];
  }
  else {
    // hyperbolic correction
    source[nbEqs-1] = 0.0;
  }

  //std::cout << "SOURCE TERM" << endl;

  //std::cout << "ST activated" << endl;
  CFreal RSun = 6.9551e8; // m
  CFreal x = currState->getCoordinates()[XX]*RSun;
  CFreal y = currState->getCoordinates()[YY]*RSun;
  CFreal z = (dim==DIM_3D) ? currState->getCoordinates()[ZZ]*RSun : 0.;
  CFreal x_dimless = currState->getCoordinates()[XX];
  CFreal y_dimless = currState->getCoordinates()[YY];
  CFreal z_dimless = (dim==DIM_3D) ? currState->getCoordinates()[ZZ] : 0.;
  CFreal r2_dimless = x_dimless*x_dimless + y_dimless*y_dimless + z_dimless*z_dimless;
  CFreal r_dimless = std::sqrt(r2_dimless);

  CFreal r2 = x*x + y*y + z*z;
  CFreal r = std::sqrt(r2);
  CFreal rho_ref = 1.67e-13;
  CFreal density = (*currState)[0]*rho_ref;
  CFreal GMsun = 1.327474512e20; // SI value
  CFreal gsun = 274.00; // at the surface m per sec2
  CFreal fgrav = -density*gsun*pow(6.9551e8,2)/r2;

  // gravity term in dimensionless units:
  CFreal fgravity_dimless = -0.82977*(*currState)[0]/r2_dimless;  
  // the numerical factor is - G*MSun*mu_0*rho_0/(B_0**2*R_Sun)

  CFreal fx = fgrav*x/r;
  CFreal fy = fgrav*y/r;
  CFreal fz = fgrav*z/r;
  CFreal force_density_ref = std::pow(2.2e-4,2)/(1.25e-6);//*6.9e8);
  //std::cout << "_____________" << endl;
  //std::cout << "r = " << r << endl;
  //std::cout << "rho = " << density << endl;
  //std::cout << "fgrav = " << fgrav << endl;
  //std::cout << "_____________" << endl;
  //source[1] = fx/force_density_ref;
  //source[2] = fy/force_density_ref;
  //source[3] = fz/force_density_ref;
  // VERSION 2
  //source[1] = fgravity_dimless*x_dimless/r_dimless;
  //source[2] = fgravity_dimless*y_dimless/r_dimless;
  //source[3] = fgravity_dimless*z_dimless/r_dimless;
  CFreal Vr_dimless = x_dimless/r_dimless*(*currState)[1] 
                      + y_dimless/r_dimless*(*currState)[2] 
                      + z_dimless/r_dimless*(*currState)[3];
  //source[7] = fgravity_dimless*Vr_dimless;



 // B-FIELD WEIGHTED 3-D HEATING MODEL - Petsov et al. 2003
    // First attempt: Instead of summing /int /psi dV over all tetrahedrons, which has to be done prior to calling this function
    // we can simply use the constant factor provided by Downs et al. 2010 of 4e-8 J/s T
    // Then the heating term is simply
    const CFreal const_Downs2010 = 4.e-8/0.03851; // to adimensional value    // Originally 4.e-8
    CFreal Qh2 = const_Downs2010*std::sqrt((*currState)[4]*(*currState)[4] + (*currState)[5]*(*currState)[5] + (*currState)[6]*(*currState)[6]);


  //const CFreal r = stateCoordsSpherical[0];
  //cout << "r = " << r_dimless << endl;
  //const CFreal theta = stateCoordsSpherical[1];
  const CFreal theta = std::atan2(std::sqrt(x_dimless*x_dimless + y_dimless*y_dimless),z_dimless);
  //cout << "theta = " << theta << endl;

  const CFreal th1 = (MathTools::MathConsts::CFrealPi()/180.0)*17.5;
  const CFreal th2 = (MathTools::MathConsts::CFrealPi()/180.0)*61.5;

  // critical angle
  CFreal critTheta = 0.;
  if (r_dimless <= 7.0) {
	critTheta = asin(sqrt(sin(th1)*sin(th1)+cos(th1)*cos(th1)*(r-1.0)/8.0));
  }
  if ((r_dimless > 7.0) && (r_dimless < 30.0)) {
	critTheta = asin(sqrt(sin(th2)*sin(th2)+cos(th2)*cos(th2)*(r-7.0)/40.0));
  }

  // volumetric heating amplitude (J/(kg*s*K))
  const CFreal q0 = _ManchesterHeatingAmplitude;
  // just a check:
  //cout << "Manchester heating amplitude = " << q0 << endl;

  // reference length = RSun
  //const CFreal lRef = _varSet->getLRef();
  const CFreal lRef = 6.9e8;
  //cout << "lRef = " << lRef << endl;

  // target temperature and the heating scale height
  CFreal T0 = 0.;
  CFreal sigma0 = 0.;
  if (fabs(0.5*MathTools::MathConsts::CFrealPi()-theta) < critTheta) {
    T0 = 1.5e6;
    sigma0 = 4.5*lRef;
  }
  if (fabs(0.5*MathTools::MathConsts::CFrealPi()-theta) > critTheta) {
    T0 = 2.63e6;
    sigma0 = 4.5*lRef*(2.0-(sin(theta)*sin(theta)/(sin(critTheta)*sin(critTheta))));
  }
  
  //_varSet->computePhysicalData(*currState, _physicalData);

  // Boltzmann constant
  const CFreal k = 1.3806503e-23;

  // mass of proton
  const CFreal mp = 1.67262158e-27;

  // mass of electron
  const CFreal me = 9.10938188e-31;

  //const CFreal nRef = _varSet->getNRef();
  //cout << "nRef = " << nRef << endl;
  //const CFreal TRef = _varSet->getTRef(); 
  //cout << "TRef = " << TRef << endl;

  //const CFreal rhoRef = nRef*(mp+me);
  const CFreal rhoRef = 1.67e-13;
  //cout << "rhoRef = " << rhoRef << endl;
  //const CFreal vRef = sqrt(2.0*k*TRef/mp);
  const CFreal vRef = 480.248;
  //cout << "VRef = " << vRef << endl;

  const CFreal T = (*currState)[7]*0.03851*1.27*mp*vRef*vRef/((*currState)[0]*1.67e-13*k);
  const CFreal T2 = (*currState)[7]*0.03851*1.27*mp*vRef*vRef/(2.0*(*currState)[0]*1.67e-13*k);
  //cout << "T = " << T << endl;
  //cout << "T with factor 2 in denominator = " << T2 << endl;

  // dimensional volumetric heating term (J/(m^3*s))
  const CFreal Q = (*currState)[0]*1.67e-13*q0*(T0-T2)*exp(-(r_dimless-1.0)*(r_dimless-1.0)/(sigma0*sigma0));
  //cout << "Q = " << Q << endl;

  const CFreal nondimconst = lRef/(rhoRef*pow(vRef,3));
  //cout << "QRef = " << nondimconst << endl;



  // Radiative loss pleitner Dec 4
  // Q_rad = n_e^2*10^(-17.73)*T^(-2./3.) = (rho/(mu*mH))^2*10^(-17.73)*T^(-2./3.)
  const CFreal mu_cor = 1.27; // Molecular weight in the corona
  const CFreal mH = 1.6733e-27; // kg
  const CFreal Q_rad = pow((*currState)[0]*1.67e-13/(mu_cor*mH),2)*pow(10,-17.73)*pow(T2,-2.0/3.0);
  




  //std::cout << "SOURCE TERM" << endl;

      // VERSION #3
      CFreal l0 = 6.9e8;
      CFreal B0 = 2.2e-4;
      CFreal mu0 = 1.2566e-6;
      CFreal rho0 = 1.67e-13;
      CFreal V0 = B0/std::sqrt(mu0*rho0);
      CFreal g0 = pow(V0,2)/l0;
      CFreal g = -(GMsun/r2)/g0;
      CFreal gx = g*x_dimless/r_dimless;
      CFreal gy = g*y_dimless/r_dimless;
      CFreal gz = g*z_dimless/r_dimless;
      
      
      // ---- DENSITY ---------------------------------------------------------
      source[0] = 0.0;






      // --- M O M E N T U M   D E N S I T Y ----------------------------------
          source[1] = 0.0;
          source[2] = 0.0;
          source[3] = 0.0;

      if (_gravity == 1){
	//cout << "Gravity enabled" << endl;
	source[1] += (*currState)[0]*gx*volumes[elementID];
	source[2] += (*currState)[0]*gy*volumes[elementID];
	source[3] += (*currState)[0]*gz*volumes[elementID];

    CFreal gtot = std::sqrt(gx*gx+gy*gy+gz*gz);
	//socket_gravity.getDataHandle()[elementID*4]   = (*currState)[0]*gx*volumes[elementID]; // comment out *vol... to print the force density
	//socket_gravity.getDataHandle()[elementID*4+1] = (*currState)[0]*gy*volumes[elementID];
	//socket_gravity.getDataHandle()[elementID*4+2] = (*currState)[0]*gz*volumes[elementID];

	socket_gravity.getDataHandle()[elementID*4]   = gtot; // Just a test
	socket_gravity.getDataHandle()[elementID*4+1] = gtot*volumes[elementID];
	socket_gravity.getDataHandle()[elementID*4+2] = (*currState)[0]*gtot;
    socket_gravity.getDataHandle()[elementID*4+3] = (*currState)[0]*gtot*volumes[elementID]; // remove vol...
      }
      
      source[4] = 0.0;
      source[5] = 0.0;
      source[6] = 0.0;
      
      CFreal Vx = (*currState)[1];
      CFreal Vy = (*currState)[2];
      CFreal Vz = (*currState)[3];
      CFreal Vdotg = Vx*gx + Vy*gy + Vz*gz; // dimensionless
      
      // --- E N E R G Y   D E N S I T Y --------------------------------------
      source[7] = 0.0;
      if (_gravity == 1){
        source[7] += (*currState)[0]*Vdotg*volumes[elementID];
	
	
      }
      
      if (_Manchester == 1) {
        //cout << "Manchester heating term enabled." << endl;
        source[7] += Q*nondimconst*volumes[elementID];
      }

      if (_RadiativeLossTerm == 1) {
        //cout << "Radiative loss term enabled." << endl;
        source[7] -= Q_rad*nondimconst*volumes[elementID];
      }

      source[8] = 0.0;
  CFLog(DEBUG_MAX, "MHDConsACASourceTerm::computeSource() => END\n");
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<Framework::BaseDataSocketSource> >
MHDConsACASourceTerm::providesSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSource> > result = ComputeSourceTermFVMCC::providesSockets();
  result.push_back(&socket_gravity);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

} // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
