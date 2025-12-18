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
  _gradP(),
  _gradBx(),
  _gradBy(),
  _gradBz(),
  _gradRho(),
  _gradVx(),
  _gradVy(),
  _gradVz(),
  _gradEpsP(),
  _gradEpsM(),
  socket_gravity("gravity"),
  socket_heating("Manchester"),
  socket_radiativeloss("RadiativeLossTerm"),
 // socket_zp("zplus"),
 // socket_zm("zminus"),
  socket_wavepressure("wavepressure"),
  socket_divBCellCenter("divBCellCenter")
{
  addConfigOptionsTo(this);
  _gravity = 0; 
  setParameter("gravity",&_gravity);
  setParameter("PevtsovHeating",&_PevtsovHeating);
  setParameter("PevtsovHeatingFactor",&_PevtsovHeatingFactor);
  _Manchester = 0;
  setParameter("Manchester",&_Manchester);
  setParameter("ManchesterHeatingAmplitude",&_ManchesterHeatingAmplitude);
  setParameter("ManchesterSigma",&_ManchesterSigma);
  setParameter("Qh4H_const", &_Qh4H_const);
  setParameter("Qlio_AR", &_Qlio_AR);
  setParameter("P0_W", &_P0_W);
  setParameter("alfven_pressure", &_alfven_pressure);
  setParameter("Qh2_activate", &_Qh2_activate);
  setParameter("Qh3_activate", &_Qh3_activate);
  setParameter("Qh4_activate", &_Qh4_activate);
  setParameter("Qh_lio_activate", &_Qh_lio_activate);
  setParameter("R1", &_R1);
  setParameter("R2R3", &_R2R3);
  
  setParameter("Q_w", &_Q_w);
  setParameter("divQ",&_divQ);
  setParameter("divQConductivity",&_divQConductivity);
  setParameter("divQalphaCollisionless",&_divQalphaCollisionless);
  setParameter("ViscosityAndResistivity",&_ViscosityAndResistivity);
  setParameter("Viscosity",&_Viscosity);
  setParameter("Resistivity",&_Resistivity);
  setParameter("wavepressure", &_wave_pressure);  
  _RadiativeLossTerm = 0;
  setParameter("RadiativeLossTerm",&_RadiativeLossTerm);
  _rotation = 0;
  setParameter("Rotation",&_rotation);
//  setParameter("zplus",&_zplus);
//  setParameter("zminus",&_zminus); 
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
  options.addConfigOption< CFreal >("Qh4H_const", "Heating coefficient for the Qh4 heating function");
  options.addConfigOption< CFreal >("Qlio_AR", "Heating coefficient for the Qhlio heating function");
  options.addConfigOption< CFreal >("P0_W","amplitude for wave pressure");
  options.addConfigOption< CFint >("alfven_pressure", "flag for switching on alfven pressure");
  options.addConfigOption< CFint >("R1", "flag for activating R1 term of elsasser euqations");
  options.addConfigOption< CFint >("R2R3", "flag to activate all terms on the right side of the Elsasse equations");
 
  options.addConfigOption< CFint >("Q_w", "add heating to the energy equation from the WTD model");
  options.addConfigOption< CFint >("Qh2_activate", "whether Qh2 heating term should be activated");
  options.addConfigOption< CFint >("Qh3_activate", "whether Qh3 heating term should be activated");
  options.addConfigOption< CFint >("Qh4_activate", "whether Qh4 heating term should be activated");
  options.addConfigOption< CFint >("Qh_lio_activate", "whether Qh_lio heating term should be activated");
  options.addConfigOption< CFreal >("ManchesterSigma","Sigma value of Manchester heating term.");
  options.addConfigOption< CFint >("divQ","Switch on heat conduction.");
  options.addConfigOption< CFreal >("divQConductivity","Conductivity of heat conduction.");
  options.addConfigOption< CFreal >("divQalphaCollisionless","alpha value of heat conductivity in collisionless regime.");
  options.addConfigOption< CFint >("ViscosityAndResistivity","Switch on viscosity and resisitivity.");
  options.addConfigOption< CFreal >("Viscosity","nu value.");
  options.addConfigOption< CFreal >("Resistivity","eta value.");
  options.addConfigOption< CFint >("RadiativeLossTerm","Switch on optically thin approximation for radiation losses.");
  options.addConfigOption< CFint >("wavepressure", "save approximated wavepressure term.");
  options.addConfigOption< CFint >("Rotation","Switch on rotation (Coriolis and centrifugal forces).");
//  options.addConfigOption< CFint >("zplus", "z plus term from WKB");
//  options.addConfigOption< CFint >("zminus", "z minus term from WKB");
  
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
  
  socket_heating.getDataHandle().resize(socket_volumes.getDataHandle().size());
  socket_radiativeloss.getDataHandle().resize(socket_volumes.getDataHandle().size());
  socket_wavepressure.getDataHandle().resize(socket_volumes.getDataHandle().size());
  socket_divBCellCenter.getDataHandle().resize(socket_volumes.getDataHandle().size());
  
  DataHandle<CFreal> normals = this->socket_normals.getDataHandle();
  _gradP.resize(PhysicalModelStack::getActive()->getDim(),0.);
  _gradBx.resize(PhysicalModelStack::getActive()->getDim(),0.);
  _gradBy.resize(PhysicalModelStack::getActive()->getDim(),0.);
  _gradBz.resize(PhysicalModelStack::getActive()->getDim(),0.);
  _gradRho.resize(PhysicalModelStack::getActive()->getDim(),0.);
  _gradVx.resize(PhysicalModelStack::getActive()->getDim(),0.);
  _gradVy.resize(PhysicalModelStack::getActive()->getDim(),0.);
  _gradVz.resize(PhysicalModelStack::getActive()->getDim(),0.);
  _gradEpsP.resize(PhysicalModelStack::getActive()->getDim(),0.);
  _gradEpsM.resize(PhysicalModelStack::getActive()->getDim(),0.);

}

//////////////////////////////////////////////////////////////////////////////

void MHDConsACASourceTerm::computeSource(Framework::GeometricEntity *const element,
					 RealVector& source,
					 RealMatrix& jacobian)
{
  CFLog(DEBUG_MAX, "MHDConsACASourceTerm::computeSource() => START\n");
  DataHandle<CFreal> normals = this->socket_normals.getDataHandle(); 
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
  const CFuint startID = elementID*dim;

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
  
  if (_gravity == 0 && _Manchester == 0 && _RadiativeLossTerm ==0) return;

  CFLog(DEBUG_MIN, "gravity = " << _gravity << ", Manchester = " << _Manchester << ", RadiativeLossTerm = " << _RadiativeLossTerm << "\n");  
  
  // Define all the units and constants inside the code
  const CFreal mu_cor = 1.27; // Molecular weight in the corona 
  const CFreal mH = 1.6733e-27; // kg
  const CFreal pRef = 0.03851; // taken from pluto (probably), N/m^2
  const CFreal k = 1.3806503e-23; // Boltzmann constant in SI
  const CFreal mp = 1.67262158e-27; // mass of proton in kg
  const CFreal me = 9.10938188e-31; // mass of electron in SI
  const CFreal RSun = 6.9551e8; // m

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
  CFreal rxy = std::sqrt(x*x+y*y);

  const CFreal B0 = 2.2e-4;
  const CFreal mu0 = 1.2566e-6;
  const CFreal rhoRef = 1.67e-13; // kg/m^3
  const CFreal vRef = B0/std::sqrt(mu0*rhoRef);
  //const CFreal vRef = 480363.085276; //m/s // from the formula B0/std::sqrt(mu0*rho0) is 481,514 m/s
  const CFreal GMsun = 1.327474512e20; // SI value
  const CFreal gsun = 274.00; // at the surface m/s^2
  const CFreal g0 = pow(vRef,2)/RSun; // this value is 332
 // const CFreal l0 = 6.9e8;
  //const CFreal g0 = pow(vRef,2)/l0; // this value is 332
  const CFreal g = -(GMsun/r2)/g0; // gravity term dimensionless
  const CFreal gx = g*x_dimless/r_dimless;
  const CFreal gy = g*y_dimless/r_dimless;
  const CFreal gz = g*z_dimless/r_dimless;
  
  CFreal T_ref = pRef*mu_cor*mH/(2.0*k*rhoRef); // Unit of T in coconut
  CFreal T = (*currState)[7]*pRef/(2*(*currState)[0]*rhoRef*k/(mu_cor*mH)); //temperature in the code
  CFreal density = (*currState)[0]*rhoRef; // temperature in the code
  CFreal Bx =(*currState)[4]*2.2e-4; // in T
  CFreal By =(*currState)[5]*2.2e-4;
  CFreal Bz =(*currState)[6]*2.2e-4; 
  CFreal Bnorm = std::sqrt(std::pow(Bx,2) + std::pow(By,2) + std::pow(Bz,2));
  CFreal Br = x/r*Bx + y/r*By + z/r*Bz;
  CFreal Btheta = (x*z)/(rxy*r)*Bx + (y*z)/(r*rxy)*By - rxy/r*Bz;
  CFreal Bphi = -y/rxy*Bx + x/rxy*By;
  CFreal Vx = (*currState)[1];
  CFreal Vy =(*currState)[2];
  CFreal Vz =(*currState)[3];

  //---------------------------------- Heating functions --------------------------------------------------------------------------

  // The Constnt value for heating function 4 from Down et al 2010. 
  const CFreal Qh4H_const = _Qh4H_const;  
  const CFreal Qlio_AR = _Qlio_AR;
  // Switches passed from the CFcase file to choose which heating profile toa ctivate 
  const CFint Qh2_bool = _Qh2_activate;
  const CFint Qh3_bool = _Qh3_activate;
  const CFint Qh4_bool = _Qh4_activate;
  const CFint Qh_lio_bool = _Qh_lio_activate;

  // Qh2 is the petsov heating term. Q = H*|B|. (1e-7/1e-6/1e-4) is to convert from cgs to SI. (*2.2e-4) is to have magnetic field in the units. 
  CFreal Qh2 = Qh4H_const*Bnorm*1e-7 / 1e-6 /1e-4;
  // Qh3 is the exponential heating function which is uniform around the Sun.
//  CFreal Qh3 = 7.28e-5 * std::exp(-(r - RSun)/(40e6)) * ( 1e-7 / 1e-6);
  CFreal Qh3 = 4.9e-5 * std::exp(-(r - RSun)/(0.7*RSun)) * ( 1e-7 / 1e-6);
  CFreal Qh_extra = 5.0e-7 * std::exp(-(r - RSun)/(0.7*RSun)) * ( 1e-7 / 1e-6);
  // Qh_reville is also the exponential function, in the form given in the Reville et al 2020 paper.
  CFreal Qh_reville = 8.0e-5 *std::pow((RSun/r),2)*std::exp(-(r - RSun)/RSun)* ( 1e-7 / 1e-4);
  // Qh4 is the combined function of qh2 and Qh3 
  CFreal Qh4 = Qh2* std::exp(-(r - RSun)/(0.7*RSun));
 // CFreal Qh4 = Qh2* std::exp(-(r - RSun)/(40e6));


  CFreal lambda_0 = 0.7*RSun;
  CFreal H_0_exp = 4.9128e-7*(1e-7/1e-6);
  CFreal H_exp = H_0_exp*std::exp(-(r-RSun)/(lambda_0));


  CFreal H_0_QS = 1.18e-5*(1e-7/1e-6);
//  CFreal Br_c = 1.0e-4;
  CFreal Br_c = 0.55e-4;
  CFreal f_QS = 0.5*(1+std::tanh((1.7-r/RSun)/0.1))*std::exp(-(r/RSun-1)/0.2);
  CFreal H_QS = H_0_QS*f_QS*(Btheta*Btheta+Bphi*Bphi)/(Bnorm*(std::abs(Br)+Br_c));


  CFreal H_0_AR = 1.87e-5*(1e-7/1e-6);
  CFreal B_0 = 1e-4;
  CFreal g_AR = 0.5*(1+std::tanh((Bnorm/(1.0e-4)-Qlio_AR)/3.97));
  CFreal H_AR = H_0_AR*g_AR*std::pow((Bnorm/B_0), 1.2);

  // CFreal H_CH = (H_exp  + H_AR); // OLD
  CFreal H_CH = (H_exp + H_QS + H_AR); // NEW -> this will change the residuals for the regression case
  // Manchester heating term
  const CFreal theta = std::atan2(std::sqrt(x_dimless*x_dimless + y_dimless*y_dimless),z_dimless);
  const CFreal th1 = (MathTools::MathConsts::CFrealPi()/180.0)*17.5;
  const CFreal th2 = (MathTools::MathConsts::CFrealPi()/180.0)*61.5;
  CFreal critTheta = 0.;
  if (r_dimless <= 7.0) {
        critTheta = std::asin(std::sqrt(std::sin(th1)*std::sin(th1) + std::cos(th1)*std::cos(th1)*(r_dimless-1.0)/8.0));
  }
  if ((r_dimless > 7.0) && (r_dimless < 30.0)) {
        critTheta = std::asin(std::sqrt(std::sin(th2)*std::sin(th2) + std::cos(th2)*std::cos(th2)*(r_dimless-7.0)/40.0));
  }

  // volumetric heating amplitude (J/(kg*s*K))
  const CFreal q0 = _ManchesterHeatingAmplitude;
  CFreal T0_manch = 0.;
  CFreal sigma0 = 0.;
 
  // target temperature and the heating scale height
  if (std::fabs(0.5*MathTools::MathConsts::CFrealPi()-theta) < critTheta) {
    T0_manch = 2.63e6;
    CFreal new_theta= std::fabs(0.5*MathTools::MathConsts::CFrealPi()-theta);
    sigma0 = 4.5*(2.0-(std::sin(new_theta)*std::sin(new_theta)/(std::sin(critTheta)*std::sin(critTheta))));
  }
  if (std::fabs(0.5*MathTools::MathConsts::CFrealPi()-theta) >= critTheta) {
    T0_manch=1.5e6;
    sigma0 = 4.5;
  }

  const CFreal Q_manchester = (*currState)[0]*1.67e-13*q0*(T0_manch-T)*exp(-(r_dimless-1.0)*(r_dimless-1.0)/(sigma0*sigma0));


 
  // Adimentionalizing constant for all heating terms:
  const CFreal nondimconst = RSun/(rhoRef*pow(vRef,3));

 
  //const CFreal nondimconst = RSun/(rhoRef*pow(vRef,3))*mu_cor*mH/k/density;
  //-------------------------------- Radiative Losses -----------------------------------------------------------------------------

  // Definition according to Rosner 1978

  CFreal Q_rad = 0.0;

  if ((T < pow(10,4.6))) {
    Q_rad = std::pow(1e-6*density/(mu_cor*mH),2)*std::pow(10,-21.85) * ( 1e-7 / 1e-6); 
  }

  else if ((T >= pow(10, 4.6)) && (T < pow(10,4.9))) {
    Q_rad = std::pow(1e-6*density/(mu_cor*mH),2)*std::pow(10,-31.0) * ( 1e-7 / 1e-6) * std::pow(T,2.0);
  } 

  else if ((T >= pow(10, 4.9)) && (T < pow(10,5.4))) {
    Q_rad = std::pow(1e-6*density/(mu_cor*mH),2)*std::pow(10,-21.2) * ( 1e-7 / 1e-6);
  }

  else if ((T >= pow(10, 5.4)) && (T < pow(10,5.75))) {
    Q_rad = std::pow(1e-6*density/(mu_cor*mH),2)*std::pow(10,-10.4) * ( 1e-7 / 1e-6) * std::pow(T,-2.0);
  }

  else if ((T >= pow(10, 5.75)) && (T < pow(10,6.3))) {
    Q_rad = std::pow(1e-6*density/(mu_cor*mH),2)*std::pow(10,-21.94) * ( 1e-7 / 1e-6);
  }

  else if ((T >= pow(10, 6.3))) {
    Q_rad = std::pow(1e-6*density/(mu_cor*mH),2)*std::pow(10,-17.73) * ( 1e-7 / 1e-6) * std::pow(T,-2.0/3.0);
  }



  // Definition acording to Athay 1986

  CFreal Q_rad_athay = 0.0;
  Q_rad_athay = std::pow(1e-6*density/(mu_cor*mH),2)*(1e-7/1e-6)*1e-22*(0.4*exp(-30.0*pow((log10(T)-4.6),2.0))+4.0*exp(-20.0*pow((log10(T)-4.9),2.0))+ 4.5*exp(-16.0*pow((log10(T)-5.35),2.0))+2.0*exp(-4.0*pow((log10(T)-6.1),2.0)));

 



    // Compute gradients for the approximated pressure
    CFreal P0_W = _P0_W;
     
    // Gradients for p = sin^2(theta)
    CFreal common_der_expr = 1.0/(x_dimless*x_dimless+y_dimless*y_dimless+z_dimless*z_dimless)/pow((x_dimless*x_dimless+y_dimless*y_dimless),1.0/3.0); 
    CFreal common_derivative = pow(common_der_expr,3.0/2.0);
   // CFreal dP_dx =( P0_W*x_dimless*pow(z_dimless,2)*common_derivative);
   // CFreal dP_dy =( P0_W*y_dimless*pow(z_dimless,2)*common_derivative);
   // CFreal dP_dz =(- P0_W*z_dimless*pow(x_dimless*x_dimless+y_dimless*y_dimless,1.0/2.0)*common_derivative);
    
  

    // The pressure function of which we are calculating gradients analytically. 
//    CFreal wave_pressure_approx = P0_W*(x_dimless*x_dimless + y_dimless*y_dimless)/pow(r2_dimless,3);
//    CFreal wave_pressure_approx = 0.0; //P0_W* 1.0/r_dimless/r_dimless;

  

      // Gradients for p = 1/r^2
//    CFreal dP_dx = -2.0*P0_W*x_dimless*pow(r_dimless,-2.0);
//    CFreal dP_dy = -2.0*P0_W*y_dimless*pow(r_dimless,-2.0);
//    CFreal dP_dz = -2.0*P0_W*z_dimless*pow(r_dimless,-2.0);

    
    // Gradients for P = sin^2(theta)/ r^2
   // CFreal dP_dx = 2.0*P0_W*x_dimless*(pow(z_dimless,2.0)-pow(y_dimless,2.0)-pow(x_dimless,2))/pow(r_dimless,6.0);
   // CFreal dP_dy = 2.0*P0_W*y_dimless*(pow(z_dimless,2.0)-pow(y_dimless,2.0)-pow(x_dimless,2))/pow(r_dimless,6.0); 
   // CFreal dP_dz = -4.0*P0_W*z_dimless*(pow(x_dimless,2.0)+pow(y_dimless,2.0))/pow(r_dimless,6);
///    CFreal dP_dx = 2.0*P0_W*x_dimless*(pow(z_dimless,2.0)-4.0*pow(y_dimless,2.0)-4.0*pow(x_dimless,2))/pow(r_dimless,8);
//    CFreal dP_dy = 2.0*P0_W*y_dimless*(pow(z_dimless,2.0)-4.0*pow(y_dimless,2.0)-4.0*pow(x_dimless,2))/pow(r_dimless,8);
//    CFreal dP_dz = -6.0*P0_W*z_dimless*(pow(x_dimless,2.0)+pow(y_dimless,2.0))/pow(r_dimless,8);

//    CFreal lambd = 4.4/std::sqrt((*currState)[4]*(*currState)[4] + (*currState)[5]*(*currState)[5] + (*currState)[6]*(*currState)[6]);
//    CFreal heating_alfwen = wave_pressure_approx*std::sqrt(8.0/(*currState)[0]*wave_pressure_approx);


    // Calcualting the gradient of variable pressure 

    CFreal lambda_perp = 0;
    const CFuint PID = 7;
    const CFuint gradPID = elementID*nbEqs + PID;
    _gradP[XX] = this->m_ux[gradPID];
    _gradP[YY] = this->m_uy[gradPID];
    _gradP[ZZ] = this->m_uz[gradPID];
    CFreal n_x = normals[startID];
    CFreal n_y = normals[startID+1];
    CFreal n_z = normals[startID+2];
    CFreal Qw = 0.0;
    
    CFreal dP_dx = 0.0;
    CFreal dP_dy = 0.0;
    CFreal dP_dz = 0.0;
    CFreal wave_pressure_approx = 0.0;
    if (nbEqs == 11) {
      source[9]=0;
      source[10]=0;
      // Calculating the source terms for Eps equations
      CFreal z0 = 30000; // This hardcoding should be removed... distances/coordinates are non-dimensional!!! 
      const CFuint BXID = 4;
      const CFuint BYID = 5;
      const CFuint BZID = 6;
      const CFuint RHOID = 0;
      const CFuint VXID = 1;
      const CFuint VYID = 2;
      const CFuint VZID = 3;
      const CFuint gradBXID = elementID*nbEqs + BXID;
      const CFuint gradBYID = elementID*nbEqs + BYID;
      const CFuint gradBZID = elementID*nbEqs + BZID;
      const CFuint gradRhoID = elementID*nbEqs + RHOID;
      const CFuint gradVXID = elementID*nbEqs + VXID;
      const CFuint gradVYID = elementID*nbEqs + VYID;
      const CFuint gradVZID = elementID*nbEqs + VZID;
      //   std::cout << "definition of grad vars " << gradBXID << gradRhoID << gradVXID << endl;
      _gradBx[XX] = this->m_ux[gradBXID];
      _gradBx[YY] = this->m_uy[gradBXID];
      _gradBx[ZZ] = this->m_uz[gradBXID];
      _gradBy[XX] = this->m_ux[gradBYID];
      _gradBy[YY] = this->m_uy[gradBYID];
      _gradBy[ZZ] = this->m_uz[gradBYID];
      _gradBz[XX] = this->m_ux[gradBZID];
      _gradBz[YY] = this->m_uy[gradBZID];
      _gradBz[ZZ] = this->m_uz[gradBZID];
      _gradRho[XX] = this->m_ux[gradRhoID];
      _gradRho[YY] = this->m_uy[gradRhoID];
      _gradRho[ZZ] = this->m_uz[gradRhoID];
      _gradVx[XX] = this->m_ux[gradVXID];
      _gradVx[YY] = this->m_uy[gradVXID];
      _gradVx[ZZ] = this->m_uz[gradVXID];
      _gradVy[XX] = this->m_ux[gradVYID];
      _gradVy[YY] = this->m_uy[gradVYID];
      _gradVy[ZZ] = this->m_uz[gradVYID];
      _gradVz[XX] = this->m_ux[gradVZID];
      _gradVz[YY] = this->m_uy[gradVZID];
      _gradVz[ZZ] = this->m_uz[gradVZID];
      //    std::cout << "d/dx bx " << _gradBx[XX] << " Bx " << Bx << " Bxref * grad " << _gradBx[XX]*2.2e-4 << endl;
      const CFreal dvaxdx = 1/std::sqrt(mu0*density)*(_gradBx[XX]*B0 - 0.5*Bx/density*_gradRho[XX]*rhoRef)/RSun;
      const CFreal dvaxdy = 1/std::sqrt(mu0*density)*(_gradBx[YY]*B0 - 0.5*Bx/density*_gradRho[YY]*rhoRef)/RSun;
      const CFreal dvaxdz = 1/std::sqrt(mu0*density)*(_gradBx[ZZ]*B0 - 0.5*Bx/density*_gradRho[ZZ]*rhoRef)/RSun;

      const CFreal dvaydx = 1/std::sqrt(mu0*density)*(_gradBy[XX]*B0 - 0.5*By/density*_gradRho[XX]*rhoRef)/RSun;
      const CFreal dvaydy = 1/std::sqrt(mu0*density)*(_gradBy[YY]*B0 - 0.5*By/density*_gradRho[YY]*rhoRef)/RSun;
      const CFreal dvaydz = 1/std::sqrt(mu0*density)*(_gradBy[ZZ]*B0 - 0.5*By/density*_gradRho[ZZ]*rhoRef)/RSun;

      const CFreal dvazdx = 1/std::sqrt(mu0*density)*(_gradBz[XX]*B0 - 0.5*Bz/density*_gradRho[XX]*rhoRef)/RSun;
      const CFreal dvazdy = 1/std::sqrt(mu0*density)*(_gradBz[YY]*B0 - 0.5*Bz/density*_gradRho[YY]*rhoRef)/RSun;
      const CFreal dvazdz = 1/std::sqrt(mu0*density)*(_gradBz[ZZ]*B0 - 0.5*Bz/density*_gradRho[ZZ]*rhoRef)/RSun;


      const CFreal dvxdx = _gradVx[XX]*vRef/RSun;
      const CFreal dvydy = _gradVy[YY]*vRef/RSun;
      const CFreal dvzdz = _gradVz[ZZ]*vRef/RSun;
      const CFreal vxdim = Vx*vRef;
      const CFreal vydim = Vy*vRef;
      const CFreal vzdim = Vz*vRef;
      const CFreal vax = Bx/std::sqrt(mu0*density);
      const CFreal vay = By/std::sqrt(mu0*density);
      const CFreal vaz = Bz/std::sqrt(mu0*density);
      const CFreal va_length = Bnorm/std::sqrt(mu0*density);
      const CFreal R1_p = 1./(4.*density) *((vxdim - vax)*_gradRho[XX] + (vydim - vay)*_gradRho[YY]+ (vzdim - vaz)*_gradRho[ZZ])*rhoRef/RSun;
      const CFreal R1_m = 1./(4.*density) *((vxdim + vax)*_gradRho[XX] + (vydim + vay)*_gradRho[YY]+ (vzdim + vaz)*_gradRho[ZZ])*rhoRef/RSun;
//      const CFreal R2_p = 1./2.*((vxdim/vax-1)*dvaxdx+(vydim/vay-1)*dvaydy+(vzdim/vaz-1)*dvazdz);
//      const CFreal R2_m = 1./2.*((vxdim/vax+1)*dvaxdx+(vydim/vay+1)*dvaydy+(vzdim/vaz+1)*dvazdz);
      
      const CFreal R2_1 = vax*dvaxdx+vay*dvaxdy+vaz*dvaxdz;
      const CFreal R2_2 = vax*dvaydx+vay*dvaydy+vaz*dvaydz;
      const CFreal R2_3 = vax*dvazdx+vay*dvazdy+vaz*dvazdz;
      
      const CFreal R2_p = 1./(2.*va_length*va_length)*(R2_1*(vxdim-vax) + R2_2*(vydim-vay) + R2_3*(vzdim-vaz));
      const CFreal R2_m = 1./(2.*va_length*va_length)*(R2_1*(vxdim+vax) + R2_2*(vydim+vay) + R2_3*(vzdim+vaz));
      const CFreal second_term_plus = dvxdx + dvydy + dvzdz + dvaxdx + dvaydy + dvazdz;
      const CFreal second_term_minus = dvxdx + dvydy + dvzdz - dvaxdx - dvaydy - dvazdz;
      lambda_perp = 0.0000005*RSun * std::sqrt(1e-5/Bnorm); // Should be 1e-4 
      CFreal third_term_plus = std::abs((*currState)[10])*vRef/(2.0*lambda_perp);
      CFreal third_term_minus = std::abs((*currState)[9])*vRef/(2.0*lambda_perp);
         
//      std::cout << "R1_p*zp " << R1_p*(*currState)[9]  << "  " << endl;
//      std::cout << "thirs term " << third_term_plus* (*currState)[9] << "  " <<endl;
      CFreal source_plus = (*currState)[9] *second_term_plus;
      CFreal source_minus = (*currState)[10] *second_term_minus;
      if (_R1==1 && _R2R3==0){
      	source_plus = (*currState)[9] * (R1_p + second_term_plus);
	source_minus = (*currState)[10] * (R1_m + second_term_minus);
      } else if (_R1==1 && _R2R3 == 1){
        source_plus = (*currState)[9] * (R1_p + second_term_plus -third_term_plus) +(*currState)[10]*R2_p;
	source_minus = (*currState)[10] * (R1_m + second_term_minus-third_term_minus)+ (*currState)[9]*R2_m;
      } 
    //  CFreal source_plus = (*currState)[9] * (R1_p + second_term_plus -third_term_plus) +(*currState)[10]*R2_p;
    //  CFreal source_minus = (*currState)[10] * (R1_m + second_term_minus-third_term_minus)+ (*currState)[9]*R2_m;
      //const CFreal source_plus = (*currState)[9] * (second_term_plus);
      //const CFreal source_minus = (*currState)[10] * (second_term_minus);
      source[9]  += source_plus*RSun/vRef*volumes[elementID]; 
      source[10] += source_minus*RSun/vRef*volumes[elementID];
      // CFLog(INFO, "AFTER source 9 " << source[9] << " source 10 before " << source[10] << "\n"); 
      //    source[9] += 1000.*volumes[elementID];
      //    source[10] += 1000.*volumes[elementID];
      //    source[9] = 0.0;
      //    source[10] = 0.0;
    
    // VERSION #3   
    CFreal epsP = (*currState)[9];
    CFreal epsM = (*currState)[10];   
    Qw =0.01*density *(std::abs(epsM)*epsP*epsP + std::abs(epsP)*epsM*epsM)*vRef*vRef*vRef/(4.0*lambda_perp);
    const CFuint EpsPID = 9;
    const CFuint EpsMID = 10;
    const CFuint gradEpsPID = elementID*nbEqs + EpsPID;
    const CFuint gradEpsMID = elementID*nbEqs + EpsMID;
    _gradEpsP[XX] = this->m_ux[gradEpsPID];
    _gradEpsP[YY] = this->m_uy[gradEpsPID];
    _gradEpsP[ZZ] = this->m_uz[gradEpsPID];
    _gradEpsM[XX] = this->m_ux[gradEpsMID];
    _gradEpsM[YY] = this->m_uy[gradEpsMID];
    _gradEpsM[ZZ] = this->m_uz[gradEpsMID];
    dP_dx = _gradRho[XX]*std::pow((epsM-epsP),2)/8.0 + (*currState)[0]/8.0*2.0*(epsM-epsP)*(_gradEpsM[XX]-_gradEpsP[XX]);
    dP_dy = _gradRho[YY]*std::pow((epsM-epsP),2)/8.0 + (*currState)[0]/8.0*2.0*(epsM-epsP)*(_gradEpsM[YY]-_gradEpsP[YY]);
    dP_dz = _gradRho[ZZ]*std::pow((epsM-epsP),2)/8.0 + (*currState)[0]/8.0*2.0*(epsM-epsP)*(_gradEpsM[ZZ]-_gradEpsP[ZZ]);
    wave_pressure_approx = density*(epsM-epsP)*(epsM-epsP)*vRef*vRef/8.0;
    }
   
    if (_Q_w == 0){
      Qw = 0;
    } 

  //  CFreal lambda_perp = 0.01*RSun * std::sqrt(5e-5/Bnorm); 
 //   Qw =density *(std::abs(epsM)*epsP*epsP + std::abs(epsP)*epsM*epsM)*vRef*vRef*vRef/(4.0*lambda_perp); 
//    const CFuint EpsPID = 9;
//    const CFuint EpsMID = 10;
//    const CFuint gradEpsPID = elementID*nbEqs + EpsPID;
//    const CFuint gradEpsMID = elementID*nbEqs + EpsMID;
//    _gradEpsP[XX] = this->m_ux[gradEpsPID];
//    _gradEpsP[YY] = this->m_uy[gradEpsPID];
 //   _gradEpsP[ZZ] = this->m_uz[gradEpsPID];
//    _gradEpsM[XX] = this->m_ux[gradEpsMID];
//    _gradEpsM[YY] = this->m_uy[gradEpsMID];
//    _gradEpsM[ZZ] = this->m_uz[gradEpsMID];
//    CFreal dP_dx = _gradRho[XX]*std::pow((epsM-epsP),2)/8.0 + (*currState)[0]/8.0*2.0*(epsM-epsP)*(_gradEpsM[XX]-_gradEpsP[XX]);
//    CFreal dP_dy = _gradRho[YY]*std::pow((epsM-epsP),2)/8.0 + (*currState)[0]/8.0*2.0*(epsM-epsP)*(_gradEpsM[YY]-_gradEpsP[YY]);
//    CFreal dP_dz = _gradRho[ZZ]*std::pow((epsM-epsP),2)/8.0 + (*currState)[0]/8.0*2.0*(epsM-epsP)*(_gradEpsM[ZZ]-_gradEpsP[ZZ]);
//    CFreal  wave_pressure_approx = density*(epsM-epsP)*(epsM-epsP)*vRef*vRef/8.0; 
    
    // ---- DENSITY ---------------------------------------------------------
    source[0] = 0.0;
    
    // --- M O M E N T U M   D E N S I T Y ----------------------------------
    source[1] = 0.0;
    source[2] = 0.0;
    source[3] = 0.0;
    
    if (_gravity == 1){
      source[1] += (*currState)[0]*gx*volumes[elementID];
      source[2] += (*currState)[0]*gy*volumes[elementID];
      source[3] += (*currState)[0]*gz*volumes[elementID];
      //	std::cout << "source momentum 1 " << source[1] << " density " << (*currState)[0] << " addition  " <<  (*currState)[0]*gx*volumes[elementID] << endl;
      socket_gravity.getDataHandle()[elementID*4]   = (*currState)[0]*gx; //*volumes[elementID];
      socket_gravity.getDataHandle()[elementID*4+1] = (*currState)[0]*gy; //*volumes[elementID];
      socket_gravity.getDataHandle()[elementID*4+2] = (*currState)[0]*gz; //*volumes[elementID];
    } 
    
    socket_wavepressure.getDataHandle()[elementID] = 0.0;
    if (_alfven_pressure == 1) {
      source[1] -= dP_dx*volumes[elementID];
      source[2] -= dP_dy*volumes[elementID];
      source[3] -= dP_dz*volumes[elementID]; 
      socket_wavepressure.getDataHandle()[elementID] = wave_pressure_approx;
    }

    // Rotation source terms (Coriolis and centrifugal forces)
    if (_rotation == 1) {
      // Sidereal rotation rate [rad/s]
      const CFreal Omega_phys = 2.97e-6;
      // Non-dimensional rotation rate
      const CFreal Omega_nd = Omega_phys * RSun / vRef;
      
      // Coriolis acceleration: a_cor = -2 * Omega x v
      // With Omega along z-axis: Omega x v = (Omega_z * Vy, -Omega_z * Vx, 0)
      const CFreal cor_x = -2.0 * Omega_nd * Vy;
      const CFreal cor_y =  2.0 * Omega_nd * Vx;
      const CFreal cor_z =  0.0;
      
      // Centrifugal acceleration: a_cent = -Omega^2 * (x, y, 0)
      const CFreal cent_x = -Omega_nd * Omega_nd * x_dimless;
      const CFreal cent_y = -Omega_nd * Omega_nd * y_dimless;
      const CFreal cent_z =  0.0;
      
      // Add to momentum equations (note: source = rho * a for conservative form)
      source[1] -= (*currState)[0] * (cor_x + cent_x) * volumes[elementID];
      source[2] -= (*currState)[0] * (cor_y + cent_y) * volumes[elementID];
      source[3] -= (*currState)[0] * (cor_z + cent_z) * volumes[elementID];
    }
      
    source[4] = 0.0;
    source[5] = 0.0;
    source[6] = 0.0;
    
    // --- E N E R G Y   D E N S I T Y --------------------------------------
    source[7] = 0.0;
    
    if (_gravity == 1){
      CFreal Vdotg = Vx*gx + Vy*gy + Vz*gz; // dimensionless
        source[7] += (*currState)[0]*Vdotg*volumes[elementID];
	socket_gravity.getDataHandle()[elementID*4+3] = (*currState)[0]*Vdotg; //*volumes[elementID];
    }

    // Rotation contribution to energy equation
    if (_rotation == 1) {
      const CFreal Omega_phys = 2.97e-6;
      const CFreal Omega_nd = Omega_phys * RSun / vRef;
      
      // Centrifugal acceleration
      const CFreal cent_x = -Omega_nd * Omega_nd * x_dimless;
      const CFreal cent_y = -Omega_nd * Omega_nd * y_dimless;
      const CFreal cent_z =  0.0;
      
      // Energy source: v dot a_cent (centrifugal work)
      const CFreal Vdotrot = Vx * cent_x + Vy * cent_y + Vz * cent_z;
      source[7] -= (*currState)[0] * Vdotrot * volumes[elementID];
    }
    
    if (_Manchester == 1) {
      
      CFreal result_flux = 0.0;
      if (Qh2_bool == 1) {
	result_flux = Qh2;
      } else if (Qh3_bool == 1) {
	result_flux = Qh3;
      } else if (Qh4_bool == 1) {
	result_flux = Qh4;
      } else if (Qh_lio_bool == 1) {
	result_flux = H_CH;
      }
      source[7] += ((result_flux+Qw)*nondimconst)*volumes[elementID];
      
      socket_heating.getDataHandle()[elementID]=result_flux+Qw;
    }
    
    if (_RadiativeLossTerm == 1) {
      
      
      source[7] -= Q_rad*nondimconst*volumes[elementID];
      
      socket_radiativeloss.getDataHandle()[elementID] = Q_rad;
    }
    
    source[8] = 0.0;
    // Eps modification 
    //    if (states->size() == 11) {
    //    source[9] = source_plus/vRef*volumes[elementID];
    //    source[10] = source_minus/vRef*volumes[elementID];
    //    socket_zp.getDataHandle()[elementID] = source_plus/vRef;
    //    socket_zm.getDataHandle()[elementID] = source_minus/vRef;
    //    }
    CFLog(DEBUG_MAX, "MHDConsACASourceTerm::computeSource() => END\n");
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<Framework::BaseDataSocketSource> >
MHDConsACASourceTerm::providesSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSource> > result = ComputeSourceTermFVMCC::providesSockets();
  result.push_back(&socket_gravity);
  result.push_back(&socket_heating);
  result.push_back(&socket_radiativeloss);
  result.push_back(&socket_wavepressure);
  result.push_back(&socket_divBCellCenter);
  
//  result.push_back(&socket_zp);
//  result.push_back(&socket_zm);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

} // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
