#include "Framework/MethodStrategyProvider.hh"
#include "Framework/MeshData.hh"
#include "Framework/DataHandle.hh"
#include "Framework/GeometricEntity.hh"
#include "MHD/MHDTerm.hh"
#include "FiniteVolumeMHD/FiniteVolumeMHD.hh"
#include "FiniteVolumeMHD/MHDConsACAHWSourceTerm.hh"
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

MethodStrategyProvider<MHDConsACAHWSourceTerm,
		       CellCenterFVMData,
		       ComputeSourceTerm<CellCenterFVMData>,
		       FiniteVolumeMHDModule>
mHDConsACAHWSTFVMCCProvider("MHDConsACAHWST");

//////////////////////////////////////////////////////////////////////////////

MHDConsACAHWSourceTerm::MHDConsACAHWSourceTerm(const std::string& name) :
  ComputeSourceTermFVMCC(name),
  _gradP(),
  _gradBx(),
  _gradBy(),
  _gradBz(),
  _gradRho(),
  _gradVx(),
  _gradVy(),
  _gradVz(),
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
  setParameter("divQ",&_divQ);
  setParameter("divQConductivity",&_divQConductivity);
  setParameter("divQalphaCollisionless",&_divQalphaCollisionless);
  setParameter("Resistivity",&_Resistivity);
  setParameter("wavepressure", &_wave_pressure);  
  _RadiativeLossTerm = 0;
  setParameter("RadiativeLossTerm",&_RadiativeLossTerm);
  _deCompE = 0;
  setParameter("deCompEorNot", &_deCompE);
//  setParameter("zplus",&_zplus);
//  setParameter("zminus",&_zminus); 
  //_projectionIDs = vector<CFuint>();
  //setParameter("ProjectionIDs",&_projectionIDs);

}

//////////////////////////////////////////////////////////////////////////////

void MHDConsACAHWSourceTerm::defineConfigOptions(Config::OptionList& options)
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
  options.addConfigOption< CFint >("Qh2_activate", "whether Qh2 heating term should be activated");
  options.addConfigOption< CFint >("Qh3_activate", "whether Qh3 heating term should be activated");
  options.addConfigOption< CFint >("Qh4_activate", "whether Qh4 heating term should be activated");
  options.addConfigOption< CFint >("Qh_lio_activate", "whether Qh_lio heating term should be activated");
  options.addConfigOption< CFreal >("ManchesterSigma","Sigma value of Manchester heating term.");
  options.addConfigOption< CFint >("divQ","Switch on heat conduction.");
  options.addConfigOption< CFreal >("divQConductivity","Conductivity of heat conduction.");
  options.addConfigOption< CFreal >("divQalphaCollisionless","alpha value of heat conductivity in collisionless regime.");
  options.addConfigOption< CFreal >("Resistivity","eta value.");
  options.addConfigOption< CFint >("RadiativeLossTerm","Switch on optically thin approximation for radiation losses.");
  options.addConfigOption< CFint >("deCompEorNot", "Calculate source term caused by decomposed energy equation.");
  options.addConfigOption< CFint >("wavepressure", "save approximated wavepressure term.");
//  options.addConfigOption< CFint >("zplus", "z plus term from WKB");
//  options.addConfigOption< CFint >("zminus", "z minus term from WKB");
  
}

//////////////////////////////////////////////////////////////////////////////

MHDConsACAHWSourceTerm::~MHDConsACAHWSourceTerm()
{
}

//////////////////////////////////////////////////////////////////////////////

void MHDConsACAHWSourceTerm::setup()
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
}

//////////////////////////////////////////////////////////////////////////////

void MHDConsACAHWSourceTerm::computeSource(Framework::GeometricEntity *const element,
					 RealVector& source,
					 RealMatrix& jacobian)
{
  CFLog(DEBUG_MAX, "MHDConsACAHWSourceTerm::computeSource() => START\n");
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
  
  // const CFreal l0 = 6.9e8;
  const CFreal B0 = 2.2e-4;
  const CFreal mu0 = 1.2566e-6;
  const CFreal rhoRef = 1.67e-13; // kg/m^3
  const CFreal vRef = B0/std::sqrt(mu0*rhoRef);
  //const CFreal vRef = 480363.085276; //m/s // from the formula B0/std::sqrt(mu0*rho0) is 481,514 m/s
  const CFreal GMsun = 1.327474512e20; // SI value
  const CFreal gsun = 274.00; // at the surface m/s^2
  const CFreal g0 = pow(vRef,2)/RSun; // this value is 332
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
  CFreal Btheta = (x*z)/(rxy*r)*Bz + (y*z)/(r*rxy)*By - rxy/r*Bz;
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

  CFreal H_CH = (H_exp + H_QS + H_AR);
//  CFreal H_CH = (H_exp  + H_AR);
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
/*
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
*/

  // Definition according to Dere et al. 2019
  std::vector<CFreal> temperature;  
  std::vector<CFreal> RadCur; 
  std::vector<CFreal> logRadCur; 
  CFuint m_sampleTem=101;
  temperature.resize(m_sampleTem, 0.0);
  RadCur.resize(m_sampleTem, 0.0);
  logRadCur.resize(m_sampleTem, 0.0);
  // erg S^{-1} cm^{-3}
     logRadCur[0]= -23.00744648;
     logRadCur[1]= -22.55439580;
     logRadCur[2]= -22.15614458;
     logRadCur[3]= -21.83268267;
     logRadCur[4]= -21.64589156;
     logRadCur[5]= -21.61618463;
     logRadCur[6]= -21.68402965;
     logRadCur[7]= -21.79048499;
     logRadCur[8]= -21.87614836;
     logRadCur[9]= -21.91009489;
     logRadCur[10]= -21.89962945;
     logRadCur[11]= -21.86012091;
     logRadCur[12]= -21.79588002;
     logRadCur[13]= -21.71669877;
     logRadCur[14]= -21.62342304;
     logRadCur[15]= -21.52143350;
     logRadCur[16]= -21.41793664;
     logRadCur[17]= -21.33068312;
     logRadCur[18]= -21.27736608;
     logRadCur[19]= -21.25181197;
     logRadCur[20]= -21.24184538;
     logRadCur[21]= -21.25806092;
     logRadCur[22]= -21.27901426;
     logRadCur[23]= -21.27164622;
     logRadCur[24]= -21.24412514;
     logRadCur[25]= -21.21467016;
     logRadCur[26]= -21.19586057;
     logRadCur[27]= -21.18309616;
     logRadCur[28]= -21.18708664;
     logRadCur[29]= -21.24108811;
     logRadCur[30]= -21.35163999;
     logRadCur[31]= -21.45099674;
     logRadCur[32]= -21.48678240;
     logRadCur[33]= -21.47625353;
     logRadCur[34]= -21.45222529;
     logRadCur[35]= -21.43297363;
     logRadCur[36]= -21.42596873;
     logRadCur[37]= -21.42021640;
     logRadCur[38]= -21.41005040;
     logRadCur[39]= -21.40120949;
     logRadCur[40]= -21.40450378;
     logRadCur[41]= -21.42250820;
     logRadCur[42]= -21.44977165;
     logRadCur[43]= -21.47755577;
     logRadCur[44]= -21.51144928;
     logRadCur[45]= -21.55909092;
     logRadCur[46]= -21.63451202;
     logRadCur[47]= -21.73754891;
     logRadCur[48]= -21.85078089;
     logRadCur[49]= -21.95467702;
     logRadCur[50]= -22.03526908;
     logRadCur[51]= -22.08990945;
     logRadCur[52]= -22.11804503;
     logRadCur[53]= -22.12436006;
     logRadCur[54]= -22.11633856;
     logRadCur[55]= -22.10072681;
     logRadCur[56]= -22.08301995;
     logRadCur[57]= -22.06701918;
     logRadCur[58]= -22.05650548;
     logRadCur[59]= -22.05551733;
     logRadCur[60]= -22.06803389;
     logRadCur[61]= -22.10072681;
     logRadCur[62]= -22.15926677;
     logRadCur[63]= -22.24033216;
     logRadCur[64]= -22.32882716;
     logRadCur[65]= -22.40782324;
     logRadCur[66]= -22.47108330;
     logRadCur[67]= -22.51855737;
     logRadCur[68]= -22.54975089;
     logRadCur[69]= -22.57024772;
     logRadCur[70]= -22.58004425;
     logRadCur[71]= -22.58335949;
     logRadCur[72]= -22.58169871;
     logRadCur[73]= -22.57511836;
     logRadCur[74]= -22.56543110;
     logRadCur[75]= -22.55439580;
     logRadCur[76]= -22.54060751;
     logRadCur[77]= -22.52724355;
     logRadCur[78]= -22.51427857;
     logRadCur[79]= -22.49894074;
     logRadCur[80]= -22.48545225;
     logRadCur[81]= -22.47108330;
     logRadCur[82]= -22.45593196;
     logRadCur[83]= -22.44009337;
     logRadCur[84]= -22.42365865;
     logRadCur[85]= -22.40671393;
     logRadCur[86]= -22.38933984;
     logRadCur[87]= -22.37059040;
     logRadCur[88]= -22.35163999;
     logRadCur[89]= -22.33161408;
     logRadCur[90]= -22.31158018;
     logRadCur[91]= -22.29073004;
     logRadCur[92]= -22.26921772;
     logRadCur[93]= -22.24795155;
     logRadCur[94]= -22.22621356;
     logRadCur[95]= -22.20411998;
     logRadCur[96]= -22.18111459;
     logRadCur[97]= -22.15864053;
     logRadCur[98]= -22.13608262;
     logRadCur[99]= -22.11294562;
     logRadCur[100]= -22.08937560;   

  for (CFuint i = 0; i < m_sampleTem; ++i){
	  CFreal Exp_Fac=4.0+0.05*i;
      temperature[i]=std::pow(10,Exp_Fac); 
      RadCur[i]=std::pow(10,logRadCur[i]); 
  }
  
  CFreal Q_rad = 0.0;
  CFreal RadCur_Cur = 0.0;
  double Coefi1;
  double Coefi2;
  CFreal Temp;
  if ((T < pow(10,4.0))) {
    Q_rad = 0.0; 
  }
  else if(T > pow(10,9.0)){
	Q_rad = std::pow(1e-6*density/(mu_cor*mH),2)*RadCur[m_sampleTem-1] * ( 1e-7 / 1e-6);  
  }
  else{
	  for (CFuint i = 0; i < m_sampleTem-1; ++i){
		  if((T-temperature[i])*(temperature[i+1]-T)>=0.){
			  Temp=temperature[i+1]-temperature[i];
	          Coefi1=(temperature[i+1]-T)/Temp;
	          Coefi2=(T-temperature[i])/Temp;
			  RadCur_Cur=RadCur[i] * Coefi1 + RadCur[i+1]* Coefi2;
			  Q_rad = std::pow(1e-6*density/(mu_cor*mH),2)*RadCur_Cur * ( 1e-7 / 1e-6);
		  }
	  }  
  }    
  // from erg S^{-1} cm^{-3} to J S^{-1} cm^{-3}: *(1e-7/1e-6)
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
    CFreal wave_pressure_approx = P0_W*(x_dimless*x_dimless + y_dimless*y_dimless)/pow(r2_dimless,3);
//    CFreal wave_pressure_approx = 0.0; //P0_W* 1.0/r_dimless/r_dimless;

  

      // Gradients for p = 1/r^2
//    CFreal dP_dx = -2.0*P0_W*x_dimless*pow(r_dimless,-2.0);
//    CFreal dP_dy = -2.0*P0_W*y_dimless*pow(r_dimless,-2.0);
//    CFreal dP_dz = -2.0*P0_W*z_dimless*pow(r_dimless,-2.0);

    
    // Gradients for P = sin^2(theta)/ r^2
   // CFreal dP_dx = 2.0*P0_W*x_dimless*(pow(z_dimless,2.0)-pow(y_dimless,2.0)-pow(x_dimless,2))/pow(r_dimless,6.0);
   // CFreal dP_dy = 2.0*P0_W*y_dimless*(pow(z_dimless,2.0)-pow(y_dimless,2.0)-pow(x_dimless,2))/pow(r_dimless,6.0); 
   // CFreal dP_dz = -4.0*P0_W*z_dimless*(pow(x_dimless,2.0)+pow(y_dimless,2.0))/pow(r_dimless,6);
    CFreal dP_dx = 2.0*P0_W*x_dimless*(pow(z_dimless,2.0)-4.0*pow(y_dimless,2.0)-4.0*pow(x_dimless,2))/pow(r_dimless,8);
    CFreal dP_dy = 2.0*P0_W*y_dimless*(pow(z_dimless,2.0)-4.0*pow(y_dimless,2.0)-4.0*pow(x_dimless,2))/pow(r_dimless,8);
    CFreal dP_dz = -6.0*P0_W*z_dimless*(pow(x_dimless,2.0)+pow(y_dimless,2.0))/pow(r_dimless,8);

    CFreal lambd = 4.4/std::sqrt((*currState)[4]*(*currState)[4] + (*currState)[5]*(*currState)[5] + (*currState)[6]*(*currState)[6]);
    CFreal heating_alfwen = wave_pressure_approx*std::sqrt(8.0/(*currState)[0]*wave_pressure_approx);


    // Calcualting the gradient of variable pressure 

    const CFuint PID = 7;
    const CFuint gradPID = elementID*nbEqs + PID;
    _gradP[XX] = this->m_ux[gradPID];
    _gradP[YY] = this->m_uy[gradPID];
    _gradP[XX] = this->m_uz[gradPID];
    CFreal n_x = normals[startID];
    CFreal n_y = normals[startID+1];
    CFreal n_z = normals[startID+2];
    if (nbEqs == 11) {
    source[9]=0;
    source[10]=0;
    // Calculating the source terms for Eps equations
    const CFreal pi = 3.141592653589793238462643383279502884197;
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
    CFreal dvaxdx = 1/std::sqrt(mu0*density)*(_gradBx[XX]*B0 - 0.5*Bx/density*_gradRho[XX]*rhoRef)/RSun;
    CFreal dvaydy = 1/std::sqrt(mu0*density)*(_gradBy[YY]*B0 - 0.5*By/density*_gradRho[YY]*rhoRef)/RSun;
    CFreal dvazdz = 1/std::sqrt(mu0*density)*(_gradBz[ZZ]*B0 - 0.5*Bz/density*_gradRho[ZZ]*rhoRef)/RSun;
    CFreal dvxdx = _gradVx[XX]*vRef/RSun;
    CFreal dvydy = _gradVy[YY]*vRef/RSun;
    CFreal dvzdz = _gradVz[ZZ]*vRef/RSun;
    CFreal vxdim = Vx*vRef;
    CFreal vydim = Vy*vRef;
    CFreal vzdim = Vz*vRef;
    CFreal vax = Bx/std::sqrt(mu0*density);
    CFreal vay = By/std::sqrt(mu0*density);
    CFreal vaz = Bz/std::sqrt(mu0*density);
    CFreal R1_p = 1/(4.*density) *((vxdim - vax)*_gradRho[XX] + (vydim - vay)*_gradRho[YY]+ (vzdim - vaz)*_gradRho[ZZ])*rhoRef/RSun;
    CFreal R1_m = 1/(4.*density) *((vxdim + vax)*_gradRho[XX] + (vydim + vay)*_gradRho[YY]+ (vzdim + vaz)*_gradRho[ZZ])*rhoRef/RSun;
  //  std::cout << "R1_p  " << R1_p << "  R1_m  " << R1_m << endl;
    CFreal second_term_plus = dvxdx + dvydy + dvzdz + dvaxdx + dvaydy + dvazdz;
    CFreal second_term_minus = dvxdx + dvydy + dvzdz - dvaxdx - dvaydy - dvazdz;
    CFreal source_plus = (*currState)[9] * (R1_p + second_term_plus);
    CFreal source_minus = (*currState)[10] * (R1_m + second_term_minus);
//   std::cout << "source 9 before " << source[9] << " source 10 before " << source[10] << endl; 
//    source[9] += source_plus*RSun/vRef*volumes[elementID];
//    source[10] += source_minus*RSun/vRef*volumes[elementID];
    source[9] = 0.0;
    source[10] = 0.0;
    }
      	    // VERSION #3   
      
      // ---- DENSITY ---------------------------------------------------------
      source[0] = 0.0;
      // --- M O M E N T U M   D E N S I T Y ----------------------------------
          source[1] = 0.0;
          source[2] = 0.0;
          source[3] = 0.0;
		  //----------------------->>Artificial added continuty equation term-----------------------------------------------
		  /*
		  if (r_dimless < 1.1){
			  CFreal p_temp = (*currState)[7];
			  CFreal pth = p_temp*pRef;
			  CFreal pmag = Bnorm*Bnorm / 2. / mu0;
			  pmag = max(pmag, 1.e-16);
			  CFreal plasmaBeta = pth / pmag;
			  CFreal fac = 10.0;
			  //CFreal Q_Betafactor = 0.5*(1.0 + std::tanh((1.0 / plasmaBeta / 10.0 - 1.0)*fac));
			  CFreal Q_Betafactor = 0.5*(1.0 + std::tanh((1.05-r_dimless)/0.02));

			  CFreal rhoBD = 2.0;
			  CFreal fac_Va = 0.05;
			  CFreal vA = Bnorm / pow((density * mu0), 0.5) / vRef;
			  CFreal rho0 = max(rhoBD*std::pow(1.0 / r_dimless, 5), rhoBD*std::pow(1.0 / 22.0, 5));
			  source[0] += -fac_Va*vA*((*currState)[0] - rho0)*Q_Betafactor*volumes[elementID];
		  
			  CFreal rhovv0 = rho0*1000.0 / 480363.085276 * (-52.1 + 106.0 * std::log(r_dimless + 0.78));
			  source[1] += -fac_Va*vA*((*currState)[0] * Vx - rhovv0*x_dimless / r_dimless)*Q_Betafactor*volumes[elementID];
			  source[2] += -fac_Va*vA*((*currState)[0] * Vy - rhovv0*y_dimless / r_dimless)*Q_Betafactor*volumes[elementID];
			  source[3] += -fac_Va*vA*((*currState)[0] * Vz - rhovv0*z_dimless / r_dimless)*Q_Betafactor*volumes[elementID];
		  }
		  */
		  //-------------------------Artificial added continuty equation term<<---------------------------------------------


      if (_gravity == 1){
	source[1] += (*currState)[0]*gx*volumes[elementID];
	source[2] += (*currState)[0]*gy*volumes[elementID];
	source[3] += (*currState)[0]*gz*volumes[elementID];
	//	std::cout << "source momentum 1 " << source[1] << " density " << (*currState)[0] << " addition  " <<  (*currState)[0]*gx*volumes[elementID] << endl;
	socket_gravity.getDataHandle()[elementID*4]   = (*currState)[0]*gx; //*volumes[elementID];
	socket_gravity.getDataHandle()[elementID*4+1] = (*currState)[0]*gy; //*volumes[elementID];
	socket_gravity.getDataHandle()[elementID*4+2] = (*currState)[0]*gz; //*volumes[elementID];
      } 
      
      socket_wavepressure.getDataHandle()[elementID] =0.0;
      if (_alfven_pressure == 1) {
        source[1] -= dP_dx*volumes[elementID];
        source[2] -= dP_dy*volumes[elementID];
        source[3] -= dP_dz*volumes[elementID]; 
        socket_wavepressure.getDataHandle()[elementID] =wave_pressure_approx;
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
        source[7] += (result_flux*nondimconst)*volumes[elementID];
	
        socket_heating.getDataHandle()[elementID]   = result_flux;
      }

      if (_RadiativeLossTerm == 1) {

        
        source[7] -= Q_rad*nondimconst*volumes[elementID];
	
        socket_radiativeloss.getDataHandle()[elementID] = Q_rad;
      }

	  //>> energy Source terms caused by decomposed energy equation
	  //CFuint deCompE = 1;
	  CFreal Q_deCompE = 0.0;
	  //if (deCompE == 1 && r_dimless>1.05){
	  //_deCompE = 1;
	  if (_deCompE == 1){
		  const CFuint BXID = 4;
		  const CFuint BYID = 5;
		  const CFuint BZID = 6;
		  const CFuint gradBXID = elementID*nbEqs + BXID;
		  const CFuint gradBYID = elementID*nbEqs + BYID;
		  const CFuint gradBZID = elementID*nbEqs + BZID;
		  _gradBx[XX] = this->m_ux[gradBXID];
		  _gradBx[YY] = this->m_uy[gradBXID];
		  _gradBx[ZZ] = this->m_uz[gradBXID];
		  _gradBy[XX] = this->m_ux[gradBYID];
		  _gradBy[YY] = this->m_uy[gradBYID];
		  _gradBy[ZZ] = this->m_uz[gradBYID];
		  _gradBz[XX] = this->m_ux[gradBZID];
		  _gradBz[YY] = this->m_uy[gradBZID];
		  _gradBz[ZZ] = this->m_uz[gradBZID];
		  std::vector<CFreal> BdotLamdaB(3, 0.0);
		  std::vector<CFreal> VdotLamdaB(3, 0.0);
		  VdotLamdaB[0] = Vx*_gradBx[XX] + Vy*_gradBx[YY] + Vz*_gradBx[ZZ];
		  VdotLamdaB[1] = Vx*_gradBy[XX] + Vy*_gradBy[YY] + Vz*_gradBy[ZZ];
		  VdotLamdaB[2] = Vx*_gradBz[XX] + Vy*_gradBz[YY] + Vz*_gradBz[ZZ];
		  BdotLamdaB[0] = (Bx*_gradBx[XX] + By*_gradBx[YY] + Bz*_gradBx[ZZ]) / B0;
		  BdotLamdaB[1] = (Bx*_gradBy[XX] + By*_gradBy[YY] + Bz*_gradBy[ZZ]) / B0;
		  BdotLamdaB[2] = (Bx*_gradBz[XX] + By*_gradBz[YY] + Bz*_gradBz[ZZ]) / B0;
		  Q_deCompE = Vx*BdotLamdaB[0] + Vy*BdotLamdaB[1] + Vz*BdotLamdaB[2] -
			  (Bx*VdotLamdaB[0] + By*VdotLamdaB[1] + Bz*VdotLamdaB[2]) / B0;
		  source[7] += Q_deCompE*volumes[elementID];
	  }
	  //<< energy Source terms caused by decomposed energy equation

      source[8] = 0.0;
      // Eps modification 
  //    if (states->size() == 11) {
  //    source[9] = source_plus/vRef*volumes[elementID];
  //    source[10] = source_minus/vRef*volumes[elementID];
  //    socket_zp.getDataHandle()[elementID] = source_plus/vRef;
  //    socket_zm.getDataHandle()[elementID] = source_minus/vRef;
  //    }
  CFLog(DEBUG_MAX, "MHDConsACAHWSourceTerm::computeSource() => END\n");
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<Framework::BaseDataSocketSource> >
MHDConsACAHWSourceTerm::providesSockets()
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
