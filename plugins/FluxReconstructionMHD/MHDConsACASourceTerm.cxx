#include "FluxReconstructionMHD/MHDConsACASourceTerm.hh"
#include "Common/CFLog.hh"
#include "Framework/GeometricEntity.hh"
#include "Framework/MeshData.hh"
#include "FluxReconstructionMHD/FluxReconstructionMHD.hh"
#include "Framework/SubSystemStatus.hh"

#include "Framework/MethodCommandProvider.hh"
#include "MathTools/MathConsts.hh"

#include "MHD/MHDTerm.hh"

#include "FluxReconstructionMethod/FluxReconstructionElementData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::MHD;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace FluxReconstructionMethod {

////////////////////////////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<MHDConsACASourceTerm, FluxReconstructionSolverData, FluxReconstructionMHDModule>
MHDConsACASourceTermFRProvider("MHDConsACASourceTerm");

///////////////////////////////////////////////////////////////////////////////////////////////////

void MHDConsACASourceTerm::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< bool >("Gravity","Switch on gravity.");
  options.addConfigOption< bool >("PevtsovHeating","Switch on PevtsovHeating term.");
  options.addConfigOption< CFreal >("PevtsovHeatingFactor","Scale PevtsovHeating term.");
  options.addConfigOption< bool >("Manchester","Switch on Manchester heating term.");
  options.addConfigOption< CFreal >("ManchesterHeatingAmplitude","Scale Manchester volumetric heating amplitude.");
  options.addConfigOption< CFreal >("ManchesterSigma","Sigma value of Manchester heating term.");
  options.addConfigOption< bool >("DivQ","Switch on heat conduction.");
  options.addConfigOption< CFreal >("DivQConductivity","Conductivity of heat conduction.");
  options.addConfigOption< CFreal >("DivQalphaCollisionless","alpha value of heat conductivity in collisionless regime.");
  options.addConfigOption< bool >("ViscosityAndResistivity","Switch on viscosity and resisitivity.");
  options.addConfigOption< CFreal >("Viscosity","nu value.");
  options.addConfigOption< CFreal >("Resistivity","eta value.");
  options.addConfigOption< bool >("RadiativeLossTerm","Switch on optically thin approximation for radiation losses.");
  options.addConfigOption< bool >("AddUpdateCoeff","Add the ST time step restriction.");
}

/////////////////////////////////////////////////////////////////////////////////////////////////////

MHDConsACASourceTerm::MHDConsACASourceTerm(const std::string& name) :
  StdSourceTerm(name),
  socket_updateCoeff("updateCoeff"),
  socket_gravity("gravity"),
  socket_Br("Br"),
  socket_Btheta("Btheta"),
  socket_Bphi("Bphi"),
  socket_Vr("Vr"),
  socket_Vtheta("Vtheta"),
  socket_Vphi("Vphi"),
  m_order()
{ 
  addConfigOptionsTo(this);
  
  m_gravity = false; 
  setParameter("Gravity",&m_gravity);
  
  m_PevtsovHeating = false;
  setParameter("PevtsovHeating",&m_PevtsovHeating);
  
  m_PevtsovHeatingFactor = 1.0;
  setParameter("PevtsovHeatingFactor",&m_PevtsovHeatingFactor);

  m_Manchester = false;
  setParameter("Manchester",&m_Manchester);
  
  m_ManchesterHeatingAmplitude = 0.0;
  setParameter("ManchesterHeatingAmplitude",&m_ManchesterHeatingAmplitude);
  
  m_ManchesterSigma = 100.0;
  setParameter("ManchesterSigma",&m_ManchesterSigma);
  
  m_divQ = false;
  setParameter("DivQ",&m_divQ);
  
  m_divQConductivity = 0.0;
  setParameter("DivQConductivity",&m_divQConductivity);
  
  m_divQalphaCollisionless = 0.0;
  setParameter("DivQalphaCollisionless",&m_divQalphaCollisionless);
  
  m_ViscosityAndResistivity = false;
  setParameter("ViscosityAndResistivity",&m_ViscosityAndResistivity);
  
  m_Viscosity = 0.0;
  setParameter("Viscosity",&m_Viscosity);
  
  m_Resistivity = 0.0;
  setParameter("Resistivity",&m_Resistivity);
  
  m_RadiativeLossTerm = false;
  setParameter("RadiativeLossTerm",&m_RadiativeLossTerm);
  
  m_addUpdateCoeff = false;
  setParameter("AddUpdateCoeff",&m_addUpdateCoeff);
}

/////////////////////////////////////////////////////////////////////////////////////////////////

MHDConsACASourceTerm::~MHDConsACASourceTerm()
{
}

//////////////////////////////////////////////////////////////////////////////////////////////////

void MHDConsACASourceTerm::setup()
{
  StdSourceTerm::setup();
  
  DataHandle< CFreal > gravity = socket_gravity.getDataHandle();
  
  DataHandle< CFreal > updateCoeff = socket_updateCoeff.getDataHandle();
  
  // resize socket
  gravity.resize(updateCoeff.size()*4);
  
  // get the local FR data
  vector< FluxReconstructionElementData* >& frLocalData = getMethodData().getFRLocalData();
  cf_assert(frLocalData.size() > 0);
  // for now, there should be only one type of element
  cf_assert(frLocalData.size() == 1);
  
  const CFPolyOrder::Type order = frLocalData[0]->getPolyOrder();
  
  m_order = static_cast<CFuint>(order);
  
  DataHandle< CFreal > Br = socket_Br.getDataHandle();
  DataHandle< CFreal > Btheta = socket_Btheta.getDataHandle();
  DataHandle< CFreal > Bphi = socket_Bphi.getDataHandle();
  
  DataHandle< CFreal > Vr = socket_Vr.getDataHandle();
  DataHandle< CFreal > Vtheta = socket_Vtheta.getDataHandle();
  DataHandle< CFreal > Vphi = socket_Vphi.getDataHandle();
  
  const CFuint nbStates = updateCoeff.size();

  // resize socket
  Br.resize(nbStates);
  Btheta.resize(nbStates);
  Bphi.resize(nbStates);
  Vr.resize(nbStates);
  Vtheta.resize(nbStates);
  Vphi.resize(nbStates);
}

//////////////////////////////////////////////////////////////////////////////////////////////////

void MHDConsACASourceTerm::unsetup()
{
  StdSourceTerm::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

std::vector< Common::SafePtr< BaseDataSocketSource > >
  MHDConsACASourceTerm::providesSockets()
{
  std::vector< Common::SafePtr< BaseDataSocketSource > > result = StdSourceTerm::providesSockets();
  result.push_back(&socket_gravity);
  result.push_back(&socket_Br);
  result.push_back(&socket_Btheta);
  result.push_back(&socket_Bphi);
  result.push_back(&socket_Vr);
  result.push_back(&socket_Vtheta);
  result.push_back(&socket_Vphi);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

std::vector< Common::SafePtr< BaseDataSocketSink > >
MHDConsACASourceTerm::needsSockets()
{
  std::vector< Common::SafePtr< BaseDataSocketSink > > result = StdSourceTerm::needsSockets();
  result.push_back(&socket_updateCoeff);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

void MHDConsACASourceTerm::addSourceTerm(RealVector& resUpdates)
{       
  CFLog(VERBOSE, "MHDConsACASourceTerm::addSourceTerm() => START\n");
  
  // set gradients
  const CFuint nbrStates = m_cellStates->size();
   
  const CFreal RSun = 6.9551e8; // m
  const CFreal rho_ref = 1.67e-13;
  const CFreal GMsun = 1.327474512e20; // SI value
  const CFreal gsun = 274.00; // at the surface m per sec2
  const CFreal force_density_ref = std::pow(2.2e-4,2)/(1.25e-6);//*6.9e8);
  const CFreal th1 = (MathTools::MathConsts::CFrealPi()/180.0)*17.5;
  const CFreal th2 = (MathTools::MathConsts::CFrealPi()/180.0)*61.5;
  
  // volumetric heating amplitude (J/(kg*s*K))
  const CFreal q0 = m_ManchesterHeatingAmplitude;

  const CFreal lRef = 6.9e8;
  
  // Boltzmann constant
  const CFreal k = 1.3806503e-23;

  // mass of proton
  const CFreal mp = 1.67262158e-27;

  // mass of electron
  const CFreal me = 9.10938188e-31;

  //const CFreal rhoRef = nRef*(mp+me);
  const CFreal rhoRef = 1.67e-13;

  const CFreal vRef = 480.248;
  
  const CFreal mu_cor = 1.27; // Molecular weight in the corona
  const CFreal mH = 1.6733e-27; // kg
  
  // VERSION #3
  const CFreal l0 = 6.9e8;
  const CFreal B0 = 2.2e-4;
  const CFreal mu0 = 1.2566e-6;
  const CFreal rho0 = 1.67e-13;
  const CFreal V0 = B0/std::sqrt(mu0*rho0);
  const CFreal g0 = pow(V0,2)/l0;
  
  for (CFuint iSol = 0; iSol < nbrStates; ++iSol)
  {      
    SafePtr<MHDTerm> model = PhysicalModelStack::getActive()->getImplementor()->getConvectiveTerm().d_castTo< MHDTerm >();

    // get state ID
    const CFuint stateID = (((*m_cellStates)[iSol]))->getLocalID();
  
    const CFreal refSpeed = model->getRefSpeed();
    const CFreal refSpeedSq = refSpeed*refSpeed;

    const CFuint dim = PhysicalModelStack::getActive()->getDim(); 

    const std::string correctionType = model->getCorrectionType();
      
//    for (CFuint i = 0; i < (m_nbrEqs-1); ++i) 
//    {
//      resUpdates[m_nbrEqs*iSol + i] = 0.0;
//    }

    if (correctionType == "Mixed") 
    {
      // mixed (hyperbolic and parabolic) correction
      const CFreal dissipCoeff = model->getDissipationCoefficient();
      const CFreal dissipCoeffSq = dissipCoeff*dissipCoeff;

      resUpdates[m_nbrEqs*iSol + m_nbrEqs-1] = -(refSpeedSq/dissipCoeffSq)*(*((*m_cellStates)[iSol]))[m_nbrEqs-1];//*volumes[elementID];
    }
    else 
    {
      // hyperbolic correction
      resUpdates[m_nbrEqs*iSol + m_nbrEqs-1] = 0.0;
    }

    const CFreal x_dimless = ((*m_cellStates)[iSol]->getCoordinates())[XX];
    const CFreal y_dimless = ((*m_cellStates)[iSol]->getCoordinates())[YY];
    const CFreal z_dimless = (dim==DIM_3D) ? ((*m_cellStates)[iSol]->getCoordinates())[ZZ] : 0.0;
    const CFreal x = x_dimless*RSun;
    const CFreal y = y_dimless*RSun;
    const CFreal z = z_dimless*RSun;
      
    const CFreal r2_dimless = x_dimless*x_dimless + y_dimless*y_dimless + z_dimless*z_dimless;
    const CFreal r_dimless = std::sqrt(r2_dimless);

    const CFreal r2 = x*x + y*y + z*z;
    const CFreal r = std::sqrt(r2);
    const CFreal density = (*((*m_cellStates)[iSol]))[0]*rho_ref;
    const CFreal density_dimless = (*((*m_cellStates)[iSol]))[0];
    const CFreal fgrav = -density*gsun*pow(6.9551e8,2)/r2;

    // gravity term in dimensionless units:
    const CFreal fgravity_dimless = -0.82977*density_dimless/r2_dimless;  
    // the numerical factor is - G*MSun*mu_0*rho_0/(B_0**2*R_Sun)

    const CFreal fx = fgrav*x/r;
    const CFreal fy = fgrav*y/r;
    const CFreal fz = fgrav*z/r;

    const CFreal Vr_dimless = x_dimless/r_dimless*(*((*m_cellStates)[iSol]))[1] + y_dimless/r_dimless*(*((*m_cellStates)[iSol]))[2] + z_dimless/r_dimless*(*((*m_cellStates)[iSol]))[3];

    // B-FIELD WEIGHTED 3-D HEATING MODEL - Petsov et al. 2003
    // First attempt: Instead of summing /int /psi dV over all tetrahedrons, which has to be done prior to calling this function
    // we can simply use the constant factor provided by Downs et al. 2010 of 4e-8 J/s T
    // Then the heating term is simply
    //const CFreal const_Downs2010 = 4.e-8/0.03851; // to adimensional value    // Originally 4.e-8
    //const CFreal Qh2 = const_Downs2010*std::sqrt((*currState)[4]*(*currState)[4] + (*currState)[5]*(*currState)[5] + (*currState)[6]*(*currState)[6]);

    const CFreal theta = std::atan2(std::sqrt(x_dimless*x_dimless + y_dimless*y_dimless),z_dimless);

    // critical angle
    CFreal critTheta = 0.0;
    if (r_dimless <= 7.0) 
    {
      critTheta = asin(sqrt(sin(th1)*sin(th1)+cos(th1)*cos(th1)*(r-1.0)/8.0));
    }
    else if ((r_dimless > 7.0) && (r_dimless < 30.0)) 
    {
      critTheta = asin(sqrt(sin(th2)*sin(th2)+cos(th2)*cos(th2)*(r-7.0)/40.0));
    }
    
    // target temperature and the heating scale height
    CFreal T0 = 0.0;
    CFreal sigma0 = 0.0;
    if (fabs(0.5*MathTools::MathConsts::CFrealPi()-theta) < critTheta) 
    {
      T0 = 1.5e6;
      sigma0 = 4.5*lRef;
    }
    else if (fabs(0.5*MathTools::MathConsts::CFrealPi()-theta) > critTheta) 
    {
      T0 = 2.63e6;
      sigma0 = 4.5*lRef*(2.0-(sin(theta)*sin(theta)/(sin(critTheta)*sin(critTheta))));
    }

    const CFreal T = (*((*m_cellStates)[iSol]))[7]*0.03851*1.27*mp*vRef*vRef/(density_dimless*1.67e-13*k);
    const CFreal T2 = (*((*m_cellStates)[iSol]))[7]*0.03851*1.27*mp*vRef*vRef/(2.0*density_dimless*1.67e-13*k);

    // dimensional volumetric heating term (J/(m^3*s))
    const CFreal Q = density_dimless*1.67e-13*q0*(T0-T2)*exp(-(r_dimless-1.0)*(r_dimless-1.0)/(sigma0*sigma0));

    const CFreal nondimconst = lRef/(rhoRef*pow(vRef,3));

    // Radiative loss pleitner Dec 4
    // Q_rad = n_e^2*10^(-17.73)*T^(-2./3.) = (rho/(mu*mH))^2*10^(-17.73)*T^(-2./3.)
    const CFreal Q_rad = pow(density_dimless*1.67e-13/(mu_cor*mH),2)*pow(10,-17.73)*pow(T2,-2.0/3.0);

    const CFreal g = -(GMsun/r2)/g0;
    const CFreal gx = g*x_dimless/r_dimless;
    const CFreal gy = g*y_dimless/r_dimless;
    const CFreal gz = g*z_dimless/r_dimless;  
      
      
    // density
    resUpdates[m_nbrEqs*iSol + 0] = 0.0;

    // V
    resUpdates[m_nbrEqs*iSol + 1] = 0.0;
    resUpdates[m_nbrEqs*iSol + 2] = 0.0;
    resUpdates[m_nbrEqs*iSol + 3] = 0.0;

    if (m_gravity)
    {
      resUpdates[m_nbrEqs*iSol + 1] += density_dimless*gx;//*volumes[elementID];
      resUpdates[m_nbrEqs*iSol + 2] += density_dimless*gy;//*volumes[elementID];	
      resUpdates[m_nbrEqs*iSol + 3] += density_dimless*gz;//*volumes[elementID];	
    }
      
    // B
    resUpdates[m_nbrEqs*iSol + 4] = 0.0;
    resUpdates[m_nbrEqs*iSol + 5] = 0.0;
    resUpdates[m_nbrEqs*iSol + 6] = 0.0;
      
    // T
    const CFreal Vx = (*((*m_cellStates)[iSol]))[1];
    const CFreal Vy = (*((*m_cellStates)[iSol]))[2];
    const CFreal Vz = (*((*m_cellStates)[iSol]))[3];
    const CFreal Vdotg = Vx*gx + Vy*gy + Vz*gz; // dimensionless
      
    resUpdates[m_nbrEqs*iSol + 7] = 0.0;
    
    if (m_gravity)
    {
      resUpdates[m_nbrEqs*iSol + 7] += density_dimless*Vdotg;//*volumes[elementID];
    }
      
    if (m_Manchester) 
    {
      resUpdates[m_nbrEqs*iSol + 7] += Q*nondimconst;//*volumes[elementID];
    }

    if (m_RadiativeLossTerm) 
    {
      resUpdates[m_nbrEqs*iSol + 7] -= Q_rad*nondimconst;//*volumes[elementID];
    }

    // phi
    resUpdates[m_nbrEqs*iSol + 8] = 0.0;
  
  
    // Set gravity socket
    if (!m_isPerturbed)
    {
      if (m_gravity)
      {
        socket_gravity.getDataHandle()[stateID*4]   = density_dimless*gx; //*volumes[elementID];	
        socket_gravity.getDataHandle()[stateID*4+1] = density_dimless*gy; //*volumes[elementID];	
        socket_gravity.getDataHandle()[stateID*4+2] = density_dimless*gz; //*volumes[elementID];
        socket_gravity.getDataHandle()[stateID*4+3] = density_dimless*Vdotg; //*volumes[elementID];
      }
      
      DataHandle< CFreal > Br = socket_Br.getDataHandle();
      DataHandle< CFreal > Btheta = socket_Btheta.getDataHandle();
      DataHandle< CFreal > Bphi = socket_Bphi.getDataHandle();
      
      DataHandle< CFreal > Vr = socket_Vr.getDataHandle();
      DataHandle< CFreal > Vtheta = socket_Vtheta.getDataHandle();
      DataHandle< CFreal > Vphi = socket_Vphi.getDataHandle();
      
      const CFreal r2D = sqrt(x_dimless*x_dimless+y_dimless*y_dimless);
      
      const CFreal Bxv = (*((*m_cellStates)[iSol]))[4];
      const CFreal Byv = (*((*m_cellStates)[iSol]))[5];
      const CFreal Bzv = (*((*m_cellStates)[iSol]))[6];
      
      Br[stateID] = x_dimless/r_dimless*Bxv + y_dimless/r_dimless*Byv + z_dimless/r_dimless*Bzv;
    
      Bphi[stateID] = (-y_dimless*Bxv + x_dimless*Byv)/r2D;
    
      Btheta[stateID] = (z_dimless*x_dimless/r2D*Bxv + z_dimless*y_dimless/r2D*Byv - r2D*Bzv)/r_dimless;
      
      Vr[stateID] = x_dimless/r_dimless*Vx + y_dimless/r_dimless*Vy + z_dimless/r_dimless*Vz;
    
      Vphi[stateID] = (-y_dimless*Vx + x_dimless*Vy)/r2D;
    
      Vtheta[stateID] = (z_dimless*x_dimless/r2D*Vx + z_dimless*y_dimless/r2D*Vy - r2D*Vz)/r_dimless;
    }
  }
  
  CFLog(VERBOSE, "MHDConsACASourceTerm::addSourceTerm() => END\n");
}

//////////////////////////////////////////////////////////////////////////////////////////////////

void  MHDConsACASourceTerm::getSToStateJacobian(const CFuint iState)
{
  // reset the jacobian
  for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
  {
    m_stateJacobian[iEq] = 0.0;
  }
  
  if (m_gravity)
  {
  const CFuint dim = PhysicalModelStack::getActive()->getDim(); 
  
  const CFreal RSun = 6.9551e8; // m
  const CFreal l0 = 6.9e8;
  const CFreal B0 = 2.2e-4;
  const CFreal mu0 = 1.2566e-6;
  const CFreal rho0 = 1.67e-13;
  const CFreal V0 = B0/std::sqrt(mu0*rho0);
  const CFreal g0 = pow(V0,2)/l0;
  
  const CFreal x_dimless = ((*m_cellStates)[iState]->getCoordinates())[XX];
  const CFreal y_dimless = ((*m_cellStates)[iState]->getCoordinates())[YY];
  const CFreal z_dimless = (dim==DIM_3D) ? ((*m_cellStates)[iState]->getCoordinates())[ZZ] : 0.0;
  const CFreal x = x_dimless*RSun;
  const CFreal y = y_dimless*RSun;
  const CFreal z = z_dimless*RSun;
  
  const CFreal r2_dimless = x_dimless*x_dimless + y_dimless*y_dimless + z_dimless*z_dimless;
  const CFreal r_dimless = std::sqrt(r2_dimless);
  
  const CFreal r2 = x*x + y*y + z*z;
  
  const CFreal GMsun = 1.327474512e20; // SI value
  
  const CFreal g = -(GMsun/r2)/g0;
  
  const CFreal gx = g*x_dimless/r_dimless;
  const CFreal gy = g*y_dimless/r_dimless;
  const CFreal gz = g*z_dimless/r_dimless;  
  
  /////// Gravity
  
  // V
  m_stateJacobian[0][1] = gx;
  m_stateJacobian[0][2] = gy;
  m_stateJacobian[0][3] = gz;
  
  // p
  const CFreal density_dimless = (*((*m_cellStates)[iState]))[0];
  const CFreal Vx = (*((*m_cellStates)[iState]))[1];
  const CFreal Vy = (*((*m_cellStates)[iState]))[2];
  const CFreal Vz = (*((*m_cellStates)[iState]))[3];
  const CFreal Vdotg = Vx*gx + Vy*gy + Vz*gz; // dimensionless
    
  m_stateJacobian[0][7] = Vdotg;
  m_stateJacobian[1][7] = density_dimless*gx;
  m_stateJacobian[2][7] = density_dimless*gy;
  m_stateJacobian[3][7] = density_dimless*gz;
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluxReconstructionMethod

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
