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
  
  const CFuint nbStates = updateCoeff.size();

  // resize socket
  Br.resize(nbStates);
  Btheta.resize(nbStates);
  Bphi.resize(nbStates);
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
      
      const CFreal r2D = sqrt(x_dimless*x_dimless+y_dimless*y_dimless);
      
      const CFreal Bxv = (*((*m_cellStates)[iSol]))[4];
      const CFreal Byv = (*((*m_cellStates)[iSol]))[5];
      const CFreal Bzv = (*((*m_cellStates)[iSol]))[6];
      
      Br[stateID] = x_dimless/r_dimless*Bxv + y_dimless/r_dimless*Byv + z_dimless/r_dimless*Bzv;
    
      Btheta[stateID] = (-y_dimless*Bxv + x_dimless*Byv)/r2D;
    
      Bphi[stateID] = (z_dimless*x_dimless/r2D*Bxv + z_dimless*y_dimless/r2D*Byv - r2D*Bzv)/r_dimless;
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
  
//  SafePtr< NavierStokes3DKLogOmega > navierStokesVarSet = m_diffVarSet.d_castTo< NavierStokes3DKLogOmega >();
//  
//  // Set the wall distance before computing the turbulent viscosity
//  navierStokesVarSet->setWallDistance(m_currWallDist[iState]);
//    
//  /// destruction term of k
//    
//  const CFreal avK     = max((*((*m_cellStates)[iState]))[5],0.0);
//  const CFreal avOmega = exp((*((*m_cellStates)[iState]))[6]);
//  const CFreal T = max(1.0e-3,(*((*m_cellStates)[iState]))[4]);
//  const CFreal p = max(1.0e-3,(*((*m_cellStates)[iState]))[0]);
//  const CFreal rho = max(navierStokesVarSet->getDensity(*((*m_cellStates)[iState])),1.0e-8);
//  const CFreal R = m_eulerVarSet->getModel()->getR();
//  const CFreal overRT = 1.0/(R*T);
//  const CFreal pOverRTT = p/(R*T*T);
//  const CFreal overOmega = 1.0/avOmega;
//  const CFreal avV = sqrt((*((*m_cellStates)[iState]))[1]*(*((*m_cellStates)[iState]))[1]+(*((*m_cellStates)[iState]))[2]*(*((*m_cellStates)[iState]))[2]+(*((*m_cellStates)[iState]))[3]*(*((*m_cellStates)[iState]))[3]);
//  
//  const CFreal avGa    = min(max((*((*m_cellStates)[iState]))[7],0.01),0.99);
//    
//  const CFreal avAlpha    = max((*((*m_cellStates)[iState]))[8],0.0);
//  
//  const CFreal mu = navierStokesVarSet->getLaminarDynViscosityFromGradientVars(*((*m_cellStates)[iState]));
//  const CFreal mut = navierStokesVarSet->getTurbDynViscosityFromGradientVars(*((*m_cellStates)[iState]), m_cellGrads[iState]);
//  
//  //Compute Strain 
//  getStrain(iState);
//  
//  // Get Vorticity
//  getVorticity(iState);
//  
//  const CFreal u = (*((*m_cellStates)[iState]))[1];
//  const CFreal v = (*((*m_cellStates)[iState]))[2];
//  const CFreal w = (*((*m_cellStates)[iState]))[3];
//  const CFreal dpx = (*(m_cellGrads[iState][0]))[XX];
//  const CFreal dpy = (*(m_cellGrads[iState][0]))[YY];
//  const CFreal dpz = (*(m_cellGrads[iState][0]))[ZZ];
//    
//  // compute local freestream values with isentropic assumption
//  const CFreal MInf = m_eulerVarSet->getModel()->getMachInf();
//  const CFreal pInf = m_eulerVarSet->getModel()->getStaticPressInf();//getPressInf();
//  const CFreal gammaIsentropic = m_eulerVarSet->getModel()->getGamma();
//  const CFreal uInf = m_eulerVarSet->getModel()->getVelInf();
//  const CFreal rhoInf = (MInf*MInf)/(uInf*uInf)*gammaIsentropic*pInf;
//  
//  const CFreal pTotTerm = 1.0+(gammaIsentropic-1.0)/2.0*MInf*MInf;
//  const CFreal pTotExponent = gammaIsentropic/(gammaIsentropic-1.0);
//  const CFreal pTotalInf = pow(pTotTerm,pTotExponent)*p;
//      
//  const CFreal rhoInfLocal = rhoInf*pow(p/pInf,1.0/gammaIsentropic);
//  const CFreal MInfLocal = sqrt(2.0/(gammaIsentropic-1.0)*(pow(pTotalInf/p,1.0/pTotExponent)-1.0));
//  const CFreal uInfLocal = sqrt(gammaIsentropic*p/rhoInfLocal)*MInfLocal;
//  const CFreal TInfLocal = p/(rhoInfLocal*R);
//    
//  const CFreal Tu = min(max(100.0 * (std::sqrt(2.0*avK/3.0))/(avV),0.027),13.0);
//      
//  (*((*m_cellStates)[iState]))[4] = TInfLocal;
//    
//  const CFreal muInfLocal = navierStokesVarSet->getLaminarDynViscosityFromGradientVars(*((*m_cellStates)[iState]));
//    
//  (*((*m_cellStates)[iState]))[4] = T;
//      
//  CFreal ReThetat;
//
//  if (!m_PGrad)
//  {
//    ReThetat = rhoInfLocal*uInfLocal*uInfLocal*0.2247141166/(avAlpha*sqrt(rhoInfLocal*muInfLocal));
//  }
//  else
//  {
//    ReThetat = getRethetatwithPressureGradient(avAlpha,rhoInfLocal,uInfLocal,muInfLocal,mu,rho,dpx,dpy,dpz);
//  }
//
//  const CFreal TuInfLocal = getTuInfLocal(ReThetat,MInfLocal,Tu);
//  
//  const CFreal dux = (*(m_cellGrads[iState][1]))[XX];
//  const CFreal duy = (*(m_cellGrads[iState][1]))[YY]; 
//  const CFreal duz = (*(m_cellGrads[iState][1]))[ZZ]; 
//  const CFreal dvx = (*(m_cellGrads[iState][2]))[XX]; 
//  const CFreal dvy = (*(m_cellGrads[iState][2]))[YY]; 
//  const CFreal dvz = (*(m_cellGrads[iState][2]))[ZZ];
//  const CFreal dwx = (*(m_cellGrads[iState][3]))[XX]; 
//  const CFreal dwy = (*(m_cellGrads[iState][3]))[YY]; 
//  const CFreal dwz = (*(m_cellGrads[iState][3]))[ZZ];
//  const CFreal dgammax = (*(m_cellGrads[iState][7]))[XX];
//  const CFreal dgammay = (*(m_cellGrads[iState][7]))[YY];
//  const CFreal dgammaz = (*(m_cellGrads[iState][7]))[ZZ];
//  const CFreal dkx = (*(m_cellGrads[iState][5]))[XX];
//  const CFreal dky = (*(m_cellGrads[iState][5]))[YY];
//  const CFreal dkz = (*(m_cellGrads[iState][5]))[ZZ];
//  const CFreal dax = (*(m_cellGrads[iState][8]))[XX];
//  const CFreal day = (*(m_cellGrads[iState][8]))[YY];
//  const CFreal daz = (*(m_cellGrads[iState][8]))[ZZ];
//  
//  const CFreal betaStar = navierStokesVarSet->getBetaStar(*((*m_cellStates)[iState]));
//  
//  const CFreal DkTerm = avOmega * betaStar;
//  
//  const CFreal Dk = -rho * avK * DkTerm;
//  
//  const CFreal dmudg = 0.5*(1-2.0*avGa)*(tanh((avGa-0.25)/0.1)+1) + 0.5/0.1*(1-avGa)*avGa/(pow(cosh((avGa-0.25)/0.1),2)) + 1.0;
//
//  //const CFreal limGamma = max(0.01,avGa);
//  const CFreal mutGaMod = (avGa+avGa*(1.0-avGa)*0.5*(1.0+tanh((avGa-0.25)/(0.1))));
//    
//  if (Dk < 0.0)
//  {
//    CFreal tempSTTerm = 0.0;
//  
//    if (!m_blockDecoupled)
//    {
//      //p
//      tempSTTerm = -avK * overRT * DkTerm;
//      m_stateJacobian[0][5] += tempSTTerm;
//      m_stateJacobian[0][4] -= tempSTTerm;
//  
//      //T
//      tempSTTerm = DkTerm * avK * pOverRTT;
//      m_stateJacobian[4][5] += tempSTTerm;
//      m_stateJacobian[4][4] -= tempSTTerm;
//    }
//  
//    //k
//    tempSTTerm = -DkTerm * rho;
//    m_stateJacobian[5][5] += tempSTTerm;
//    m_stateJacobian[5][4] -= tempSTTerm;
//  
//    //logOmega
//    tempSTTerm = -DkTerm * avK * rho;
//    m_stateJacobian[6][5] += tempSTTerm;
//    m_stateJacobian[6][4] -= tempSTTerm;
//  }
//  
//  /// destruction term of logOmega
//  
//  const CFreal beta = navierStokesVarSet->getBeta(*((*m_cellStates)[iState]));
//  
//  const CFreal DomegaTerm = avOmega * beta;
//  
//  const CFreal Domega = -DomegaTerm * rho;
//  
//  if (Domega < 0.0)
//  {
//    if (!m_blockDecoupled)
//    {
//      //p
//      m_stateJacobian[0][6] = -DomegaTerm * overRT;
//  
//      //T
//      m_stateJacobian[4][6] = DomegaTerm * pOverRTT;
//    }
//  
//    //logOmega
//    m_stateJacobian[6][6] = -DomegaTerm * rho;
//  }
//  
//  /// production term of k
//
//  const CFreal coeffTauMu = navierStokesVarSet->getModel().getCoeffTau();
//  //const CFreal mutTerm = coeffTauMu*((4./3.)*((dux-dvy)*(dux-dvy)+(dux*dvy)-(dux+dvy-m_vOverRadius)*m_vOverRadius)+(duy+dvx)*(duy+dvx));
//  const CFreal gammaTerm = avGa + avGa*(1.0-avGa);
//  
//  CFreal mutTerm;
//  CFreal Pk;
//  
//  if (m_neglectSSTVTerm)
//  {
//    mutTerm = coeffTauMu*m_vorticity*m_vorticity;
//
//    Pk = mutTerm*mut;  
//  }
//  else
//  {
//    mutTerm = coeffTauMu*m_vorticity*m_vorticity;
//
//    Pk = (mutTerm*mut - (2./3.)*(avK * rho)*(dux+dvy+dwz)); 
//  }
//  
//  //const CFreal Pk = (mutTerm*mut - (2./3.)*(avK * rho * gammaTerm)*(dux+dvy+m_vOverRadius));
//  
//  const CFreal twoThirdduxduy = (2./3.)*(dux+dvy+dwz);
//  
//  if (Pk > 0.0)
//  {
//    CFreal tempSTTerm = 0.0;
//       
//    if (!m_blockDecoupled)
//    {
//      //p
//      tempSTTerm = mutTerm*avK*overRT*overOmega*mutGaMod;
//      if (!m_neglectSSTVTerm) tempSTTerm += -twoThirdduxduy*avK*overRT * gammaTerm;
//      
//      m_stateJacobian[0][5] += tempSTTerm;
//      m_stateJacobian[0][4] -= tempSTTerm;
//  
//      //T
//      tempSTTerm = -mutTerm*avK*overOmega*mutGaMod*pOverRTT;
//      if (!m_neglectSSTVTerm) tempSTTerm += twoThirdduxduy*avK*pOverRTT * gammaTerm;
//      
//      m_stateJacobian[4][5] += tempSTTerm;
//      m_stateJacobian[4][4] -= tempSTTerm;
// 
//      //gamma
//      tempSTTerm = mutTerm*mut*dmudg/mutGaMod;
//      if (!m_neglectSSTVTerm) tempSTTerm += -twoThirdduxduy * avK * rho*2.0*(1.0-avGa);
//      
//      m_stateJacobian[7][5] += tempSTTerm;
//      m_stateJacobian[7][4] -= tempSTTerm;
//    }
//  
//    //k
//    tempSTTerm = mutTerm*rho*overOmega*mutGaMod;
//    if (!m_neglectSSTVTerm) tempSTTerm += -twoThirdduxduy*rho * gammaTerm;
//    
//    m_stateJacobian[5][5] += tempSTTerm;
//    m_stateJacobian[5][4] -= tempSTTerm;
//  
//    //logOmega
//    tempSTTerm = -mutTerm*rho*avK*overOmega*mutGaMod;
//    m_stateJacobian[6][5] += tempSTTerm;
//    m_stateJacobian[6][4] -= tempSTTerm;
//  }
//  
//  /// production of logOmega
//  
//  const CFreal gamma = navierStokesVarSet->getGammaCoef();
//  const CFreal blendingCoefF1 = navierStokesVarSet->getBlendingCoefficientF1();
//  const CFreal sigmaOmega2 = navierStokesVarSet->getSigmaOmega2();
//  
//  const CFreal pOmegaFactor = (1. - blendingCoefF1) * 2. * sigmaOmega2 * ((*(m_cellGrads[iState][5]))[XX]*(*(m_cellGrads[iState][6]))[XX] + (*(m_cellGrads[iState][5]))[YY]*(*(m_cellGrads[iState][6]))[YY] + (*(m_cellGrads[iState][5]))[ZZ]*(*(m_cellGrads[iState][6]))[ZZ]);
//  
//  const CFreal Pomega = (gamma*rho/mut) * Pk * overOmega + rho * overOmega * pOmegaFactor;
//  
//  if (Pomega > 0.0)
//  {
//    if (!m_blockDecoupled)
//    {
//      //p
//      m_stateJacobian[0][6] += gamma*(mutTerm*overRT*overOmega - twoThirdduxduy*overRT/mutGaMod*gammaTerm) + pOmegaFactor*overRT*overOmega;
//  
//      //T
//      m_stateJacobian[4][6] += gamma*pOverRTT*(-mutTerm*overOmega + twoThirdduxduy/mutGaMod*gammaTerm) - pOmegaFactor*pOverRTT*overOmega;
//    }
// 
//    //logOmega
//    m_stateJacobian[6][6] +=  -gamma*rho*mutTerm*overOmega - pOmegaFactor*rho*overOmega;
//  }
//  
//  // P gamma
//  const CFreal cfg1 = 1.735;
//  const CFreal cfg2 = 5.45;
//  const CFreal cfg3 = 0.95375;
//  const CFreal cfg4 = 2.2;
//    
//  CFreal fg = 1.0;
//    
//  if (avGa < 0.45) fg -= exp(-cfg1*tan(cfg2*avGa-cfg3)-cfg4);
//    
//  const CFreal fMnsigmaTerm = 1.0+0.58*pow(MInfLocal,0.6);
//  const CFreal fMnsigma = 1.0/(fMnsigmaTerm*fMnsigmaTerm);
//    
//  CFreal dpds = 1.0/avV*(u*dpx + v*dpy + w*dpz);
//
//  CFreal PRC = 1.0;
//    
//  if (m_PGrad)
//  {
//    const CFreal lambdaTerm1 = -mu/(rho*rho*uInfLocal*uInfLocal*uInfLocal);
//
//    const CFreal dpdsLimit = -m_lambdaLim/(lambdaTerm1*ReThetat*ReThetat);
//    dpds = min(max(dpds,-dpdsLimit),dpdsLimit);
//    
//    const CFreal kPGrad = -muInfLocal/(rhoInfLocal*rhoInfLocal*uInfLocal*uInfLocal*uInfLocal)*fabs(1.0-MInfLocal*MInfLocal)*dpds;
//    
//      
//    if (kPGrad < 0)
//    {
//      const CFreal PRCfactor1 = 474.0*pow(TuInfLocal,-2.9);
//      const CFreal PRCexponent = 1.0-exp(2.0e6*kPGrad);
//      PRC = pow(PRCfactor1,PRCexponent);
//    }
//    else
//    {
//      const CFreal PRCexponent = -3227.0*pow(kPGrad,0.5985);
//      PRC = pow(10.0,PRCexponent);
//    }
//  }
//    
//  const CFreal fk = 1.0+fg*(PRC-1.0); 
//   
//  const CFreal nsigma = 1.25e-11*pow(TuInfLocal,7.0/4.0)*fk*fMnsigma;
//  
//  const CFreal betaGA = sqrt(nsigma)*uInfLocal*rhoInfLocal/muInfLocal;
//  
//  const CFreal gTerm = (1-avGa)*sqrt(-log(1.0-avGa));
//  
//  const CFreal brc2 = betaGA*rho*avV*2.0;
//  
//  CFreal prodTerm_Ga = fg*gTerm*brc2;
//  
//  if (prodTerm_Ga > 0.0)
//  {    
//    const CFreal dgTermdg = (2.0*log(1-avGa)+1.0)/(2.0*sqrt(-log(1-avGa))); 
//    
//    CFreal dfgdg = 0.0;
//    
//    if (avGa<0.45) dfgdg = cfg1*cfg2/pow(cos(cfg2*avGa-cfg3),2) * exp(cfg1*tan(cfg3-cfg2*avGa)-cfg4);
//    
//    //const CFreal dbdg = 0.5*uInfLocal*rhoInfLocal/muInfLocal/sqrt(nsigma) * 1.25e-11*pow(TuInfLocal,7.0/4.0) * fMnsigma * (PRC-1.0) * dfgdg;
//     
//    //gamma
//    //m_stateJacobian[6][6] += brc2 * fg * dgTermdg + brc2 * gTerm * dfgdg;// + gTerm * fg * 2.0 * rho * avV * dbdg;
//   
//    if (!m_blockDecoupled)
//    {
//      const CFreal gTermFgBeta2 = gTerm * fg * betaGA *2.0;
//      
//      //p
//      m_stateJacobian[0][7] +=  gTermFgBeta2 * overRT * avV;
//  
//      //T
//      m_stateJacobian[4][7] +=  -gTermFgBeta2 * pOverRTT * avV;
//      
//      //u
//      m_stateJacobian[1][7] +=  gTermFgBeta2 * rho * u/avV;
//              
//      //v
//      m_stateJacobian[2][7] +=  gTermFgBeta2 * rho * v/avV;    
//      
//      //w
//      m_stateJacobian[3][7] +=  gTermFgBeta2 * rho * w/avV;  
//    }
//  }
//  
//  // D gamma
//  const CFreal cEg = 20.0/0.57;
//  
//  const CFreal fMuGamma = 1.0-exp(-256.0*(m_currWallDist[iState]*uInfLocal*rhoInfLocal/muInfLocal)*(m_currWallDist[iState]*uInfLocal*rhoInfLocal/muInfLocal));
//  const CFreal fMMuGamma = (1.0+0.26*(gammaIsentropic-1.0)/2.0*MInfLocal*MInfLocal)*sqrt(1+0.38*pow(MInfLocal,0.6));
//  
//  //const CFreal gammaLim = std::min(std::max(0.01,avGa),0.99);
//  const CFreal muGamma = 0.57*pow(-log(1.0-avGa),-5.0/6.0*(1.0-avGa))*fMuGamma*fMMuGamma*mu;
//    
//  //const CFreal dudn = -1.0/(avV*avV)*(u*u*duy-v*v*dvx+u*v*(dvy-dux));//-duy; //
//  //const CFreal dgammadn = 1.0/avV*(v*dgammax - u*dgammay);//-dgammay; //
//  
//  const CFreal overVel = 1.0/avV;
//    const CFreal dVeldx = overVel*(u*dux+v*dvx+w*dwx);
//    const CFreal dVeldy = overVel*(u*duy+v*dvy+w*dwy);
//    const CFreal dVeldz = overVel*(u*duz+v*dvz+w*dwz);
//    
////    const CFreal dVeldx = overVel*(u*dux+v*dvx);
////    const CFreal dVeldy = overVel*(u*duy+v*dvy);
////    const CFreal dVeldz = overVel*(u*duz+v*dvz);
//    
//    CFreal dVeldnx;
//    CFreal dVeldny;
//    CFreal dVeldnz;
//    CFreal dVeldnSize;
//    
//    // temporary solution for cone case, code below based on largest du/dn, i.e. direction for shear vector, should work, but u grads in 3D seem too oscillatory in P1
//    const CFreal sinA = 0.139173101;//0.1218693434;    
//const CFreal cosA = 0.9902680687;//0.9925461516;
//const CFreal xvec = -sinA;
//    const CFreal zcurr = ((*m_cellStates)[iState]->getCoordinates())[ZZ];
//    const CFreal ycurr = ((*m_cellStates)[iState]->getCoordinates())[YY];
//    const CFreal nsizecurr = 1.0/sqrt(zcurr*zcurr+ycurr*ycurr);
//    const CFreal dVeldn1 = (dVeldy*ycurr*cosA+zcurr*dVeldz*cosA)*nsizecurr+xvec*dVeldx;
//    dVeldnSize = fabs(dVeldn1);
//    dVeldnx = xvec;
//    dVeldny = ycurr*nsizecurr*cosA;//dVeldn1*ycurr*nsizecurr/dVeldnSize;
//    dVeldnz = zcurr*nsizecurr*cosA;//dVeldn1*zcurr*nsizecurr/dVeldnSize;
//    
//    const CFreal dudn = -dVeldnSize;//1.0/(avV*avV)*(u*u*duy-v*v*dvx+u*v*(dvy-dux));//-duy; //
//    
//    const CFreal dgammadn = (dVeldnSize < 1.0e-7*avV) ? 0.0 : -(dVeldnx*dgammax + dVeldny*dgammay + dVeldnz*dgammaz);//-dgammay; //
//  
//  const CFreal dgTerm = -cEg * avV/(uInfLocal*uInfLocal) * dudn * dgammadn;
//    
//  CFreal  destructionTerm_Ga  = dgTerm * muGamma;
//    
//  //destructionTerm_Ga = max(-10.0*fabs(prodTerm_Ga),destructionTerm_Ga);
//  //if (avGa<0.4) destructionTerm_Ga = max(-fabs(prodTerm_Ga),destructionTerm_Ga);
//  
//  if (destructionTerm_Ga < 0.0 && m_addDGDA)
//  {
//    const CFreal mgdg = muGamma * 5.0/6.0 * (log(-log(1-avGa)) + 1.0/log(1-avGa));
//    
//    // gamma
//    //m_stateJacobian[6][6] +=  dgTerm * mgdg;
//  
//    if (!m_blockDecoupled)
//    {
//      // u
//      m_stateJacobian[1][7] +=  destructionTerm_Ga * u/(avV*avV);
//      
//      // v
//      m_stateJacobian[2][7] +=  destructionTerm_Ga * v/(avV*avV);  
//      
//      // w
//      m_stateJacobian[3][7] +=  destructionTerm_Ga * w/(avV*avV); 
//    }
//  }
//  
//  // production  Alpha
//  
//  const CFreal cpa1 = 0.03;
//  //const CFreal cpa2 = 50.0;
//
//  CFreal Rethetac = getRethetat(Tu, false);
//    
//  ///@todo check if this is local M or local MInf
//  Rethetac *= sqrt(1.0+0.38*pow(MInfLocal,0.6));
//    
//  CFreal lambda = 0.0;
//
//  if (m_PGrad)
//  {
//    lambda = -Rethetac*Rethetac*mu/(rho*rho*uInfLocal*uInfLocal*uInfLocal)*dpds;
//  }
//    
//  const CFreal alphac = rhoInfLocal*uInfLocal*uInfLocal*pow(0.09+lambda,0.62)/(Rethetac*sqrt(rhoInfLocal*muInfLocal));
//    
//  const CFreal t = (500.0 * mu )/(rho * avV * avV);
//  
//  //const CFreal dkdn = 1.0/avV*(v*dkx - u*dky);//-dky; //  
//  const CFreal dkdn = (dVeldnSize < 1.0e-7*avV) ? 0.0 : -(dVeldnx*dkx + dVeldny*dky + dVeldnz*dkz);
//    
//  const CFreal  Rew         = (rho * m_currWallDist[iState] * m_currWallDist[iState] * avOmega)/(mu);   
//  const CFreal  Fwake1      = (1.0e-5 * Rew)*(1.0e-5 * Rew);   
//  const CFreal  Fwake       = exp(-Fwake1);
//  const CFreal  thetaBL     = (ReThetat*mu)/(rho*avV);
//     
//  const CFreal  delta       = (375.0 * m_vorticity * m_currWallDist[iState] * thetaBL)/(avV);
//    
//  const CFreal  coefFtheta0 = (m_currWallDist[iState]/delta)*(m_currWallDist[iState]/delta)*(m_currWallDist[iState]/delta)*(m_currWallDist[iState]/delta);
//  const CFreal  coefFtheta1 = exp(-coefFtheta0);
//  const CFreal  Ftheta1     = Fwake * coefFtheta1;
//  const CFreal  Ftheta3     = 1.0-(((avGa-0.01)/(0.99-0.01))*((avGa-0.01)/(0.99-0.01)));
//  const CFreal  Ftheta4     = std::max(Ftheta1,Ftheta3);
//  CFreal  Fthetat     = std::min(Ftheta4,1.0);
//  
//  const CFreal dkdnFactor = fabs(dkdn)/max(avK,1.0e-6);
//
//  if (dkdnFactor>100.0) Fthetat = 1.0;
//    
//  CFreal prodTerm_alpha = cpa1 * (rho/t) * (alphac - avAlpha) * (1.0 - Fthetat);
//
//  // compute destruction term of alpha
//  const CFreal cna = 0.4;
//  const CFreal dadn = 1.0/avV*(v*dax - u*day);
//    
//  CFreal destructionTerm_alpha = -cna*rho*avV*dadn;
//  
//  ///@todo Fthetat and alpha_c comtribution
//  if (prodTerm_alpha > 0.0 || !m_limPRe)
//  {
//    // alpha
//    m_stateJacobian[8][8] +=  -cpa1 * (rho/t) * (1.0 - Fthetat);
//  
//    if (!m_blockDecoupled)
//    {
//      const CFreal tTerm = 2.0*rho*avV*avV/(500.0*mu);
//        
//      // p
//      m_stateJacobian[0][8] +=  cpa1 * tTerm * overRT * (alphac - avAlpha) * (1.0 - Fthetat);
//  
//      // T
//      m_stateJacobian[4][8] +=  -cpa1 * tTerm * pOverRTT * (alphac - avAlpha) * (1.0 - Fthetat);
//  
//      // u
//      m_stateJacobian[1][8] +=  cpa1 * 2.0*rho*rho*u/(500.0*mu) * (alphac - avAlpha) * (1.0 - Fthetat);
//  
//      // v
//      m_stateJacobian[2][8] +=  cpa1 * 2.0*rho*rho*v/(500.0*mu) * (alphac - avAlpha) * (1.0 - Fthetat);
//      
//      // w
//      m_stateJacobian[3][8] +=  cpa1 * 2.0*rho*rho*w/(500.0*mu) * (alphac - avAlpha) * (1.0 - Fthetat);
//    }
//  }
//  
//  if (destructionTerm_alpha < 0.0 && !m_blockDecoupled)
//  {
//      const CFreal daTerm = -cna*dadn;
//        
//      // p
//      m_stateJacobian[0][8] +=  daTerm*avV*overRT;
//  
//      // T
//      m_stateJacobian[4][8] +=  -daTerm*avV*pOverRTT;
//  
//      // u
//      m_stateJacobian[1][8] +=  daTerm*rho*u/avV;
//  
//      // v
//      m_stateJacobian[2][8] +=  daTerm*rho*v/avV;
//      
//      // w
//      m_stateJacobian[3][8] +=  daTerm*rho*w/avV;
//  }
//
//  if (m_addUpdateCoeff)
//  {
//    DataHandle<CFreal> updateCoeff = socket_updateCoeff.getDataHandle();
//  
//    CFreal update = 0.0;
//  
//    const CFreal PkUpdate = max(Pk,0.0);
//  
//    const CFreal DkUpdate = min(Dk,-1.0e-10);
//  
//    const CFreal G = sqrt(-PkUpdate*rho/DkUpdate);
//  
//    const CFreal Gk = pow(G,1.0/gamma);
//  
//    const CFreal kUpdateTerm = (PkUpdate + DkUpdate)/max(avK,0.001)/(Gk-1) + m_stateJacobian[5][5];
//  
//    update = max(0.0,kUpdateTerm);
//  
//    //if (update>1000) CFLog(INFO,"kUpdate: " << update << "\n");
//  
//    const CFreal logOmega = (*((*m_cellStates)[iState]))[6];
//  
//    const CFreal Gomega = pow(G,1.0/logOmega);
//  
//    const CFreal PomegaUpdate = max(Pomega,0.0);
//  
//    const CFreal DomegaUpdate = min(Domega,0.0);
//  
//    const CFreal omegaUpdateTerm = (PomegaUpdate + DomegaUpdate)/max(logOmega,0.1)/(Gomega-1) + m_stateJacobian[6][6];
//  
//    update = max(update,omegaUpdateTerm);
//  
//    // get the local ID of the current sol pnt
//    const CFuint solID = (*m_cellStates)[iState]->getLocalID();
//  
//    //if (update>1000) CFLog(INFO,"old: " << updateCoeff[solID] << ", new: " << update << ", J: " << m_solPntJacobDets[iState] << "\n\n");
//  
//    update *= (2.0*m_order+1)*m_solPntJacobDets[iState];//*m_solPntJacobDets[iState];
//  
//    if (update>0.5*updateCoeff[solID]) CFLog(INFO,"Large influence ST Dt: oldCoeff: " << updateCoeff[solID] << ", ST addition: " << update << ", J: " << m_solPntJacobDets[iState] << "\n\n");
//    
//    // add the wave speed update previously computed
//    updateCoeff[solID] += update;
//  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluxReconstructionMethod

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
