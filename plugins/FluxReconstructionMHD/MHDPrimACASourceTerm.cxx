#include "FluxReconstructionMHD/MHDPrimACASourceTerm.hh"
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

MethodCommandProvider<MHDPrimACASourceTerm, FluxReconstructionSolverData, FluxReconstructionMHDModule>
MHDPrimACASourceTermFRProvider("MHDPrimACASourceTerm");

///////////////////////////////////////////////////////////////////////////////////////////////////

void MHDPrimACASourceTerm::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< bool >("Gravity","Switch on gravity.");
  options.addConfigOption< bool >("Rotation","Switch on rotation.");
}

/////////////////////////////////////////////////////////////////////////////////////////////////////

MHDPrimACASourceTerm::MHDPrimACASourceTerm(const std::string& name) :
  StdSourceTerm(name),
  socket_updateCoeff("updateCoeff"),
  socket_gravity("gravity"),
  m_order()
{ 
  addConfigOptionsTo(this);
  
  m_gravity = false; 
  setParameter("Gravity",&m_gravity);

  m_rotation = true;
  setParameter("Rotation",&m_rotation);

}

/////////////////////////////////////////////////////////////////////////////////////////////////

MHDPrimACASourceTerm::~MHDPrimACASourceTerm()
{
}

//////////////////////////////////////////////////////////////////////////////////////////////////

void MHDPrimACASourceTerm::setup()
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

}

//////////////////////////////////////////////////////////////////////////////////////////////////

void MHDPrimACASourceTerm::unsetup()
{
  StdSourceTerm::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

std::vector< Common::SafePtr< BaseDataSocketSource > >
  MHDPrimACASourceTerm::providesSockets()
{
  std::vector< Common::SafePtr< BaseDataSocketSource > > result = StdSourceTerm::providesSockets();
  result.push_back(&socket_gravity);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

std::vector< Common::SafePtr< BaseDataSocketSink > >
MHDPrimACASourceTerm::needsSockets()
{
  std::vector< Common::SafePtr< BaseDataSocketSink > > result = StdSourceTerm::needsSockets();
  result.push_back(&socket_updateCoeff);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

void MHDPrimACASourceTerm::addSourceTerm(RealVector& resUpdates)
{       
  CFLog(VERBOSE, "MHDPrimACASourceTerm::addSourceTerm() => START\n");
  
  // set gradients
  const CFuint nbrStates = m_cellStates->size();
   
  const CFreal RSun = 6.9551e8; // m
  const CFreal GMsun = 1.327474512e20; // SI value
  const CFreal gsun = 274.00; // at the surface m per sec2

  const CFreal lRef = 6.9e8;
  
  // Boltzmann constant
  const CFreal k = 1.3806503e-23;

  // mass of proton
  const CFreal mp = 1.67262158e-27;

  // mass of electron
  const CFreal me = 9.10938188e-31;

  //const CFreal rhoRef = nRef*(mp+me);
  const CFreal rhoRef = 1.67e-13;
  
  const CFreal mu_cor = 1.27; // Molecular weight in the corona
  const CFreal mH = 1.6733e-27; // kg
  
  const CFreal l0 = 6.9e8;
  const CFreal B0 = 2.2e-4;
  const CFreal mu0 = 1.2566e-6;
  const CFreal rho0 = 1.67e-13;
  const CFreal V0 = B0/std::sqrt(mu0*rho0);
  const CFreal g0 = pow(V0,2)/l0;

  // --- Precompute constants once (outside iSol loop) ---
  const CFreal Omega_phys = 2.97e-6; //2.866e-6; //2.0*MathTools::MathConsts::CFrealPi()/(25.0*24.0*3600.0);    // sidereal rotation [rad/s]
  const CFreal v_ref      = 4.8e5;                         // m/s
  const CFreal Omega_nd   = Omega_phys * RSun / v_ref;    // ~4.15e-3 (nondim)
  
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
    const CFreal density = (*((*m_cellStates)[iSol]))[0]*rhoRef;
    const CFreal density_dimless = (*((*m_cellStates)[iSol]))[0];

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
      
    //Rotation
    if (m_rotation)
    {
      // Coriolis accel (a_cor = -2 Ω × v)
      const CFreal cor_x = -2.0 * (Omega_nd * Vy);
      const CFreal cor_y = +2.0 * (Omega_nd * Vx);
      const CFreal cor_z =  0.0;

      // Centrifugal accel (a_cent = -Ω^2 (x,y,0))
      const CFreal cent_x = - (Omega_nd*Omega_nd) * x_dimless;
      const CFreal cent_y = - (Omega_nd*Omega_nd) * y_dimless;
      const CFreal cent_z =  0.0;

      // momentum
      resUpdates[m_nbrEqs*iSol + 1] -= ( density_dimless * cor_x
                                    + density_dimless * cent_x);
      resUpdates[m_nbrEqs*iSol + 2] -= ( density_dimless * cor_y
                                    + density_dimless * cent_y);
      resUpdates[m_nbrEqs*iSol + 3] -= ( density_dimless * cor_z
                                    + density_dimless * cent_z);

      // energy: v·(a_cent)
      const CFreal Vdotrot = Vx*(cent_x)
                          + Vy*(cent_y)
                          + Vz*(cent_z);
      resUpdates[m_nbrEqs*iSol + 7] -= density_dimless * Vdotrot;    
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
      
    }
  }
  
  CFLog(VERBOSE, "MHDPrimACASourceTerm::addSourceTerm() => END\n");
}

//////////////////////////////////////////////////////////////////////////////////////////////////

void  MHDPrimACASourceTerm::getSToStateJacobian(const CFuint iState)
{
  // reset the jacobian
  for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
  {
    m_stateJacobian[iEq] = 0.0;
  }
  
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
  
  const CFreal density_dimless = (*((*m_cellStates)[iState]))[0];
  const CFreal Vx = (*((*m_cellStates)[iState]))[1];
  const CFreal Vy = (*((*m_cellStates)[iState]))[2];
  const CFreal Vz = (*((*m_cellStates)[iState]))[3];
  
  // Rotation constants (same as in addSourceTerm)
  const CFreal Omega_phys = 2.866e-6; // sidereal rotation [rad/s]
  const CFreal v_ref      = 4.8e5;     // m/s
  const CFreal Omega_nd   = Omega_phys * RSun / v_ref; // nondimensional
  
  /////// Gravity
  if (m_gravity)
  {
    // V
    m_stateJacobian[0][1] = gx;
    m_stateJacobian[0][2] = gy;
    m_stateJacobian[0][3] = gz;
    
    // p
    const CFreal Vdotg = Vx*gx + Vy*gy + Vz*gz; // dimensionless
      
    m_stateJacobian[0][7] = Vdotg;
    m_stateJacobian[1][7] = density_dimless*gx;
    m_stateJacobian[2][7] = density_dimless*gy;
    m_stateJacobian[3][7] = density_dimless*gz;
  }
  
  /////// Rotation (Coriolis and Centrifugal forces)
  if (m_rotation)
  {
    // Coriolis force derivatives: F_cor = -2*rho*Omega x v
    // dF_cor_x/dVy = -2*rho*Omega_z = -2*rho*Omega_nd
    // dF_cor_y/dVx = +2*rho*Omega_z = +2*rho*Omega_nd
    m_stateJacobian[1][2] += -2.0 * density_dimless * Omega_nd; // dF_y/dVx (Coriolis)
    m_stateJacobian[2][1] += 2.0 * density_dimless * Omega_nd;  // dF_x/dVy (Coriolis)
    
    // Derivative w.r.t. density for Coriolis force
    const CFreal cor_x = -2.0 * (Omega_nd * Vy);
    const CFreal cor_y = +2.0 * (Omega_nd * Vx);
    
    // Centrifugal force derivatives: F_cent = -rho*Omega^2*(x,y,0)
    // These are independent of velocity, only depend on density
    const CFreal cent_x = - (Omega_nd*Omega_nd) * x_dimless;
    const CFreal cent_y = - (Omega_nd*Omega_nd) * y_dimless;

    m_stateJacobian[0][1] += -cor_x+cent_x;
    m_stateJacobian[0][2] += -cor_y+cent_y;
    
    // Energy equation derivatives from rotational work
    // E_rot = -rho * v·a_cent = -rho * (Vx*cent_x + Vy*cent_y)
    const CFreal Vdotrot = Vx*cent_x + Vy*cent_y;
    
    // Derivatives w.r.t. density
    m_stateJacobian[0][7] += -Vdotrot; // dE_rot/drho
    
    // Derivatives w.r.t. velocity components
    m_stateJacobian[1][7] += -density_dimless * cent_x; // dE_rot/dVx
    m_stateJacobian[2][7] += -density_dimless * cent_y; // dE_rot/dVy
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluxReconstructionMethod

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
