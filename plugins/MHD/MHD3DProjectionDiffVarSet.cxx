#include "MHD/MHD3DProjectionDiffVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace MHD {

//////////////////////////////////////////////////////////////////////////////

void MHD3DProjectionDiffVarSet::setup()
{
  MHDProjectionDiffVarSet::setup();
}

//////////////////////////////////////////////////////////////////////////////

RealVector& MHD3DProjectionDiffVarSet::getFlux(const RealVector& state,
					  const vector<RealVector*>& gradients,
					  const RealVector& normal,
					  const CFreal& radius)
{
  setGradientState(state);
  // AL: for now no viscosity
  (getModel().getPhysicalData())[MHDProjectionDiffTerm::MU] = 0.;
      
  /*
  computeTransportProperties(state, gradients, normal);
  computeStressTensor(state, gradients, radius);
  
  const RealVector& nsData = getModel().getPhysicalData();
  const CFreal qFlux = -getModel().getCoeffQ()*nsData[NSTerm::LAMBDA]*
    MathFunctions::innerProd(*gradients[_TID], normal);
  
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  _flux.slice(_uID, dim) = _tau*normal;
  
  _flux[_TID] = MathFunctions::innerProd(_flux.slice(_uID,dim), _gradState.slice(_uID,dim)) - qFlux;
  */
   
   CFreal qFlux = 0.; // PETER q*n
   CFreal qFlux_me = 0;
   CFreal mu = 1.27;       // Mean molecular weight
   CFreal mH = 1.67e-27;   // Mass hydrogen
   CFreal kB = 1.38e-23;

   CFreal Bx = state[4]*2.2e-4; // in T
   CFreal By = state[5]*2.2e-4;
   CFreal Bz = state[6]*2.2e-4;

   const CFreal V0 = 480363.085276; //480248.389661; // m/s
   const CFreal l0 = 6.95e8; // mi
   const CFreal l0_cgs = 6.95e10; // solar radius in cm
   const CFreal P0 = 0.03851; // N/m**2
   const CFreal rho0 = 1.67e-13; // kg/m**3
   const CFreal me = 9.10938188e-31;
   CFreal T0 = P0*mH/(2*rho0*kB); // T0 as normalization for temperature
   CFreal T0_new = state[7]*P0*mH/(2*state[0]*rho0*kB); // Temperature in the domain as a full variable
   CFreal Vx = state[1]*V0;
   CFreal Vy = state[2]*V0;
   CFreal Vz = state[3]*V0;

   // energy source term base unit S0:
   const CFreal S0 = rho0*std::pow(V0,3)/l0;
   // Heat flux base unit q0: q = q_adim*q0:
   CFreal q0 = l0*S0;

   const CFreal kpar_cgs = 9.0e-12*std::pow(state[7]*T0/state[0],2.5); // in erg*cm^2/(s K); numerical value from [Mikic et al 1999]
   // Unit conversion of thermal conductivity:
   // erg/ cm^3/(s K) = 10**-7/10**-6 J*m/(s K)

   const CFreal kpar_SI_new = 9.0e-7*std::pow(T0_new, 2.5)*1.0e-5; //k seems to be ~ erg/cm^1/ s K), so we apply conversion to SI.i
   const CFreal Bnorm = std::sqrt(std::pow(Bx,2) + std::pow(By,2) + std::pow(Bz,2)); // in Teslas
   const CFreal Vnorm = std::sqrt(std::pow(Vx,2) + std::pow(Vy,2) + std::pow(Vz,2)); // in m/s
   const CFreal commFactor = kpar_SI_new/(Bnorm*Bnorm); // in SI

   // The unit of Tesla**2 cancels from the expression of the heat flux; one could alternatively skip
   // the multiplication *2.2e-4 in the definition of Bx, By, Bz entirely here to avoid large/small numbers in the numerics.

   CFreal dTdx = (*gradients[7])[XX]*T0/l0; // in Kelvin/meters
   CFreal dTdy = (*gradients[7])[YY]*T0/l0; // in K/m
   CFreal dTdz = (*gradients[7])[ZZ]*T0/l0; // in K/m


   // Compute the dimensional heat flux vector components_

   CFreal qFlux_x = Bx*Bx*dTdx + Bx*By*dTdy + Bx*Bz*dTdz;
   CFreal qFlux_y = By*Bx*dTdx + By*By*dTdy + By*Bz*dTdz;
   CFreal qFlux_z = Bz*Bx*dTdx + Bz*By*dTdy + Bz*Bz*dTdz;
   
   CFreal alpha_cond = 0.;
   //CFreal mach = Vnorm* std::sqrt(mu*mH/(2*kB*T0_new));
   CFreal mach = Vnorm* std::sqrt(me/(2*kB*T0_new)); // mach numbe should be for electron, but number density should be for proton
   if (mach <= 0.0249) {
       alpha_cond = std::pow(mach, -0.2112);
   } else if ((mach > 0.0249) && (mach <= 0.3146)) {
       alpha_cond = 0.436* std::pow(mach,-0.436);
   } else if ((mach > 0.3146)) {
       alpha_cond = 0.035*std::pow(mach,-2.617);
   }

// old thermal conduction, collision limit at 10 solar radius, Hollweg 1978

   alpha_cond = 1.0;
//   CFreal smooth_factor = 1.0/(1.0+std::pow((radius*l0-l0),2.0)/std::pow((10.0*l0-l0),2.0));
   if (radius <= 10.0) {
     qFlux = -commFactor*(qFlux_x*normal[XX] + qFlux_y*normal[YY] + qFlux_z*normal[ZZ]); // q cdot n, this is for old condu 
   } else {
    qFlux =3.0/2.0*alpha_cond*state[0]*rho0/mu/mH*kB*T0_new*(Vx*normal[XX]+Vy*normal[YY]+Vz*normal[ZZ]);
   }

  
// old thermal conduction, smoothed transition at 10 R_sun, combination of Hollweg 1978 and Reville 2020

   
//   CFreal smooth_factor = 1.0/(1.0+std::pow((radius*l0-l0),2.0)/std::pow((10.0*l0-l0),2.0));
//   CFreal qflux_collisional = 0.0;
//   CFreal qflux_collisionless = 0.0;
//   qflux_collisional = -commFactor*(qFlux_x*normal[XX] + qFlux_y*normal[YY] + qFlux_z*normal[ZZ]);
//   qflux_collisionless = 3.0/2.0*alpha_cond*state[0]*rho0/mu/mH*kB*T0_new*Vnorm;
//   qFlux = (smooth_factor*qflux_collisional + (1.0-smooth_factor)*qflux_collisionless);


// thermal conduction without magnetic field, Reville 2020

//   alpha_cond = 1.0;
//   CFreal q_p = 3.0/2.0*alpha_cond*state[0]*rho0/mu/mH*kB*T0_new*Vnorm;
//   CFreal q_s = -kpar_SI_new*(dTdx*normal[XX] + dTdy*normal[YY]+dTdz*normal[ZZ]);
//   CFreal smooth_factor = 1.0/(1.0+std::pow((radius*l0-l0),4.0)/std::pow((10.0*l0-l0),4.0)); //definition in Reville et al 2020.
 //  qFlux =smooth_factor*q_s + (1.0-smooth_factor)*q_p;
   //qFlux =mu*mH/state[0]/rho0/kB* smooth_factor*q_s + (1.0-smooth_factor)*q_p;    



 _flux[7] = -qFlux/q0; 
 
 //std::cout << _flux[7] << endl; // The heat flux qFlux has a negative sign: - kappa nabla T. @Andrea: is this the same minus sign
                         // as in your template here _flux ... = - qFlux; or is this an additional minus sign coming
                         // from the geometry and orientation of the faces???
  return _flux;
}

//////////////////////////////////////////////////////////////////////////////

RealMatrix& MHD3DProjectionDiffVarSet::getFlux(const RealVector& state,
                                          const vector<RealVector*>& gradients,
                                          const CFreal& radius)
{
  cf_assert(PhysicalModelStack::getActive()->getDim() == 3);
  
  setGradientState(state);

  /*
  // dummy normal here
  computeTransportProperties(state, gradients, _normal);
  computeStressTensor(state, gradients, radius);
  
  RealVector& nsData = getModel().getPhysicalData();
  const RealVector& gradT = *gradients[_TID];
  const CFreal avU = _gradState[_uID];
  const CFreal avV = _gradState[_vID];
  const CFreal avW = _gradState[_wID];
  const CFreal tauXX = _tau(XX, XX);
  const CFreal tauXY = _tau(XX, YY);
  const CFreal tauYY = _tau(YY, YY);
  const CFreal tauXZ = _tau(XX, ZZ);
  const CFreal tauYZ = _tau(YY, ZZ);
  const CFreal tauZZ = _tau(ZZ, ZZ);
  const RealVector qFlux = -getModel().getCoeffQ()*nsData[NSTerm::LAMBDA]*gradT;
  
  _fluxVec(_uID,XX) = tauXX;
  _fluxVec(_uID,YY) = tauXY;
  _fluxVec(_uID,ZZ) = tauXZ;
  
  _fluxVec(_vID,XX) = tauXY;
  _fluxVec(_vID,YY) = tauYY;
  _fluxVec(_vID,ZZ) = tauYZ;
  
  _fluxVec(_wID,XX) = tauXZ;
  _fluxVec(_wID,YY) = tauYZ;
  _fluxVec(_wID,ZZ) = tauZZ;
  
  _fluxVec(_TID,XX) = tauXX*avU + tauXY*avV + tauXZ*avW - qFlux[XX];
  _fluxVec(_TID,YY) = tauXY*avU + tauYY*avV + tauYZ*avW - qFlux[YY];
  _fluxVec(_TID,ZZ) = tauXZ*avU + tauYZ*avV + tauZZ*avW - qFlux[ZZ];
  */

  throw Common::NotImplementedException(FromHere(), "MHD3DProjectionDiffVarSet::getFlux()");
  return _fluxVec;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace MHD

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
