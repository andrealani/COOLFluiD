#include "MHD/MHD.hh"
#include "MHD/MHD3DProjectionDiffPrimHW.hh"
#include "Environment/ObjectProvider.hh"
#include "MHD/MHDProjectionTerm.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace MHD {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<MHD3DProjectionDiffPrimHW, DiffusiveVarSet, MHDModule, 2> 
mhd3DProjectionDiffPrimHWProvider("MHD3DProjectionDiffPrimHW");
      
//////////////////////////////////////////////////////////////////////////////

MHD3DProjectionDiffPrimHW::MHD3DProjectionDiffPrimHW(const std::string& name,
						 Common::SafePtr<Framework::PhysicalModelImpl> model) :
  MHD3DProjectionDiffVarSet(name, model),
  _eulerModel(model->getConvectiveTerm().d_castTo<MHDProjectionTerm>()),
  _tempX()
{
  vector<std::string> names(9);
  names[0] = "rho";
  names[1] = "u";
  names[2] = "v";
  names[3] = "w";
  names[4] = "Bx";
  names[5] = "By";
  names[6] = "Bz";
  names[7] = "p";
  names[8] = "phi";
  setVarNames(names);
}

//////////////////////////////////////////////////////////////////////////////

MHD3DProjectionDiffPrimHW::~MHD3DProjectionDiffPrimHW()
{
}

//////////////////////////////////////////////////////////////////////////////

void MHD3DProjectionDiffPrimHW::setGradientVars(const vector<RealVector*>& states,
						RealMatrix& values,
						const CFuint stateSize)
{
  /*
    Here you have to compute multiple states [rho u v w Bx By Bz T phi] from the given states [rho u v w Bx By Bz p phi]
    states and values are matrices storing multiple states
  */


CFreal mu = 1.27;       // Mean molecular weight
CFreal mH = 1.67e-27;   // Mass hydrogen
CFreal kB = 1.38e-23;


//cout << "state size = " << stateSize << endl;

// First copy over all states
for (CFuint i = 0; i < 9; ++i) {
    for (CFuint j = 0; j < stateSize ; ++j) {
		values(i,j) = (*states[j])[i];
		
	  if (i == 7) {
	      values(i,j) = (*states[j])[i]/(*states[j])[0]; // T_adim = P[adim]/rho[adim]
	             } else {
          values(i,j) = (*states[j])[i];
                 }
      //std::cout << "i = " << i << "\n";
      //std::cout << "j = " << j << "\n";
    }
}


// Now overwrite values(7,j):
//for (CFuint j = 0; j < 6; ++j) {
   //values(7,j) = (*states[j])[7]*mu*mH/(2*(*states[j])[0]*kB;
   // Adimensional:
//   values(7,j) = (*states[j])[7]/(*states[j])[0]; // T_adim = P[adim]/rho[adim]
//}







/*
for (CFuint i = 0; i < nbValues; ++i) {
    for (CFuint j = 0; j < stateSize; ++j) {
      
      if (j==7) {
       values(i,j) = (*states[j])[i]*mu*mH/(2*(*states[0])[i]*kB);
      } else {
        values(i,j) = (*states[j])[i];
      }

    }
}
*/

  /// examples for Navier-Stokes
  // from [p u v w T] to [p u v w T] 
  /*for (CFuint i = 0; i < nbValues; ++i) {
    for (CFuint j = 0; j < stateSize; ++j) {
      values(i,j) = (*states[j])[i];
    }
    }*/

  // example
  // from [rho rhoU rhoV rhoW rhoE] to [p u v w T] 
  /*const CFreal R = _eulerModel->getR();
  const CFreal ovCv = (_eulerModel->getGamma() - 1.)/R;
  
  for (CFuint i = 0; i < stateSize; ++i) {
    const RealVector& state = *states[i];
    const CFreal ovRho = 1./state[0]; 
    values(1,i) = state[1]*ovRho;
    values(2,i) = state[2]*ovRho;
    values(3,i) = state[3]*ovRho;
    
    const CFreal V2 = values(1,i)*values(1,i) + 
      values(2,i)*values(2,i) + values(3,i)*values(3,i);
    values(4,i) = (state[4]*ovRho - 0.5*V2)*ovCv;
    values(0,i) = R*state[0]*values(4,i);
    }*/
  
}

//////////////////////////////////////////////////////////////////////////////

void MHD3DProjectionDiffPrimHW::setGradientVarGradients(const vector<RealVector*>& states,
						      const vector< vector<RealVector*> >& stateGradients,
						      vector< vector<RealVector*> >& gradVarGradients,
						      const CFuint stateSize)
{
 throw Common::NotImplementedException
   (FromHere(), "MHD3DProjectionDiffPrimHW::setGradientVarGradients()");
}
      
//////////////////////////////////////////////////////////////////////////////

void MHD3DProjectionDiffPrimHW::setStateGradients(const vector<RealVector*>& states,
						const vector< vector<RealVector*> >& gradVarGradients,
						vector< vector<RealVector*> >& stateGradients,
						const CFuint stateSize)
{
  /// PETER
  throw Common::NotImplementedException
   (FromHere(), "MHD3DProjectionDiffPrimHW::setStateGradients()");
}

//////////////////////////////////////////////////////////////////////////////

CFreal MHD3DProjectionDiffPrimHW::getDynViscosity(const RealVector& state,
						const vector<RealVector*>& gradients)
{
  throw Common::NotImplementedException
    (FromHere(), "MHD3DProjectionDiffPrimHW::getDynViscosity()");
  return 0;
}

//////////////////////////////////////////////////////////////////////////////

CFreal MHD3DProjectionDiffPrimHW::getDensity(const RealVector& state)
{
  //throw Common::NotImplementedException
  //  (FromHere(), "MHD3DProjectionDiffPrimHW::getDensity()");
  //return 0;
  return state[0];
}

//////////////////////////////////////////////////////////////////////////////

void MHD3DProjectionDiffPrimHW::setGradientState(const RealVector& state)
{
  cf_assert(_gradState.size() == state.size());

  /// examples for Navier-Stokes
  /// from [p u v w T] to [p u v w T] 
  // _gradState = state;

  /// from [rho rhoU rhoV rhoW rhoE] to [p u v w T] 
  /*
    const CFreal R = _eulerModel->getR();
    const CFreal cv = R/(_eulerModel->getGamma() - 1.);
    
    // _gradState = [p u v w T]
    _gradState[1] = state[1]/state[0];
    _gradState[2] = state[2]/state[0];
    _gradState[3] = state[3]/state[0];
    const CFreal V2 = _gradState[1]*_gradState[1] +
    _gradState[2]*_gradState[2] +
    _gradState[3]*_gradState[3];
    
    _gradState[4] = (state[4] - 0.5*state[0]*V2)/(state[0]*cv);
    _gradState[0] = R*state[0]*_gradState[4];
  */
  
  CFreal mu = 1.27;       // Mean molecular weight
  CFreal mH = 1.67e-27;   // Mass hydrogen
  CFreal kB = 1.38e-23;
  
  // Again, copy over the full state vector:
  _gradState = state;
  // And then overwrite _gratState[7] = T
  //_gradState[7] = state[7]*mu*mH/(2.0*state[0]*kB);
  // Adimensional:
  _gradState[7] = state[7]/state[0];
  // _gradState[7] = _gradState[7]/(2.0*state[0]) - state[7]/(2.0*state[0]*state[0])*_gradState[0];
}
      
//////////////////////////////////////////////////////////////////////////////

void MHD3DProjectionDiffPrimHW::computeFluxJacobian(const RealVector& state,
						  const RealVector& gradientJacob,
						  const RealVector& normal,
						  const CFreal& radius,
						  RealMatrix& fluxJacob)
{
  throw Common::NotImplementedException
    (FromHere(), "MHD3DProjectionDiffPrimHW::computeFluxJacobian()");
}

//////////////////////////////////////////////////////////////////////////////

void MHD3DProjectionDiffPrimHW::setComposition(const RealVector& state,
					       const bool isPerturb,
					       const CFuint iVar)
{
}

//////////////////////////////////////////////////////////////////////////////

RealVector& MHD3DProjectionDiffPrimHW::getFlux(const RealVector& state,
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
   const CFreal epsilon= 1.0e-16; // you should not add eps if Bnorm is > 0.
   const CFreal kpar_SI_new = 9.0e-7*std::pow(T0_new, 2.5)*1.0e-5; //k seems to be ~ erg/cm^1/ s K), so we apply conversion to SI.i
   const CFreal Bnorm = std::sqrt(std::pow(Bx,2) + std::pow(By,2) + std::pow(Bz,2)) + epsilon; // in Teslas
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
   // Mark 2025.03.10 by HP
   CFreal p_temp = state[7];
   CFreal pth = p_temp*P0;
   CFreal mu0 = 1.2566e-6;
   CFreal pmag = Bnorm*Bnorm / 2. / mu0;
   pmag = max(pmag,1.e-16); 
   CFreal plasmaBeta = pth / pmag;
   //CFreal Q_Betafactor= 1.0+std::tanh(1.0/plasmaBeta/100.0); 
   CFreal fac=1.0; //1.0/20.0;
   CFreal Q_Strenghfactor=1.0; //+max(std::tanh((Bnorm-1.5e-3)/5.e-4),0.0);
   //CFreal Q_Betafactor= 1.0+std::tanh(std::tanh((1.0/plasmaBeta/10.0-1.0)*fac));  
   CFreal Q_Betafactor= 1.0+std::tanh(fac/plasmaBeta/100.0); 
   CFreal qFlux_1 = 0.;
   CFreal qFlux_2 = 0.; 
   CFreal Va=Bnorm / pow((state[0]*rho0 * mu0), 0.5);
   //Q_Betafactor = 1.0;
   // Mark 2025.03.10 by HP	
   alpha_cond = 1.0;
   //   CFreal smooth_factor = 1.0/(1.0+std::pow((radius*l0-l0),2.0)/std::pow((10.0*l0-l0),2.0));
   qFlux_1=-commFactor*(qFlux_x*normal[XX] + qFlux_y*normal[YY] + qFlux_z*normal[ZZ]);
   qFlux_2=3.0/2.0*alpha_cond*state[0]*rho0/mu/mH*kB*T0_new*(Vx*normal[XX]+Vy*normal[YY]+Vz*normal[ZZ]);
   if (radius <= 10.0) {
     //qFlux = -commFactor*(qFlux_x*normal[XX] + qFlux_y*normal[YY] + qFlux_z*normal[ZZ]); // q cdot n, this is for old condu 
     //qFlux=qFlux_1*min(Va/Vnorm,1.0)+qFlux_2*(1.0-min(Va/Vnorm,1.0));
     //qFlux=qFlux_1*min(Va/Vnorm,1.0)*Q_Betafactor+qFlux_2*(1.0-min(Va/Vnorm,1.0));
     //qFlux=(qFlux_1*min(Va/Vnorm,1.0)+qFlux_2*(1.0-min(Va/Vnorm,1.0)))*Q_Betafactor
     qFlux=(qFlux_1*min(Va/Vnorm,1.0)+qFlux_2*(1.0-min(Va/Vnorm,1.0)))*Q_Betafactor*Q_Strenghfactor;
     //qFlux = qFlux_1;
   } else {
     //qFlux =3.0/2.0*alpha_cond*state[0]*rho0/mu/mH*kB*T0_new*(Vx*normal[XX]+Vy*normal[YY]+Vz*normal[ZZ]);
     //qFlux=qFlux_1*min(Va/Vnorm,1.0)+qFlux_2*(1.0-min(Va/Vnorm,1.0));
     //qFlux = qFlux*Q_Betafactor;
     qFlux=(qFlux_1*min(Va/Vnorm,1.0)+qFlux_2*(1.0-min(Va/Vnorm,1.0)))*Q_Betafactor*Q_Strenghfactor;
     qFlux =qFlux*(1.0-min((radius-10.0)/2.5,1.0))+min((radius-10.0)/2.5,1.0)*qFlux_2;
     //qFlux = qFlux_2;
     //qFlux=qFlux_1*min(Va/Vnorm,1.0)*Q_Betafactor+qFlux_2*(1.0-min(Va/Vnorm,1.0));
   }
   //qFlux = qFlux*Q_Betafactor;
   
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

} // namespace MHD

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
