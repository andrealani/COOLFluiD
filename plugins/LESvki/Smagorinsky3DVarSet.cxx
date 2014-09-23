#include "Smagorinsky3DVarSet.hh"
#include "NavierStokes/EulerTerm.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Physics::NavierStokes;
//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace LESvki {

//////////////////////////////////////////////////////////////////////////////

void Smagorinsky3DVarSet::setup()
{
  NavierStokes::NavierStokesVarSet::setup();
}

//////////////////////////////////////////////////////////////////////////////

RealVector& Smagorinsky3DVarSet::getFlux(const RealVector& state,
                                          const vector<RealVector*>& gradients,
                                          const RealVector& normal,
                                          const CFreal& radius)
{
  
  CFuint dim = PhysicalModelStack::getActive()->getDim();
  setGradientState(state);
  //CF_DEBUG_OBJ("getFlux with normal");
  CFreal vOverRadius = 0.0;
  if (radius > MathTools::MathConsts::CFrealEps()) {
    // if the face is a boundary face, the radius could be 0
    // check against eps instead of 0. for safety
    vOverRadius = _gradState[2]/radius;
  }
 
  const RealVector& gradU = *gradients[1];
  const RealVector& gradV = *gradients[2];
  const RealVector& gradW = *gradients[3];  
  const RealVector& gradT = *gradients[4];
  //CF_DEBUG_OBJ(vOverRadius);
  const CFreal twoThirdDivV = 2./3.*(gradU[XX] + gradV[YY] +  gradW[ZZ]); // + vOverRadius); // Why vOverRadius??
  const CFreal avU = _gradState[1];
  const CFreal avV = _gradState[2];
  const CFreal avW = _gradState[3];  

  RealVector& nsData = getModel().getPhysicalData();

  // adimensional dynamical viscosity
  if (_useBackUpValues || _freezeDiffCoeff) {
    nsData[NSTerm::MU] = _dynViscCoeff;
    nsData[NSTerm::LAMBDA] = _thermCondCoeff;
  }
  else {
    // adimensional dynamical viscosity
    nsData[NSTerm::MU] = getDynViscosity(state, gradients)*getModel().getArtDiffCoeff();

    // adimensional thermal conductivity
    nsData[NSTerm::LAMBDA] = getThermConductivity
      (state, nsData[NSTerm::MU]);


    if (_setBackUpValues) {
      _dynViscCoeff = nsData[NSTerm::MU];
      _thermCondCoeff = nsData[NSTerm::LAMBDA];
    }
  }

  // CF_DEBUG_OBJ(nsData[NSTerm::MU]);
  const CFreal coeffTauMu = getModel().getCoeffTau()*nsData[NSTerm::MU];
  const CFreal tauXX = coeffTauMu*(2.*gradU[XX] - twoThirdDivV);
  const CFreal tauXY = coeffTauMu*(gradU[YY] + gradV[XX]);
  const CFreal tauXZ = coeffTauMu*(gradU[ZZ] + gradW[XX]);  
  const CFreal tauYZ = coeffTauMu*(gradW[YY] + gradV[ZZ]);  
  const CFreal tauYY = coeffTauMu*(2.*gradV[YY] - twoThirdDivV);
  const CFreal tauZZ = coeffTauMu*(2.*gradW[ZZ] - twoThirdDivV);  
  
  const CFreal nx = normal[XX];
  const CFreal ny = normal[YY];
  const CFreal nz = normal[ZZ];
  
  const CFreal qFlux = -getModel().getCoeffQ()*nsData[NSTerm::LAMBDA]*
    (gradT[XX]*nx + gradT[YY]*ny + gradT[ZZ]*nz);
  
  //Smagorinksy Subgrid Stress
  const CFreal fltRho = state[0]; // Filtered density
  const CFreal D = 3.0*radius; // Filter width Vreman: dx
  //CF_DEBUG_OBJ(D); 
  const CFreal SXX = gradU[XX];
  const CFreal SXY = 0.5*(gradU[YY]+gradV[XX]);
  const CFreal SXZ = 0.5*(gradU[ZZ]+gradW[XX]);
  const CFreal SYZ = 0.5*(gradW[YY]+gradV[ZZ]);  
  const CFreal SYY = gradV[YY];  
  const CFreal SZZ = gradW[ZZ];
  
 
  //Smagorinsky model:
  const CFreal Cs = 0.1; // Smagorinsky constant 0.1-0.2
  const CFreal mut = getTurbDynViscosityFromGradientVars(state, gradients);
  const CFreal nuSGS = pow(D,2)*mut; // Smagorinsky Subgrid viscosity
  //   CF_DEBUG_OBJ(nuSGS);

  //WALES and Smagorinsky model: 
  const CFreal tauSGSdXX = -2.*fltRho*nuSGS*(gradU[XX]-(gradU[XX]+gradV[YY]+gradW[ZZ])/dim);
  const CFreal tauSGSdXY = -2.*fltRho*nuSGS*(gradU[YY]+gradV[XX])/2.;
  const CFreal tauSGSdXZ = -2.*fltRho*nuSGS*(gradW[XX]+gradU[ZZ])/2.;  
  const CFreal tauSGSdYZ = -2.*fltRho*nuSGS*(gradW[YY]+gradV[ZZ])/2.;  
  const CFreal tauSGSdYY = -2.*fltRho*nuSGS*(gradV[YY]-(gradU[XX]+gradV[YY]+gradW[ZZ])/dim);
  const CFreal tauSGSdZZ = -2.*fltRho*nuSGS*(gradW[ZZ]-(gradU[XX]+gradV[YY]+gradW[ZZ])/dim); 
  
  
  // deviatoric part of Sugrid Stress tensor
  
  //Subgrid heat flux:
  const CFreal Pr_t = 0.6; // given in Comte/Lesieur. Turbulent Prandtl number
  const CFreal Cp = getModel().getCpOverPrandtl()*getModel().getPrandtl(); // specific heat at constant pressure
  const CFreal Qx = - fltRho * Cp * nuSGS / Pr_t *gradT[XX];
  const CFreal Qy = - fltRho * Cp * nuSGS / Pr_t *gradT[YY]; 
  const CFreal Qz = - fltRho * Cp * nuSGS / Pr_t *gradT[ZZ];   
  //Subgrid heat flux (definition of Comte/Lesieur)
 
  
  _flux[1] = tauXX*nx + tauXY*ny + tauXZ*nz - tauSGSdXX*nx - tauSGSdXY*ny - tauSGSdXZ*nz;
  _flux[2] = tauXY*nx + tauYY*ny + tauYZ*nz - tauSGSdXY*nx - tauSGSdYY*ny - tauSGSdYZ*nz;
  _flux[3] = tauXZ*nx + tauYZ*ny + tauZZ*nz - tauSGSdXZ*nx - tauSGSdYZ*ny - tauSGSdZZ*nz; 
  _flux[4] = (tauXX*avU + tauXY*avV + tauXZ*avW)*nx + (tauXY*avU + tauYY*avV + tauYZ*avW)*ny + (tauXZ*avU + tauYZ*avV + tauZZ*avW)*nz - qFlux - (Qx*nx + Qy*ny + Qz*nz);
  
  return _flux;
}

//////////////////////////////////////////////////////////////////////////////

RealMatrix& Smagorinsky3DVarSet::getFlux(const RealVector& state,
                                          const vector<RealVector*>& gradients,
                                          const CFreal& radius)
{
  setGradientState(state);
  CFuint dim = PhysicalModelStack::getActive()->getDim();
  CF_DEBUG_OBJ("getFlux without normal");
  CFreal vOverRadius = 0.0;
  if (radius > MathTools::MathConsts::CFrealEps()) {
    // if the face is a boundary face, the radius could be 0
    // check against eps instead of 0. for safety
    vOverRadius = _gradState[2]/radius;
  }

  const RealVector& gradU = *gradients[1];
  const RealVector& gradV = *gradients[2];
  const RealVector& gradW = *gradients[3];  
  const RealVector& gradT = *gradients[4];
  
  const CFreal twoThirdDivV = 2./3.*(gradU[XX] + gradV[YY] +  gradW[ZZ]);
  const CFreal avU = _gradState[1];
  const CFreal avV = _gradState[2];
  const CFreal avW = _gradState[3]; 

  RealVector& nsData = getModel().getPhysicalData();

  // adimensional dynamical viscosity
  if (_useBackUpValues || _freezeDiffCoeff) {
    nsData[NSTerm::MU] = _dynViscCoeff;
    nsData[NSTerm::LAMBDA] = _thermCondCoeff;
  }
  else {
    // adimensional dynamical viscosity
    nsData[NSTerm::MU] = getDynViscosity(state, gradients);

    // adimensional thermal conductivity
    nsData[NSTerm::LAMBDA] = getThermConductivity
      (state, nsData[NSTerm::MU]);


    if (_setBackUpValues) {
      _dynViscCoeff = nsData[NSTerm::MU];
      _thermCondCoeff = nsData[NSTerm::LAMBDA];
    }
  }

  CFreal Dynturbvisc = getTurbDynViscosityFromGradientVars(state, gradients);
  const CFreal coeffTauMu = getModel().getCoeffTau()*nsData[NSTerm::MU];
  const CFreal tauXX = coeffTauMu*(2.*gradU[XX] - twoThirdDivV);
  const CFreal tauXY = coeffTauMu*(gradU[YY] + gradV[XX]);
  const CFreal tauXZ = coeffTauMu*(gradU[ZZ] + gradW[XX]);  
  const CFreal tauYZ = coeffTauMu*(gradW[YY] + gradV[ZZ]);  
  const CFreal tauYY = coeffTauMu*(2.*gradV[YY] - twoThirdDivV);
  const CFreal tauZZ = coeffTauMu*(2.*gradW[ZZ] - twoThirdDivV); 
  const RealVector qFlux = -getModel().getCoeffQ()*nsData[NSTerm::LAMBDA]*gradT;
  
  //Smagorinksy Subgrid Stress
  const CFreal fltRho = state[0]; // Filtered density
  const CFreal D = radius; // Filter width
  const CFreal Cs = 0.18; // Smagorinsky constant
  
  const CFreal SXX = gradU[XX];
  const CFreal SXY = 0.5*(gradU[YY]+gradV[XX]);
  const CFreal SXZ = 0.5*(gradU[ZZ]+gradW[XX]);
  const CFreal SYZ = 0.5*(gradW[YY]+gradV[ZZ]);  
  const CFreal SYY = gradV[YY];  
  const CFreal SZZ = gradW[ZZ];

  const CFreal S = sqrt(2.*pow(SXX,2)+2.*pow(SYY,2)+2.*pow(SZZ,2)+pow(SXY,2)+pow(SXZ,2)+pow(SYZ,2));
  // sqrt(2)*Frobenius_Norm(Strain_rate_tensor)=sqrt(2SijSij)
  
  const CFreal nuSGS = pow(Cs,2)*pow(D,2)*S; // Smagorinsky Subgrid viscosity
  //const CFreal nuSGS = Dynturbvisc*pow(D,2); // Smagorinsky Subgrid viscosity
  const CFreal tauSGSdXX = -2.*fltRho*nuSGS*(gradU[XX]-(gradU[XX]+gradV[YY]+gradW[ZZ])/dim);
  const CFreal tauSGSdXY = -2.*fltRho*nuSGS*(gradU[YY]+gradV[XX])/2.;
  const CFreal tauSGSdXZ = -2.*fltRho*nuSGS*(gradW[XX]+gradU[ZZ])/2.;  
  const CFreal tauSGSdYZ = -2.*fltRho*nuSGS*(gradW[YY]+gradV[ZZ])/2.;  
  const CFreal tauSGSdYY = -2.*fltRho*nuSGS*(gradV[YY]-(gradU[XX]+gradV[YY]+gradW[ZZ])/dim);
  const CFreal tauSGSdZZ = -2.*fltRho*nuSGS*(gradW[ZZ]-(gradU[XX]+gradV[YY]+gradW[ZZ])/dim); 
  // deviatoric part of Sugrid Stress tensor
  
  //Subgrid heat flux
  const CFreal Pr_t = 0.6; // given in Comte/Lesieur. Turbulent Prandtl number
  const CFreal Cp = getModel().getCpOverPrandtl()*getModel().getPrandtl(); // specific heat at constant pressure
  const CFreal Qx = - fltRho * Cp * nuSGS / Pr_t *gradT[XX];
  const CFreal Qy = - fltRho * Cp * nuSGS / Pr_t *gradT[YY]; 
  const CFreal Qz = - fltRho * Cp * nuSGS / Pr_t *gradT[ZZ];   
  //Subgrid heat flux (definition of Comte/Lesieur)

  _fluxVec(1,XX) = tauXX - tauSGSdXX;
  _fluxVec(1,YY) = tauXY - tauSGSdXY;
  _fluxVec(1,ZZ) = tauXZ - tauSGSdXZ; 

  _fluxVec(2,XX) = tauXY - tauSGSdXY;
  _fluxVec(2,YY) = tauYY - tauSGSdYY;
  _fluxVec(2,ZZ) = tauYZ - tauSGSdYZ;
 
  _fluxVec(3,XX) = tauXZ - tauSGSdXZ;
  _fluxVec(3,YY) = tauYZ - tauSGSdYZ;
  _fluxVec(3,ZZ) = tauZZ - tauSGSdZZ; 

  _fluxVec(4,XX) = tauXX*avU + tauXY*avV + tauXZ*avW - qFlux[XX] - Qx;
  _fluxVec(4,YY) = tauXY*avU + tauYY*avV + tauYZ*avW - qFlux[YY] - Qy;
  _fluxVec(4,ZZ) = tauXZ*avU + tauYZ*avV + tauZZ*avW - qFlux[ZZ] - Qz;  

  return _fluxVec;
}

//////////////////////////////////////////////////////////////////////////////

void Smagorinsky3DVarSet::getAxiSourceTerm(const RealVector& physicalData,
					    const RealVector& state,
					    const vector<RealVector*>& gradients,
					    const CFreal& radius,
					    RealVector& source)
{
 
}

//////////////////////////////////////////////////////////////////////////////

CFreal Smagorinsky3DVarSet::getHeatFlux(const RealVector& state,
					 const vector<RealVector*>& gradients,
					 const RealVector& normal)
{
 
  const CFreal mu = getDynViscosity(state, gradients);
  const CFreal lambda = getThermConductivity(state, mu);
  const RealVector& gradT = *gradients[3];
  
  return (-getModel().getCoeffQ()*lambda*(gradT[XX]*normal[XX] + gradT[YY]*normal[YY]) + gradT[ZZ]*normal[ZZ]);
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NavierStokes

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
