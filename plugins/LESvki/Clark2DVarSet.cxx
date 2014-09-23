#include "Clark2DVarSet.hh"
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

void Clark2DVarSet::setup()
{
  NavierStokes::NavierStokesVarSet::setup();
}

//////////////////////////////////////////////////////////////////////////////

RealVector& Clark2DVarSet::getFlux(const RealVector& state,
                                          const vector<RealVector*>& gradients,
                                          const RealVector& normal,
                                          const CFreal& radius)
{
  
  const CFuint dim = 3.0;//PhysicalModelStack::getActive()->getDim();
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
  const RealVector& gradT = *gradients[3];
  //CF_DEBUG_OBJ(vOverRadius);
  const CFreal twoThirdDivV = 2./3.*(gradU[XX] + gradV[YY]); // + vOverRadius); // Why vOverRadius??
  const CFreal avU = _gradState[1];
  const CFreal avV = _gradState[2];

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

   const CFreal coeffTauMu = getModel().getCoeffTau()*nsData[NSTerm::MU];
  const CFreal tauXX = coeffTauMu*(2.*gradU[XX] - twoThirdDivV);
  const CFreal tauXY = coeffTauMu*(gradU[YY] + gradV[XX]);
  const CFreal tauYY = coeffTauMu*(2.*gradV[YY] - twoThirdDivV);
  const CFreal nx = normal[XX];
  const CFreal ny = normal[YY];
  const CFreal qFlux = -getModel().getCoeffQ()*nsData[NSTerm::LAMBDA]*
    (gradT[XX]*nx + gradT[YY]*ny);


  const CFreal fltRho = state[0]; // Filtered density
  const CFreal D = 3.0*radius; // Filter width Vreman: dx
   //Smagorinsky model:
  const CFreal Cs = 0.1; // Smagorinsky constant 0.1-0.2
  const CFreal mut = getTurbDynViscosityFromGradientVars(state, gradients);
  const CFreal nuSGS = pow(D,2)*mut; // Smagorinsky Subgrid viscosity


 //  Clark model:
  const CFreal a11 = 1/12 * pow(D,2) * (gradU[XX]*gradU[XX]+gradU[YY]*gradU[YY]);
  const CFreal a12 = 1/12 * pow(D,2) * (gradU[XX]*gradV[XX]+gradU[YY]*gradV[YY]);
  const CFreal a22 = 1/12 * pow(D,2) * (gradV[XX]*gradV[XX]+gradV[YY]*gradV[YY]); 
  
  const CFreal tauSGSdXX = fltRho * a11 -2.*fltRho*nuSGS*(gradU[XX]-(gradU[XX]+gradV[YY])/dim);
  const CFreal tauSGSdXY = fltRho * a12 -2.*fltRho*nuSGS*(gradU[YY]+gradV[XX])/2.;
  const CFreal tauSGSdYY = fltRho * a22 -2.*fltRho*nuSGS*(gradV[YY]-(gradU[XX]+gradV[YY])/dim);

  
  //Subgrid heat flux:
  const CFreal Pr_t = 0.6; // given in Comte/Lesieur. Turbulent Prandtl number
  const CFreal Cp = getModel().getCpOverPrandtl()*getModel().getPrandtl(); // specific heat at constant pressure
  const CFreal Qx = - fltRho * Cp * nuSGS / Pr_t *gradT[XX];
  const CFreal Qy = - fltRho * Cp * nuSGS / Pr_t *gradT[YY]; 
  //Subgrid heat flux (definition of Comte/Lesieur)
  
  _flux[1] = tauXX*nx + tauXY*ny - tauSGSdXX*nx - tauSGSdXY*ny;
  _flux[2] = tauXY*nx + tauYY*ny - tauSGSdXY*nx - tauSGSdYY*ny;;
  _flux[3] = (tauXX*avU + tauXY*avV)*nx + (tauXY*avU + tauYY*avV)*ny - qFlux - (Qx*nx + Qy*ny);
  
  return _flux;
}

//////////////////////////////////////////////////////////////////////////////

RealMatrix& Clark2DVarSet::getFlux(const RealVector& state,
                                          const vector<RealVector*>& gradients,
                                          const CFreal& radius)
{
  setGradientState(state);
  CF_DEBUG_OBJ("getFlux without normal");
  CFreal vOverRadius = 0.0;
  if (radius > MathTools::MathConsts::CFrealEps()) {
    // if the face is a boundary face, the radius could be 0
    // check against eps instead of 0. for safety
    vOverRadius = _gradState[2]/radius;
  }

  const RealVector& gradU = *gradients[1];
  const RealVector& gradV = *gradients[2];
  const RealVector& gradT = *gradients[3];
  const CFreal twoThirdDivV = 2./3.*(gradU[XX] + gradV[YY] + vOverRadius);
  const CFreal avU = _gradState[1];
  const CFreal avV = _gradState[2];

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
  const CFreal tauYY = coeffTauMu*(2.*gradV[YY] - twoThirdDivV);
  const RealVector qFlux = -getModel().getCoeffQ()*nsData[NSTerm::LAMBDA]*gradT;
  
  //Smagorinksy Subgrid Stress
  const CFreal fltRho = state[0]; // Filtered density
  const CFreal D = radius; // Filter width
  const CFreal Cs = 0.18; // Smagorinsky constant
  const CFreal S = sqrt(2.*pow(gradU[XX],2)+2.*pow(gradV[YY],2)+pow(gradU[YY]+gradV[XX],2)); 
  // sqrt(2)*Frobenius_Norm(Strain_rate_tensor)=sqrt(2SijSij)
  const CFreal nuSGS = pow(Cs,2)*pow(D,2)*S; // Smagorinsky Subgrid viscosity
  //const CFreal nuSGS = Dynturbvisc*pow(D,2); // Smagorinsky Subgrid viscosity
  const CFreal tauSGSdXX = -2.*fltRho*nuSGS*(gradU[XX]-(gradU[XX]+gradV[YY])/3.);
  const CFreal tauSGSdXY = -2.*fltRho*nuSGS*(gradU[YY]+gradV[XX])/2.;
  const CFreal tauSGSdYY = -2.*fltRho*nuSGS*(gradV[YY]-(gradU[XX]+gradV[YY])/3.);
  // deviatoric part of Sugrid Stress tensor
  
  //Subgrid heat flux
  const CFreal Pr_t = 0.6; // given in Comte/Lesieur. Turbulent Prandtl number
  const CFreal Cp = getModel().getCpOverPrandtl()*getModel().getPrandtl(); // specific heat at constant pressure
  const CFreal Qx = fltRho * Cp * nuSGS / Pr_t *gradT[XX];
  const CFreal Qy = fltRho * Cp * nuSGS / Pr_t *gradT[YY]; 
  //Subgrid heat flux (definition of Comte/Lesieur)

  _fluxVec(1,XX) = tauXX - tauSGSdXX;
  _fluxVec(1,YY) = tauXY - tauSGSdXY;

  _fluxVec(2,XX) = tauXY - tauSGSdXY;
  _fluxVec(2,YY) = tauYY - tauSGSdYY;

  _fluxVec(3,XX) = tauXX*avU + tauXY*avV - qFlux[XX] - Qx;
  _fluxVec(3,YY) = tauXY*avU + tauYY*avV - qFlux[YY] - Qy;

  return _fluxVec;
}

//////////////////////////////////////////////////////////////////////////////

void Clark2DVarSet::getAxiSourceTerm(const RealVector& physicalData,
					    const RealVector& state,
					    const vector<RealVector*>& gradients,
					    const CFreal& radius,
					    RealVector& source)
{
 
}

//////////////////////////////////////////////////////////////////////////////

CFreal Clark2DVarSet::getHeatFlux(const RealVector& state,
					 const vector<RealVector*>& gradients,
					 const RealVector& normal)
{
 const CFreal mu = getDynViscosity(state, gradients);
  const CFreal lambda = getThermConductivity(state, mu);
  const RealVector& gradT = *gradients[3];
  
  return (-getModel().getCoeffQ()*lambda*(gradT[XX]*normal[XX] + gradT[YY]*normal[YY]));
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NavierStokes

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
