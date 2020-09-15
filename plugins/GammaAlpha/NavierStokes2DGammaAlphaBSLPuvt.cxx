#include "GammaAlpha.hh"
#include "NavierStokes2DGammaAlphaBSLPuvt.hh"
#include "Environment/ObjectProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace GammaAlpha {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<NavierStokes2DGammaAlphaBSLPuvt, DiffusiveVarSet,
	       GammaAlphaModule, 2>
ns2DGammaAlphaBSLPuvtProvider("NavierStokes2DGammaAlphaBSLPuvt");

//////////////////////////////////////////////////////////////////////////////

NavierStokes2DGammaAlphaBSLPuvt::NavierStokes2DGammaAlphaBSLPuvt
(const std::string& name, SafePtr<PhysicalModelImpl> model) :
  NavierStokes2DGammaAlphaPuvt(name, model)
{
}

//////////////////////////////////////////////////////////////////////////////

NavierStokes2DGammaAlphaBSLPuvt::~NavierStokes2DGammaAlphaBSLPuvt()
{
}

//////////////////////////////////////////////////////////////////////////////

void NavierStokes2DGammaAlphaBSLPuvt::computeBlendingCoefFromGradientVars(const RealVector& state, const RealVector& gradK, const RealVector& gradOmega)
{

  const CFreal rho = std::max(0.,getDensity(state));
  const CFreal K = std::max(0., state[4]);
  const CFreal Omega = std::exp(state[5]);
  const CFreal overOmega = 1./Omega;
  const CFreal mu = getLaminarDynViscosityFromGradientVars(state);

  const CFreal distance = std::max(_wallDistance, 1.e-12);

  //Menet SST model from 2003 : The definition of CD_kw uses 1.e-10 rather than 1.-20
  const CFreal CD_kw = std::max(2.* rho * _sigmaOmega2 * (gradK[XX]*gradOmega[XX] + gradK[YY]*gradOmega[YY]), 1.e-10);//* overOmega

///@todo here be careful with the adimensionalization...check this!!!
  const CFreal arg1_1 = sqrt(K)/(0.09*Omega*distance);
  const CFreal arg1_2 = (500. * mu)/(rho*distance*distance*Omega);
  const CFreal arg1_3 = (4.*rho*_sigmaOmega2*K)/(CD_kw * distance * distance);

  const CFreal arg1 = min(max(arg1_1, arg1_2), arg1_3);
  const CFreal F1org=  tanh(arg1*arg1*arg1*arg1);
 // Corrected Belnding Coefficient
  const CFreal Ry= rho * distance * std::sqrt(K)/mu;
  const CFreal F2 = std::exp(-Ry/120);
  const CFreal F3bis = std::pow(F2,8);
  const CFreal F3    = exp(-F3bis);
 
 _blendingCoef1 = std::max(F1org,F3);
}

//////////////////////////////////////////////////////////////////////////////

CFreal NavierStokes2DGammaAlphaBSLPuvt::getTurbDynViscosityFromGradientVars(const RealVector& state, const vector<RealVector*>& gradients)
{
  cf_assert(_wallDistance >= 0.);
  //cf_assert(gradients.size() > 0);

  CFreal mut = 0.;
  if((_wallDistance > 0.) && (gradients.size() > 0))
  {
    CFreal mut = NavierStokes2DGammaAlphaPuvt::getTurbDynViscosityFromGradientVars(state, gradients);
    
    const CFreal avGa = std::min(std::max(state[6],0.01),0.99);

    // compute effective eddy viscosity
    const CFreal cmut1 = 0.5;
    const CFreal cmut2 = 0.25;
    const CFreal cmut3 = 0.1;
    mut *= (avGa+avGa*(1.0-avGa)*cmut1*(1.0+tanh((avGa-cmut2)/(cmut3))));
  }
    
  return mut;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace GammaAlpha

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

