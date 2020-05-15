#include "GReKO.hh"
#include "NavierStokes2DGReKLogOBSLPuvt.hh"
#include "Environment/ObjectProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace GReKO {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<NavierStokes2DGReKLogOBSLPuvt, DiffusiveVarSet,
	       GReKOModule, 2>
ns2DGReKLogOBSLPuvtProvider("NavierStokes2DGReKLogOBSLPuvt");

//////////////////////////////////////////////////////////////////////////////

NavierStokes2DGReKLogOBSLPuvt::NavierStokes2DGReKLogOBSLPuvt
(const std::string& name, SafePtr<PhysicalModelImpl> model) :
  NavierStokes2DGReKLogOPuvt(name, model)
{
}

//////////////////////////////////////////////////////////////////////////////

NavierStokes2DGReKLogOBSLPuvt::~NavierStokes2DGReKLogOBSLPuvt()
{
}

//////////////////////////////////////////////////////////////////////////////

void NavierStokes2DGReKLogOBSLPuvt::computeBlendingCoefFromGradientVars(const RealVector& state, const RealVector& gradK, const RealVector& gradOmega)
{

  const CFreal rho = getDensity(state);
  const CFreal K = std::max(0., state[4]);
  const CFreal Omega = std::exp(state[5]);
  const CFreal overOmega = 1./Omega;
  const CFreal mu = getLaminarDynViscosityFromGradientVars(state);

  const CFreal distance = std::max(_wallDistance, 1.e-12);

  //Menet SST model from 2003 : The definition of CD_kw uses 1.e-10 rather than 1.-20
  const CFreal CD_kw = std::max(2.* rho * _sigmaOmega2 * (gradK[XX]*gradOmega[XX] + gradK[YY]*gradOmega[YY]), 1.e-20);//* overOmega

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

    } // namespace GReKO

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

