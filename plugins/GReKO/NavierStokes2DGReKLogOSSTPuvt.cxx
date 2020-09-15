#include "GReKO.hh"
#include "NavierStokes2DGReKLogOSSTPuvt.hh"
#include "Environment/ObjectProvider.hh"
#include "Framework/MeshData.hh"

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

Environment::ObjectProvider<NavierStokes2DGReKLogOSSTPuvt, DiffusiveVarSet,
	       GReKOModule, 2>
ns2DGReKLogOSSTPuvtProvider("NavierStokes2DGReKLogOSSTPuvt");

//////////////////////////////////////////////////////////////////////////////

NavierStokes2DGReKLogOSSTPuvt::NavierStokes2DGReKLogOSSTPuvt
(const std::string& name, SafePtr<PhysicalModelImpl> model) :
  NavierStokes2DGReKLogOBSLPuvt(name, model)
{
}

//////////////////////////////////////////////////////////////////////////////

NavierStokes2DGReKLogOSSTPuvt::~NavierStokes2DGReKLogOSSTPuvt()
{
}

//////////////////////////////////////////////////////////////////////////////

void NavierStokes2DGReKLogOSSTPuvt::setModelCoefficients()
{

  NavierStokes2DGReKLogOPuvt::setModelCoefficients();

  //Modified k-omega coeficients for SST
  _sigmaK1 = 0.85;
 // _sigmaK2 = 1.0;
 // _sigmaOmega2 = 0.856;
 // _beta2 = 0.0828;
  
 // Menter SST model from 2003  
 // _gamma1= 5/9;
 // _gamma2= 0.44;
}

//////////////////////////////////////////////////////////////////////////////

CFreal NavierStokes2DGReKLogOSSTPuvt::getTurbDynViscosityFromGradientVars(const RealVector& state, const vector<RealVector*>& gradients)
{

  cf_assert(_wallDistance >= 0.);
///@todo put this here??
  cf_assert(gradients.size() > 0);
//if (_wallDistance == 0.) CFLog(INFO, "wallDist: " << _wallDistance << "\n");
  CFreal mut = 0.;
  if((_wallDistance > 0.) && (gradients.size() > 0))
  {
    const CFreal a1 = 0.31;

    const CFreal rho = getDensity(state);
    const CFreal K = std::max(state[4],0.0);
    const CFreal Omega = std::exp(state[5]);
    
    //const CFreal vorticity = fabs((*(gradients[2]))[XX] - (*(gradients[1]))[YY]);
    const CFreal gradU_X = (*(gradients[1]))[XX];
    const CFreal gradU_Y = (*(gradients[1]))[YY];
    const CFreal gradV_X = (*(gradients[2]))[XX];
    const CFreal gradV_Y = (*(gradients[2]))[YY];
    const CFreal gradSum = (gradU_Y + gradV_X);
    const CFreal strain = std::sqrt(2.*(gradU_X*gradU_X + 0.5*gradSum*gradSum + gradV_Y*gradV_Y));
    
    const CFreal mu = getLaminarDynViscosityFromGradientVars(state);

    ///@todo here there should be an adimensionalization coef
    const CFreal arg2_1 = 2. * sqrt(K) / (0.09 * Omega * _wallDistance);
    const CFreal arg2_2 = (500. * mu) / (rho * _wallDistance * _wallDistance * Omega);
    const CFreal arg2 = std::max(arg2_1,arg2_2);

    const CFreal F2 = tanh(arg2*arg2);

    mut = (a1 * rho * K )/ std::max(a1*Omega , strain*F2);
  }

  return mut;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace GReKO

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

