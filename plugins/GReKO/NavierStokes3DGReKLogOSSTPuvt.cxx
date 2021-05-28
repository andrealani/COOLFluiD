#include "GReKO.hh"
#include "NavierStokes3DGReKLogOSSTPuvt.hh"
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

Environment::ObjectProvider<NavierStokes3DGReKLogOSSTPuvt, DiffusiveVarSet,
	       GReKOModule, 2>
ns3DGReKLogOSSTPuvtProvider("NavierStokes3DGReKLogOSSTPvt");

//////////////////////////////////////////////////////////////////////////////

NavierStokes3DGReKLogOSSTPuvt::NavierStokes3DGReKLogOSSTPuvt
(const std::string& name, SafePtr<PhysicalModelImpl> model) :
  NavierStokes3DGReKLogOBSLPuvt(name, model)
{
}

//////////////////////////////////////////////////////////////////////////////

NavierStokes3DGReKLogOSSTPuvt::~NavierStokes3DGReKLogOSSTPuvt()
{
}

//////////////////////////////////////////////////////////////////////////////

void NavierStokes3DGReKLogOSSTPuvt::setModelCoefficients()
{

  NavierStokes3DGReKLogOPuvt::setModelCoefficients();

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

CFreal NavierStokes3DGReKLogOSSTPuvt::getTurbDynViscosityFromGradientVars(const RealVector& state, const vector<RealVector*>& gradients)
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
    const CFreal K = std::max(state[5],0.0);
    const CFreal Omega = std::exp(state[6]);
    
    //const CFreal vorticity = fabs((*(gradients[2]))[XX] - (*(gradients[1]))[YY]);
    const CFreal gradU_X = (*(gradients[1]))[XX];
    const CFreal gradU_Y = (*(gradients[1]))[YY];
    const CFreal gradU_Z = (*(gradients[1]))[ZZ];
    const CFreal gradV_X = (*(gradients[2]))[XX];
    const CFreal gradV_Y = (*(gradients[2]))[YY];
    const CFreal gradV_Z = (*(gradients[2]))[ZZ];
    const CFreal gradW_X = (*(gradients[3]))[XX];
    const CFreal gradW_Y = (*(gradients[3]))[YY];
    const CFreal gradW_Z = (*(gradients[3]))[ZZ];
    
    const CFreal gradSumXY = (gradU_Y + gradV_X);
    const CFreal gradSumXZ = (gradU_Z + gradW_X);
    const CFreal gradSumYZ = (gradW_Y + gradV_Z);
      
    const CFreal strain = std::sqrt(2.*(gradU_X*gradU_X + gradV_Y*gradV_Y + gradW_Z*gradW_Z) + gradSumXY*gradSumXY + gradSumXZ*gradSumXZ + gradSumYZ*gradSumYZ);
    
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

