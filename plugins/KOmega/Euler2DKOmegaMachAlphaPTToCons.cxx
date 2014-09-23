#include "NavierStokes/NavierStokes.hh"
#include "Euler2DKOmegaMachAlphaPTToCons.hh"
#include "Framework/PhysicalModel.hh"
#include "Environment/ObjectProvider.hh"
#include "EulerPhysicalModel.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace NavierStokes {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<Euler2DKOmegaMachAlphaPTToCons, VarSetTransformer, NavierStokesModule, 1> 
euler2DMachAlphaPTToConsProvider("Euler2DKOmegaMachAlphaPTToCons");

//////////////////////////////////////////////////////////////////////////////

Euler2DKOmegaMachAlphaPTToCons::Euler2DKOmegaMachAlphaPTToCons(Common::SafePtr<Framework::PhysicalModelImpl> model) :
  VarSetTransformer(model),
 _model(model->getConvectiveTerm().d_castTo<EulerTerm>())
{
}

//////////////////////////////////////////////////////////////////////////////

Euler2DKOmegaMachAlphaPTToCons::~Euler2DKOmegaMachAlphaPTToCons()
{
}

//////////////////////////////////////////////////////////////////////////////
void Euler2DKOmegaMachAlphaPTToCons::setup(const CFuint maxNbTransStates)
{
  VarSetTransformer::setup(maxNbTransStates);
}
//////////////////////////////////////////////////////////////////////////////
void Euler2DKOmegaMachAlphaPTToCons::transform(const State& state, State& result)
{
  const CFreal Mach = state[0];
  // convert angle in radiants
  const CFreal alpha = state[1]*MathTools::MathConsts::CFrealPi()/180.;
  const CFreal p = _model->getPressureFromState(state[2]);
  const CFreal T = state[3];
  const CFreal a = sqrt(_model->getGamma()*_model->getR()*T);
  const CFreal u = cos(alpha)*(Mach * a);
  const CFreal v = sin(alpha)*(Mach * a);
  const CFreal rho = _model->getDensity(p,T);
  const CFreal rhoV2 = rho*(u*u + v*v);
  const CFreal U = std::sqrt(u*u+v*v);  
  const CFreal K = 1e-6*U*U; 
  const CFreal Omega = 10.*U;
  
  result[0] = rho;
  result[1] = rho*u;
  result[2] = rho*v;
  result[3] = p/(_model->getGamma() - 1.) + 0.5*rhoV2;
  result[4] = K;
  result[5] = Omega;
}
//////////////////////////////////////////////////////////////////////////////
      
void Euler2DKOmegaMachAlphaPTToCons::transformFromRef(const RealVector& data,
						State& result)
{
  const CFreal rho = data[EulerTerm::RHO];
  
  result[0] = rho;
  result[1] = rho*data[EulerTerm::VX];
  result[2] = rho*data[EulerTerm::VY];
  result[3] = rho*data[EulerTerm::E];
  const CFuint iK = _model->getFirstScalarVar(0);
  result[4] = rho*data[iK];
  result[5] = rho*data[iK+1];
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NavierStokes

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
