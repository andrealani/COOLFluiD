#include "GammaAlpha.hh"
#include "Euler2DGammaAlphaPuvtToCons.hh"
#include "Framework/PhysicalModel.hh"
#include "Environment/ObjectProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace GammaAlpha {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<Euler2DGammaAlphaPuvtToCons, VarSetTransformer, GammaAlphaModule, 1>
euler2DGammaAlphaPuvtToConsProvider("Euler2DGammaAlphaPuvtToCons");

//////////////////////////////////////////////////////////////////////////////

Euler2DGammaAlphaPuvtToCons::Euler2DGammaAlphaPuvtToCons(Common::SafePtr<Framework::PhysicalModelImpl> model) :
  VarSetTransformer(model),
  _model(model->getConvectiveTerm().
	 d_castTo<EulerGammaAlphaTerm>())
{
}

//////////////////////////////////////////////////////////////////////////////

Euler2DGammaAlphaPuvtToCons::~Euler2DGammaAlphaPuvtToCons()
{
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DGammaAlphaPuvtToCons::transform(const State& state, State& result)
{
  const CFreal R = _model->getR();
  const CFreal p = state[0];
  const CFreal u = state[1];
  const CFreal v = state[2];
  const CFreal T = state[3];
  const CFreal K = state[4];
  const CFreal Omega = state[5];
  const CFreal Ga = state[6];
  const CFreal alpha = state[7];
  const CFreal rho = p/(R*T);
  const CFreal rhoV2 = rho*(u*u + v*v);

  result[0] = rho;
  result[1] = rho*u;
  result[2] = rho*v;
  result[3] = p/(_model->getGamma() - 1.) + 0.5*rhoV2;
  result[4] = rho*K;
  result[5] = rho*Omega;
  result[6] = rho*Ga;
  result[7] = rho*alpha;
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DGammaAlphaPuvtToCons::transformFromRef(const RealVector& data, State& result)
{
  const CFreal rho = data[EulerGammaAlphaTerm::RHO];
  
  result[0] = rho;
  result[1] = rho*data[EulerGammaAlphaTerm::VX];
  result[2] = rho*data[EulerGammaAlphaTerm::VY];
  result[3] = rho*data[EulerGammaAlphaTerm::H] - data[EulerGammaAlphaTerm::P];
  
  const CFuint iK = _model->getFirstScalarVar(0);
  result[4] = rho*data[iK];
  result[5] = rho*data[iK+1];
  result[6] = rho*data[iK+2];
  result[7] = rho*data[iK+3];
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace GammaAlpha

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
