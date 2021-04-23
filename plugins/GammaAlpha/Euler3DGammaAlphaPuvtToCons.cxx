#include "GammaAlpha.hh"
#include "Euler3DGammaAlphaPuvtToCons.hh"
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

Environment::ObjectProvider<Euler3DGammaAlphaPuvtToCons, VarSetTransformer, GammaAlphaModule, 1>
euler3DGammaAlphaPuvtToConsProvider("Euler3DGammaAlphaPuvtToCons");

//////////////////////////////////////////////////////////////////////////////

Euler3DGammaAlphaPuvtToCons::Euler3DGammaAlphaPuvtToCons(Common::SafePtr<Framework::PhysicalModelImpl> model) :
  VarSetTransformer(model),
  _model(model->getConvectiveTerm().
	 d_castTo<EulerGammaAlphaTerm>())
{
}

//////////////////////////////////////////////////////////////////////////////

Euler3DGammaAlphaPuvtToCons::~Euler3DGammaAlphaPuvtToCons()
{
}

//////////////////////////////////////////////////////////////////////////////

void Euler3DGammaAlphaPuvtToCons::transform(const State& state, State& result)
{
  const CFreal R = _model->getR();
  const CFreal p = state[0];
  const CFreal u = state[1];
  const CFreal v = state[2];
  const CFreal w = state[3];
  const CFreal T = state[4];
  const CFreal K = state[5];
  const CFreal Omega = state[6];
  const CFreal Ga = state[7];
  const CFreal alpha = state[8];
  const CFreal rho = p/(R*T);
  const CFreal rhoV2 = rho*(u*u + v*v + w*w);

  result[0] = rho;
  result[1] = rho*u;
  result[2] = rho*v;
  result[3] = rho*w;
  result[4] = p/(_model->getGamma() - 1.) + 0.5*rhoV2;
  result[5] = rho*K;
  result[6] = rho*Omega;
  result[7] = rho*Ga;
  result[8] = rho*alpha;
}

//////////////////////////////////////////////////////////////////////////////

void Euler3DGammaAlphaPuvtToCons::transformFromRef(const RealVector& data, State& result)
{
  const CFreal rho = data[EulerGammaAlphaTerm::RHO];
  
  result[0] = rho;
  result[1] = rho*data[EulerGammaAlphaTerm::VX];
  result[2] = rho*data[EulerGammaAlphaTerm::VY];
  result[3] = rho*data[EulerGammaAlphaTerm::VZ];
  result[4] = rho*data[EulerGammaAlphaTerm::H] - data[EulerGammaAlphaTerm::P];
  
  const CFuint iK = _model->getFirstScalarVar(0);
  result[5] = rho*data[iK];
  result[6] = rho*data[iK+1];
  result[7] = rho*data[iK+2];
  result[8] = rho*data[iK+3];
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace GammaAlpha

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
