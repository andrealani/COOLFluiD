#include "NavierStokes/NavierStokes.hh"
#include "Euler3DMachAlphaPTToCons.hh"
#include "Framework/PhysicalModel.hh"
#include "Environment/ObjectProvider.hh"
#include "EulerPhysicalModel.hh"
#include "MathTools/MathConsts.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace NavierStokes {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<Euler3DMachAlphaPTToCons, VarSetTransformer, NavierStokesModule, 1> 
euler3DMachAlphaPTToConsProvider("Euler3DMachAlphaPTToCons");

//////////////////////////////////////////////////////////////////////////////

Euler3DMachAlphaPTToCons::Euler3DMachAlphaPTToCons
(Common::SafePtr<Framework::PhysicalModelImpl> model) :
  VarSetTransformer(model),
 _model(model->getConvectiveTerm().d_castTo<EulerTerm>())
{
}

//////////////////////////////////////////////////////////////////////////////

Euler3DMachAlphaPTToCons::~Euler3DMachAlphaPTToCons()
{
}

//////////////////////////////////////////////////////////////////////////////

void Euler3DMachAlphaPTToCons::setup(const CFuint maxNbTransStates)
{
  VarSetTransformer::setup(maxNbTransStates);
}

//////////////////////////////////////////////////////////////////////////////

void Euler3DMachAlphaPTToCons::transform(const State& state,
					 State& result)
{
  const CFreal Mach = state[0];
  // convert angle in radiants
  const CFreal alpha = state[1]*MathTools::MathConsts::CFrealPi()/180.;
  const CFreal beta = state[2]*MathTools::MathConsts::CFrealPi()/180.;
  const CFreal p = _model->getPressureFromState(state[3]);
  const CFreal T = state[4];
  const CFreal a = sqrt(_model->getGamma()*_model->getR()*T);
  const CFreal tgA = tan(alpha);
  const CFreal tgB = tan(beta);
  const CFreal u = (Mach*a)/sqrt(1. + tgA*tgA + tgB*tgB);
  const CFreal v = u*tan(alpha);
  const CFreal w = u*tan(beta);
  const CFreal rho = _model->getDensity(p,T);
  const CFreal rhoV2 = rho*(u*u + v*v + w*w);
  
  result[0] = rho;
  result[1] = rho*u;
  result[2] = rho*v;
  result[3] = rho*w;
  result[4] = p/(_model->getGamma() - 1.) + 0.5*rhoV2;
}

//////////////////////////////////////////////////////////////////////////////

void Euler3DMachAlphaPTToCons::transformFromRef(const RealVector& data,
						State& result)
{
  const CFreal rho = data[EulerTerm::RHO];

  result[0] = rho;
  result[1] = rho*data[EulerTerm::VX];
  result[2] = rho*data[EulerTerm::VY];
  result[3] = rho*data[EulerTerm::VZ];
  result[4] = rho*data[EulerTerm::E];
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NavierStokes

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
