#include "Euler3DSAPvtToCons.hh"
#include "Framework/PhysicalModel.hh"
#include "Environment/ObjectProvider.hh"
#include "SA.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace SA {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<Euler3DSAPvtToCons, VarSetTransformer, SAModule, 1> 
euler3DSAPvtToConsProvider("Euler3DSAPvtToCons");

//////////////////////////////////////////////////////////////////////////////

Euler3DSAPvtToCons::Euler3DSAPvtToCons
(Common::SafePtr<Framework::PhysicalModelImpl> model) :
  VarSetTransformer(model),
  _model(model->getConvectiveTerm().
	 d_castTo<MultiScalarTerm<EulerTerm> >())
{
}

//////////////////////////////////////////////////////////////////////////////

Euler3DSAPvtToCons::~Euler3DSAPvtToCons()
{
}

//////////////////////////////////////////////////////////////////////////////

void Euler3DSAPvtToCons::transform(const State& state, State& result)
{

  const CFreal R = _model->getR();
  const CFreal p = state[0];
  const CFreal u = state[1];
  const CFreal v = state[2];
  const CFreal w = state[3];
  const CFreal T = state[4];
  const CFreal Nutil = state[5];
  const CFreal rho = p/(R*T);
  const CFreal rhoV2 = rho*(u*u + v*v + w*w);

  result[0] = rho;
  result[1] = rho*u;
  result[2] = rho*v;
  result[3] = rho*w;
  result[4] = p/(_model->getGamma() - 1.) + 0.5*rhoV2; //rhoE
  result[5] = rho*Nutil;
}

//////////////////////////////////////////////////////////////////////////////

void Euler3DSAPvtToCons::transformFromRef(const RealVector& data, State& result)
{
  const CFreal rho = data[EulerSATerm::RHO];

  result[0] = rho;
  result[1] = rho*data[EulerSATerm::VX];
  result[2] = rho*data[EulerSATerm::VY];
  result[3] = rho*data[EulerSATerm::VZ];
  result[4] = rho*data[EulerSATerm::H] - data[EulerSATerm::P]; //rhoE
  
  const CFuint iNutil = _model->getFirstScalarVar(0);
  result[5] = rho*data[iNutil];
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace SA

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
