#include "NavierStokes/NavierStokes.hh"
#include "Euler2DConsToPuvt.hh"
#include "Framework/PhysicalModel.hh"
#include "Environment/ObjectProvider.hh"
#include "EulerTerm.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace NavierStokes {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<Euler2DConsToPuvt, VarSetTransformer, NavierStokesModule, 1> 
euler2DConsToPuvtProvider("Euler2DConsToPuvt");

//////////////////////////////////////////////////////////////////////////////

Euler2DConsToPuvt::Euler2DConsToPuvt(Common::SafePtr<Framework::PhysicalModelImpl> model) :
  VarSetTransformer(model),
  _model(model->getConvectiveTerm().d_castTo<EulerTerm>())
{
}

//////////////////////////////////////////////////////////////////////////////

Euler2DConsToPuvt::~Euler2DConsToPuvt()
{
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DConsToPuvt::transform(const State& state, State& result)
{
  const CFreal rho = state[0];
  const CFreal rhoE = state[3];
  const CFreal u = state[1]/rho;
  const CFreal v = state[2]/rho;
  const CFreal V2 = u*u + v*v;
  const CFreal ptot = (_model->getGamma() - 1.)*(rhoE - 0.5*rho*V2);
  const CFreal p = ptot - _model->getPressInf();
  
  result[0] = p;
  result[1] = u;
  result[2] = v;
  result[3] = ptot/(_model->getR()*rho);
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DConsToPuvt::transformFromRef(const RealVector& data, State& result)
{
  result[0] = data[EulerTerm::P];
  result[1] = data[EulerTerm::VX];
  result[2] = data[EulerTerm::VY];
  result[3] = data[EulerTerm::T];
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NavierStokes

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
