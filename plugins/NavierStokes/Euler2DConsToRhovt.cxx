#include "NavierStokes/NavierStokes.hh"
#include "Euler2DConsToRhovt.hh"
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

Environment::ObjectProvider<Euler2DConsToRhovt, VarSetTransformer, NavierStokesModule, 1> 
euler2DConsToRhovtProvider("Euler2DConsToRhovt");

//////////////////////////////////////////////////////////////////////////////

Euler2DConsToRhovt::Euler2DConsToRhovt(Common::SafePtr<Framework::PhysicalModelImpl> model) :
  VarSetTransformer(model),
  _model(model->getConvectiveTerm().d_castTo<EulerTerm>())
{
}

//////////////////////////////////////////////////////////////////////////////

Euler2DConsToRhovt::~Euler2DConsToRhovt()
{
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DConsToRhovt::transform(const State& state, State& result)
{
  const CFreal rho = state[0];
  const CFreal rhoE = state[3];
  const CFreal u = state[1]/rho;
  const CFreal v = state[2]/rho;
  const CFreal V2 = u*u + v*v;
  const CFreal p = (_model->getGamma() - 1.)*(rhoE - 0.5*rho*V2);
  
  result[0] = rho;
  result[1] = u;
  result[2] = v;
  result[3] = p/(_model->getR()*rho);
}
      
//////////////////////////////////////////////////////////////////////////////

void Euler2DConsToRhovt::transformFromRef(const RealVector& data, State& result)
{  
  result[0] = data[EulerTerm::RHO];
  result[1] = data[EulerTerm::VX];
  result[2] = data[EulerTerm::VY];
  result[3] = data[EulerTerm::T];
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NavierStokes

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
