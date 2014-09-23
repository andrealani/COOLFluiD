#include "NavierStokes/NavierStokes.hh"
#include "Euler1DConsToRoe.hh"
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

Environment::ObjectProvider<Euler1DConsToRoe, VarSetTransformer, NavierStokesModule, 1> euler1DConsToRoeProvider("Euler1DConsToRoe");

//////////////////////////////////////////////////////////////////////////////

Euler1DConsToRoe::Euler1DConsToRoe(Common::SafePtr<Framework::PhysicalModelImpl> model) :
  VarSetTransformer(model),
 _model(model->getConvectiveTerm().d_castTo<EulerTerm>())
{
}

//////////////////////////////////////////////////////////////////////////////

Euler1DConsToRoe::~Euler1DConsToRoe()
{
}

//////////////////////////////////////////////////////////////////////////////

void Euler1DConsToRoe::transform(const State& state, State& result)
{ 
  const CFreal sqrho = sqrt(state[0]);
  result[0] = sqrho;
  result[1] = state[1]/sqrho;
  
  const CFreal sqrhoU2 = result[1]*result[1];
  const CFreal gamma = _model->getGamma();
  const CFreal gammaMinus1 = gamma - 1.;
  result[2] = (gamma*state[2] - 0.5*gammaMinus1*sqrhoU2)/sqrho;
  
}

//////////////////////////////////////////////////////////////////////////////

void Euler1DConsToRoe::transformFromRef(const RealVector& data, State& result) 
{
  const CFreal sqRho = sqrt(data[EulerTerm::RHO]);
  
  result[0] = sqRho;
  result[1] = sqRho*data[EulerTerm::VX];
  result[2] = sqRho*data[EulerTerm::H];
  
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NavierStokes

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
