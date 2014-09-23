#include "NavierStokes/NavierStokes.hh"
#include "Euler2DConsToRoe.hh"
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

Environment::ObjectProvider<Euler2DConsToRoe, VarSetTransformer, NavierStokesModule, 1> euler2DConsToRoeProvider("Euler2DConsToRoe");

//////////////////////////////////////////////////////////////////////////////

Euler2DConsToRoe::Euler2DConsToRoe(Common::SafePtr<Framework::PhysicalModelImpl> model) :
  VarSetTransformer(model),
 _model(model->getConvectiveTerm().d_castTo<EulerTerm>())
{
}

//////////////////////////////////////////////////////////////////////////////

Euler2DConsToRoe::~Euler2DConsToRoe()
{
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DConsToRoe::transform(const State& state, State& result)
{   
  const CFreal sqrho = sqrt(state[0]);
  result[0] = sqrho;
  result[1] = state[1]/sqrho;
  result[2] = state[2]/sqrho;
  
  const CFreal sqrhoU2 = result[1]*result[1]+ result[2]*result[2];
  const CFreal gamma = _model->getGamma();
  const CFreal gammaMinus1 = gamma - 1.;
  result[3] = (gamma*state[3] - 0.5*gammaMinus1*sqrhoU2)/sqrho;
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DConsToRoe::transformFromRef(const RealVector& data, State& result) 
{
  const CFreal sqRho = sqrt(data[EulerTerm::RHO]);
  
  result[0] = sqRho;
  result[1] = sqRho*data[EulerTerm::VX];
  result[2] = sqRho*data[EulerTerm::VY];
  result[3] = sqRho*data[EulerTerm::H];
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NavierStokes

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
