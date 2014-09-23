#include "NavierStokes/NavierStokes.hh"
#include "Euler2DRhoMachAlphaPToCons.hh"
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

Environment::ObjectProvider<Euler2DRhoMachAlphaPToCons, VarSetTransformer, NavierStokesModule, 1>
 euler2DRhoMachAlphaPToConsProvider("Euler2DRhoMachAlphaPToCons");

//////////////////////////////////////////////////////////////////////////////

Euler2DRhoMachAlphaPToCons::Euler2DRhoMachAlphaPToCons(Common::SafePtr<Framework::PhysicalModelImpl> model) :
  VarSetTransformer(model),
 _model(model->getConvectiveTerm().d_castTo<EulerTerm>())
{
}

//////////////////////////////////////////////////////////////////////////////

Euler2DRhoMachAlphaPToCons::~Euler2DRhoMachAlphaPToCons()
{
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DRhoMachAlphaPToCons::setup(const CFuint maxNbTransStates)
{
  VarSetTransformer::setup(maxNbTransStates);
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DRhoMachAlphaPToCons::transform(const State& state, State& result)
{
  const CFreal rho  = state[0];
  const CFreal Mach = state[1];
  // convert angle in radiants
  const CFreal alpha = state[2]*MathTools::MathConsts::CFrealPi()/180.;
  const CFreal p = _model->getPressureFromState(state[3]);
  const CFreal a = sqrt(_model->getGamma()*p/rho);
  const CFreal u = cos(alpha)*(Mach * a);
  const CFreal v = sin(alpha)*(Mach * a);
  const CFreal rhoV2 = rho*(u*u + v*v);
  
  result[0] = rho;
  result[1] = rho*u;
  result[2] = rho*v;
  result[3] = p/(_model->getGamma() - 1.) + 0.5*rhoV2;
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DRhoMachAlphaPToCons::transformFromRef(const RealVector& data, State& result)
{
  const CFreal rho = data[EulerTerm::RHO];

  result[0] = rho;
  result[1] = rho*data[EulerTerm::VX];
  result[2] = rho*data[EulerTerm::VY];
  result[3] = rho*data[EulerTerm::E];
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NavierStokes

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
