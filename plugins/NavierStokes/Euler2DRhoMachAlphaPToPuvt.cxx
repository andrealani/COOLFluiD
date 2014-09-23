#include "NavierStokes/NavierStokes.hh"
#include "Euler2DRhoMachAlphaPToPuvt.hh"
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

Environment::ObjectProvider<Euler2DRhoMachAlphaPToPuvt, VarSetTransformer, NavierStokesModule, 1>
 euler2DRhoMachAlphaPToPuvtProvider("Euler2DRhoMachAlphaPToPuvt");

//////////////////////////////////////////////////////////////////////////////

Euler2DRhoMachAlphaPToPuvt::Euler2DRhoMachAlphaPToPuvt(Common::SafePtr<Framework::PhysicalModelImpl> model) :
  VarSetTransformer(model),
 _model(model->getConvectiveTerm().d_castTo<EulerTerm>())
{
}

//////////////////////////////////////////////////////////////////////////////

Euler2DRhoMachAlphaPToPuvt::~Euler2DRhoMachAlphaPToPuvt()
{
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DRhoMachAlphaPToPuvt::setup(const CFuint maxNbTransStates)
{
  VarSetTransformer::setup(maxNbTransStates);
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DRhoMachAlphaPToPuvt::transform(const State& state, State& result)
{
  const CFreal rho  = state[0];
  const CFreal Mach = state[1];
  
  const CFreal p = _model->getPressureFromState(state[3]);
  
  const CFreal a = sqrt(_model->getGamma()*p/rho);
  // convert angle in radiants
  const CFreal alpha = state[2]*MathTools::MathConsts::CFrealPi()/180.;
  const CFreal u = cos(alpha)*(Mach * a);
  const CFreal v = sin(alpha)*(Mach * a);
  
  const CFreal Rg = _model->getR();  
  const CFreal T  = p/(Rg*rho);
  
  result[0] = p;
  result[1] = u;
  result[2] = v;
  result[3] = T;
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DRhoMachAlphaPToPuvt::transformFromRef(const RealVector& data, State& result)
{
  result[0] =  data[EulerTerm::P];
  result[1] =  data[EulerTerm::VX];
  result[2] =  data[EulerTerm::VY];
  result[3] =  data[EulerTerm::T];
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NavierStokes

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
