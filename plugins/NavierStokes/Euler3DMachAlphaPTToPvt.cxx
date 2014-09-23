#include "NavierStokes/NavierStokes.hh"
#include "Euler3DMachAlphaPTToPvt.hh"
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

Environment::ObjectProvider<Euler3DMachAlphaPTToPvt, VarSetTransformer, NavierStokesModule, 1> 
euler3DMachAlphaPTToPvtProvider("Euler3DMachAlphaPTToPvt");

//////////////////////////////////////////////////////////////////////////////

Euler3DMachAlphaPTToPvt::Euler3DMachAlphaPTToPvt
(Common::SafePtr<Framework::PhysicalModelImpl> model) :
  VarSetTransformer(model),
 _model(model->getConvectiveTerm().d_castTo<EulerTerm>())
{
}

//////////////////////////////////////////////////////////////////////////////

Euler3DMachAlphaPTToPvt::~Euler3DMachAlphaPTToPvt()
{
}

//////////////////////////////////////////////////////////////////////////////

void Euler3DMachAlphaPTToPvt::setup(const CFuint maxNbTransStates)
{
  VarSetTransformer::setup(maxNbTransStates);
}

//////////////////////////////////////////////////////////////////////////////

void Euler3DMachAlphaPTToPvt::transform(const State& state, State& result)
{
  const CFreal Mach = state[0];
  // convert angle in radiants
  const CFreal alpha = state[1]*MathTools::MathConsts::CFrealPi()/180.;
  const CFreal beta = state[2]*MathTools::MathConsts::CFrealPi()/180.;
  const CFreal p = state[3];
  const CFreal T = state[4];
  const CFreal a = sqrt(_model->getGamma()*_model->getR()*T);
  const CFreal tgA = tan(alpha);
  const CFreal tgB = tan(beta);
  const CFreal u = (Mach*a)/sqrt(1. + tgA*tgA + tgB*tgB);
  const CFreal v = u*tan(alpha);
  const CFreal w = u*tan(beta);
  
  result[0] = p;
  result[1] = u;
  result[2] = v;
  result[3] = w;
  result[4] = T;
}

//////////////////////////////////////////////////////////////////////////////

void Euler3DMachAlphaPTToPvt::transformFromRef(const RealVector& data, State& result)
{
  const CFreal p = data[EulerTerm::P];
  const CFreal T = p/(data[EulerTerm::RHO]*_model->getR());
  
  result[0] = p;
  result[1] = data[EulerTerm::VX];
  result[2] = data[EulerTerm::VY];
  result[3] = data[EulerTerm::VZ];
  result[4] = T;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NavierStokes

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
