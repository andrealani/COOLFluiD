#include "NavierStokes/NavierStokes.hh"
#include "Euler2DGReKOMachAlphaPTToPuvt.hh"
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

Environment::ObjectProvider<Euler2DGReKOMachAlphaPTToPuvt, VarSetTransformer, NavierStokesModule, 1> 
eulerGReKO2DMachAlphaPTToPuvtProvider("Euler2DGReKOMachAlphaPTToPuvt");

//////////////////////////////////////////////////////////////////////////////

Euler2DGReKOMachAlphaPTToPuvt::Euler2DGReKOMachAlphaPTToPuvt(Common::SafePtr<Framework::PhysicalModelImpl> model) :
  VarSetTransformer(model),
 _model(model->getConvectiveTerm().d_castTo<EulerGReKOTerm>())
{
}

//////////////////////////////////////////////////////////////////////////////

Euler2DGReKOMachAlphaPTToPuvt::~Euler2DGReKOMachAlphaPTToPuvt()
{
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DGReKOMachAlphaPTToPuvt::setup(const CFuint maxNbTransStates)
{
  VarSetTransformer::setup(maxNbTransStates);
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DGReKOMachAlphaPTToPuvt::transform(const State& state, State& result)
{
  const CFreal Mach = state[0];
  // convert angle in radiants
  const CFreal alpha = state[1]*MathTools::MathConsts::CFrealPi()/180.;
  const CFreal p = state[2];
  const CFreal T = state[3];
  const CFreal a = sqrt(_model->getGamma()*_model->getR()*T);
  const CFreal u = cos(alpha)*(Mach * a);
  const CFreal v = sin(alpha)*(Mach * a);
   
  result[0] = p;
  result[1] = u;
  result[2] = v;
  result[3] = T;
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DGReKOMachAlphaPTToPuvt::transformFromRef(const RealVector& data,State& result)
{
  const CFreal p = data[EulerTerm::P];
  const CFreal T = p/(data[EulerTerm::RHO]*_model->getR());
  
  result[0] = p;
  result[1] = data[EulerTerm::VX];
  result[2] = data[EulerTerm::VY];
  result[3] = T;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NavierStokes

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
