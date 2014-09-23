#include "KOmega.hh"
#include "Euler2DKOmegaMachAlphaPTToPuvt.hh"
#include "Framework/PhysicalModel.hh"
#include "Environment/ObjectProvider.hh"
#include "MathTools/MathFunctions.hh"
#include "MathTools/MathConsts.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace KOmega {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<Euler2DKOmegaMachAlphaPTToPuvt, VarSetTransformer,KOmegaModule, 1> 
euler2DKOmegaMachAlphaPTToPuvtProvider("Euler2DKOmegaMachAlphaPTToPuvt");

//////////////////////////////////////////////////////////////////////////////

Euler2DKOmegaMachAlphaPTToPuvt::Euler2DKOmegaMachAlphaPTToPuvt(Common::SafePtr<Framework::PhysicalModelImpl> model) :
  VarSetTransformer(model),
 _model(model->getConvectiveTerm().d_castTo<EulerKOmegaTerm>())
{
}

//////////////////////////////////////////////////////////////////////////////

Euler2DKOmegaMachAlphaPTToPuvt::~Euler2DKOmegaMachAlphaPTToPuvt()
{
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DKOmegaMachAlphaPTToPuvt::setup(const CFuint maxNbTransStates)
{
  VarSetTransformer::setup(maxNbTransStates);
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DKOmegaMachAlphaPTToPuvt::transform(const State& state, State& result)
{
  const CFreal Mach = state[0];
  // convert angle in radiants
  const CFreal alpha = state[1]*MathTools::MathConsts::CFrealPi()/180.;
  const CFreal p = state[2];
  const CFreal T = state[3];
  const CFreal a = sqrt(_model->getGamma()*_model->getR()*T);
  const CFreal u = cos(alpha)*(Mach * a);
  const CFreal v = sin(alpha)*(Mach * a);
  
  const CFreal U = std::sqrt(u*u+v*v);  
  const CFreal K = 1e-6*U*U; 
  const CFreal Omega = 10.*U;
  result[0] = p;
  result[1] = u;
  result[2] = v;
  result[3] = T;
  result[4] = K;
  result[5] = Omega;
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DKOmegaMachAlphaPTToPuvt::transformFromRef(const RealVector& data,State& result)
{
  const CFreal p = _model->getPressureFromState(data[EulerTerm::P]);
  const CFreal T = p/(data[EulerTerm::RHO]*_model->getR());
  
  result[0] = data[EulerTerm::P];
  result[1] = data[EulerTerm::VX];
  result[2] = data[EulerTerm::VY];
  result[3] = T;
  const CFuint iK = _model->getFirstScalarVar(0);
  result[4] = data[iK] ;
  result[5] = data[iK+1];
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NavierStokes

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
