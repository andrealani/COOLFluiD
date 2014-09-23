#include "NavierStokes/NavierStokes.hh"
#include "Euler2DConsToPrimInPrim.hh"
#include "EulerPhysicalModel.hh"
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

Environment::ObjectProvider<Euler2DConsToPrimInPrim, VarSetMatrixTransformer, NavierStokesModule, 1> euler2DConsToPrimInPrimProvider("Euler2DConsToPrimInPrim");

//////////////////////////////////////////////////////////////////////////////

Euler2DConsToPrimInPrim::Euler2DConsToPrimInPrim(Common::SafePtr<Framework::PhysicalModelImpl> model) :
  VarSetMatrixTransformer(model),
  _model(model->getConvectiveTerm().d_castTo<EulerTerm>())
{
}

//////////////////////////////////////////////////////////////////////////////

Euler2DConsToPrimInPrim::~Euler2DConsToPrimInPrim()
{
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DConsToPrimInPrim::setMatrix(const RealVector& state)
{
  const CFreal invRho = 1./state[0];
  const CFreal u = state[1];
  const CFreal v = state[2];
  const CFreal gammaMinus1 = _model->getGamma() - 1.;

  _transMatrix(0,0) = 1.0;
  _transMatrix(1,0) = -u*invRho;
  _transMatrix(1,1) = invRho;
  _transMatrix(2,0) = -v*invRho;
  _transMatrix(2,2) = invRho;
  _transMatrix(3,0) = (u*u + v*v)*0.5*gammaMinus1;
  _transMatrix(3,1) = -gammaMinus1*u;
  _transMatrix(3,2) = -gammaMinus1*v;
  _transMatrix(3,3) = gammaMinus1;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NavierStokes

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
