#include "NavierStokes/NavierStokes.hh"
#include "Euler3DConsToPvtInPvt.hh"
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

Environment::ObjectProvider<Euler3DConsToPvtInPvt, VarSetMatrixTransformer, NavierStokesModule, 1> euler3DConsToPvtInPvtProvider("Euler3DConsToPvtInPvt");

//////////////////////////////////////////////////////////////////////////////

Euler3DConsToPvtInPvt::Euler3DConsToPvtInPvt(Common::SafePtr<Framework::PhysicalModelImpl> model) :
  VarSetMatrixTransformer(model),
 _model(model->getConvectiveTerm().d_castTo<EulerTerm>())
{
}

//////////////////////////////////////////////////////////////////////////////

Euler3DConsToPvtInPvt::~Euler3DConsToPvtInPvt()
{
}

//////////////////////////////////////////////////////////////////////////////

void Euler3DConsToPvtInPvt::setMatrix(const RealVector& state)
{
  // not consistent for incompressible
  const CFreal R = _model->getR();
  const CFreal p = state[0];
  const CFreal invRho = R*state[3]/p;
  const CFreal u = state[1];
  const CFreal v = state[2];
  const CFreal w = state[3];
  const CFreal gammaMinus1 = _model->getGamma() - 1.;

  _transMatrix(0,0) = (u*u + v*v + w*w)*0.5*gammaMinus1;
  _transMatrix(0,1) = -gammaMinus1*u;
  _transMatrix(0,2) = -gammaMinus1*v;
  _transMatrix(0,3) = -gammaMinus1*w;
  _transMatrix(0,4) = gammaMinus1;

  _transMatrix(1,0) = -u*invRho;
  _transMatrix(1,1) = invRho;

  _transMatrix(2,0) = -v*invRho;
  _transMatrix(2,2) = invRho;

  _transMatrix(3,0) = -w*invRho;
  _transMatrix(3,3) = invRho;

  _transMatrix(4,0) = (_transMatrix(0,0) - p*invRho)*invRho/R;
  _transMatrix(4,1) = _transMatrix(0,1)*invRho/R;
  _transMatrix(4,2) = _transMatrix(0,2)*invRho/R;
  _transMatrix(4,3) = _transMatrix(0,3)*invRho/R;
  _transMatrix(4,4) = _transMatrix(0,4)*invRho/R;

}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NavierStokes

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
