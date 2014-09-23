#include "NavierStokes/NavierStokes.hh"
#include "Euler2DRoeToConsInRef.hh"
#include "EulerPhysicalModel.hh"
#include "Framework/PhysicalModel.hh"
#include "Environment/ObjectProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace NavierStokes {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<Euler2DRoeToConsInRef, VarSetMatrixTransformer, NavierStokesModule, 1> euler2DRoeToConsInRefProvider("Euler2DRoeToConsInRef");

//////////////////////////////////////////////////////////////////////////////

Euler2DRoeToConsInRef::Euler2DRoeToConsInRef(Common::SafePtr<Framework::PhysicalModelImpl> model) :
  VarSetMatrixTransformer(model),
  _model(model->getConvectiveTerm().d_castTo<EulerTerm>())
{
}

//////////////////////////////////////////////////////////////////////////////

Euler2DRoeToConsInRef::~Euler2DRoeToConsInRef()
{
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DRoeToConsInRef::setMatrixFromRef()
{
  cf_assert(_model.isNotNull());
  const RealVector& linearData = _model->getPhysicalData();

  const CFreal z0 = sqrt(linearData[EulerTerm::RHO]);
  const CFreal z1 = z0*linearData[EulerTerm::VX];
  const CFreal z2 = z0*linearData[EulerTerm::VY];
  const CFreal z3 = z0*linearData[EulerTerm::H];
  const CFreal gamma = _model->getGamma();
  const CFreal gammaMinus1 = gamma - 1.;
  const CFreal gammaMinus1DivGamma = gammaMinus1/gamma;

  _transMatrix(0,0) = 2.0*z0;
  _transMatrix(1,0) = z1;
  _transMatrix(1,1) = z0;
  _transMatrix(2,0) = z2;
  _transMatrix(2,2) = z0;
  _transMatrix(3,0) = z3/gamma;
  _transMatrix(3,1) = z1*gammaMinus1DivGamma;
  _transMatrix(3,2) = z2*gammaMinus1DivGamma;
  _transMatrix(3,3) = z0/gamma;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NavierStokes

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
