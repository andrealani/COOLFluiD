#include "NavierStokes/NavierStokes.hh"
#include "Euler2DRoeToSymmInRef.hh"
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

Environment::ObjectProvider<Euler2DRoeToSymmInRef, VarSetMatrixTransformer, NavierStokesModule, 1> euler2DRoeToSymmInRefProvider("Euler2DRoeToSymmInRef");

//////////////////////////////////////////////////////////////////////////////

Euler2DRoeToSymmInRef::Euler2DRoeToSymmInRef(Common::SafePtr<Framework::PhysicalModelImpl> model) :
  VarSetMatrixTransformer(model),
  _model(model->getConvectiveTerm().d_castTo<EulerTerm>())
{
}

//////////////////////////////////////////////////////////////////////////////

Euler2DRoeToSymmInRef::~Euler2DRoeToSymmInRef()
{
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DRoeToSymmInRef::setMatrixFromRef()
{
  cf_assert(_model.isNotNull());

  const RealVector& linearData = _model->getPhysicalData();

  const CFreal z0 = sqrt(linearData[EulerTerm::RHO]);
  const CFreal z1 = z0*linearData[EulerTerm::VX];
  const CFreal z2 = z0*linearData[EulerTerm::VY];
  const CFreal z3 = z0*linearData[EulerTerm::H];
  const CFreal invAvZ0 = 1./z0;
  const CFreal invAvZ02 = invAvZ0*invAvZ0;
  const CFreal avRho = z0*z0;
  const CFreal avU   = z1*invAvZ0;
  const CFreal avV   = z2*invAvZ0;
  const CFreal speed = sqrt(avU*avU + avV*avV);
  const CFreal cosD = avU/speed;
  const CFreal sinD = avV/speed;
  const CFreal avA = linearData[EulerTerm::A];
  const CFreal gamma = _model->getGamma();
  const CFreal gammaMinus1 = gamma - 1.;
  const CFreal gammaMinus1DivGamma = gammaMinus1/gamma;
  const CFreal k1 = gammaMinus1DivGamma/(avA*avRho);

  _transMatrix(0,0) = k1*z3;
  _transMatrix(0,1) = -k1*z1;
  _transMatrix(0,2) = -k1*z2;
  _transMatrix(0,3) = k1*z0;
  _transMatrix(1,0) = -(z1*cosD + z2*sinD)*invAvZ02;
  _transMatrix(1,1) = cosD*invAvZ0;
  _transMatrix(1,2) = sinD*invAvZ0;
  _transMatrix(2,0) = (z1*sinD - z2*cosD)*invAvZ02;
  _transMatrix(2,1) = -sinD*invAvZ0;
  _transMatrix(2,2) = cosD*invAvZ0;
  _transMatrix(3,0) = -2.*z0*avA*avA + gammaMinus1DivGamma*
    z3;
  _transMatrix(3,1) = -gammaMinus1DivGamma*z1;
  _transMatrix(3,2) = -gammaMinus1DivGamma*z2;
  _transMatrix(3,3) = gammaMinus1DivGamma*z0;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NavierStokes

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
