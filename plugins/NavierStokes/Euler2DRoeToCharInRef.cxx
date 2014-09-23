#include "NavierStokes/NavierStokes.hh"
#include "Euler2DRoeToCharInRef.hh"
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

Environment::ObjectProvider<Euler2DRoeToCharInRef, VarSetMatrixTransformer, NavierStokesModule, 1> euler2DRoeToCharInRefProvider("Euler2DRoeToCharInRef");

//////////////////////////////////////////////////////////////////////////////

Euler2DRoeToCharInRef::Euler2DRoeToCharInRef(Common::SafePtr<Framework::PhysicalModelImpl> model) :
  VarSetMatrixTransformer(model),
  _model(model->getConvectiveTerm().d_castTo<EulerTerm>())
{
}

//////////////////////////////////////////////////////////////////////////////

Euler2DRoeToCharInRef::~Euler2DRoeToCharInRef()
{
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DRoeToCharInRef::setMatrixFromRef()
{
  cf_assert(_model.isNotNull());

  const RealVector& linearData = _model->getPhysicalData();

  const CFreal z0 = sqrt(linearData[EulerTerm::RHO]);
  const CFreal z1 = z0*linearData[EulerTerm::VX];
  const CFreal z2 = z0*linearData[EulerTerm::VY];
  const CFreal z3 = z0*linearData[EulerTerm::H];
  const CFreal invAvZ0 = 1./ z0;
  const CFreal invAvZ02 = invAvZ0*invAvZ0;
  const CFreal avRho = z0* z0;
  const CFreal avU   = z1*invAvZ0;
  const CFreal avV   = z2*invAvZ0;
  const CFreal cosD = avU/linearData[EulerTerm::V];
  const CFreal sinD = avV/linearData[EulerTerm::V];
  const CFreal avA = linearData[EulerTerm::A];
  const CFreal M     = linearData[EulerTerm::V]/avA;
  const CFreal eps   = 0.05;
  const CFreal beta  = sqrt(max(eps*eps,std::abs(M*M - 1.)));
  const CFreal McosD = M*cosD;
  const CFreal MsinD = M*sinD;
  const CFreal invARho  = 1./(avRho*avA);
  const CFreal gamma = _model->getGamma();
  const CFreal gammaMinus1 = gamma - 1.;
  const CFreal gammaMinus1DivGamma = gammaMinus1/gamma;
  const CFreal k0 = invARho*gammaMinus1DivGamma;
  const CFreal k1 = beta*k0;

  _transMatrix(0,0) = M*invAvZ02*(sinD*z1 - cosD*z2)
    + k1*z3;
  _transMatrix(0,1) = -MsinD*invAvZ0 - k1*z1;
  _transMatrix(0,2) = McosD*invAvZ0 - k1*z2;
  _transMatrix(0,3) = k1*z0;

  _transMatrix(1,0) = M*invAvZ02*(-sinD*z1 + cosD*z2)
    + k1*z3;
  _transMatrix(1,1) = MsinD*invAvZ0 - k1*z1;
  _transMatrix(1,2) = -McosD*invAvZ0 - k1*z2;
  _transMatrix(1,3) = k1*z0;

  _transMatrix(2,0) = M*invAvZ02*(-cosD*z1 - sinD*z2)
    + k0*z3;
  _transMatrix(2,1) = McosD*invAvZ0 - k0*z1;
  _transMatrix(2,2) = MsinD*invAvZ0 - k0*z2;
  _transMatrix(2,3) = k0*z0;

  _transMatrix(3,0) = -2.0*z0*avA*avA +
    gammaMinus1DivGamma*z3;
  _transMatrix(3,1) = -gammaMinus1DivGamma*z1;
  _transMatrix(3,2) = -gammaMinus1DivGamma*z2;
  _transMatrix(3,3) = gammaMinus1DivGamma*z0;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NavierStokes

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
