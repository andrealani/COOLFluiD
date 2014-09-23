#include "NavierStokes/NavierStokes.hh"
#include "Euler2DSymmToConsInRef.hh"
#include "EulerPhysicalModel.hh"
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

Environment::ObjectProvider<Euler2DSymmToConsInRef, VarSetMatrixTransformer, NavierStokesModule, 1> euler2DSymmToConsInRefProvider("Euler2DSymmToConsInRef");

//////////////////////////////////////////////////////////////////////////////

Euler2DSymmToConsInRef::Euler2DSymmToConsInRef(Common::SafePtr<Framework::PhysicalModelImpl> model) :
  VarSetMatrixTransformer(model),
 _model(model->getConvectiveTerm().d_castTo<EulerTerm>())
{
}

//////////////////////////////////////////////////////////////////////////////

Euler2DSymmToConsInRef::~Euler2DSymmToConsInRef()
{
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DSymmToConsInRef::setMatrixFromRef()
{
  cf_assert(_model.isNotNull());

  const RealVector& linearData = _model->getPhysicalData();

  const CFreal avRho = linearData[EulerTerm::RHO];
  const CFreal avU = linearData[EulerTerm::VX];
  const CFreal avV = linearData[EulerTerm::VY];
  const CFreal speed = linearData[EulerTerm::V];
  const CFreal cosD = avU/speed;
  const CFreal sinD = avV/speed;
  const CFreal avA   = linearData[EulerTerm::A];
  const CFreal rhoDivA = avRho/avA;
  const CFreal invA2   = 1./(avA*avA);
  const CFreal invGammaMinus1 = 1./(_model->getGamma() - 1.);

  _transMatrix(0,0) = rhoDivA;
  _transMatrix(0,3) = -invA2;
  _transMatrix(1,0) = rhoDivA*avU;
  _transMatrix(1,1) = avRho*cosD;
  _transMatrix(1,2) = -avRho*sinD;
  _transMatrix(1,3) = -invA2*avU;
  _transMatrix(2,0) = rhoDivA*avV;
  _transMatrix(2,1) = avRho*sinD;
  _transMatrix(2,2) = avRho*cosD;
  _transMatrix(2,3) = -invA2*avV;
  _transMatrix(3,0) = rhoDivA*speed*speed*0.5 + avRho*avA*invGammaMinus1;
  _transMatrix(3,1) = avRho*avU*cosD + avRho*avV*sinD;
  _transMatrix(3,2) = -avRho*avU*sinD + avRho*avV*cosD;
  _transMatrix(3,3) = -0.5*speed*speed*invA2;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NavierStokes

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
