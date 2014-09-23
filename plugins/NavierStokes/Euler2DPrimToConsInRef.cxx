#include "NavierStokes/NavierStokes.hh"
#include "Euler2DPrimToConsInRef.hh"
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

Environment::ObjectProvider<Euler2DPrimToConsInRef, VarSetMatrixTransformer, NavierStokesModule, 1> euler2DPrimToConsInRefProvider("Euler2DPrimToConsInRef");

//////////////////////////////////////////////////////////////////////////////

Euler2DPrimToConsInRef::Euler2DPrimToConsInRef(Common::SafePtr<Framework::PhysicalModelImpl> model) :
  VarSetMatrixTransformer(model),
  _model(model->getConvectiveTerm().d_castTo<EulerTerm>())
{
}

//////////////////////////////////////////////////////////////////////////////

Euler2DPrimToConsInRef::~Euler2DPrimToConsInRef()
{
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DPrimToConsInRef::setMatrixFromRef()
{
  const RealVector& linearData = _model->getPhysicalData();

  const CFreal avRho = linearData[EulerTerm::RHO];
  const CFreal avU = linearData[EulerTerm::VX];
  const CFreal avV = linearData[EulerTerm::VY];

  _transMatrix(0,0) = 1.0;
  _transMatrix(1,0) = avU;
  _transMatrix(1,1) = avRho;
  _transMatrix(2,0) = avV;
  _transMatrix(2,2) = avRho;
  _transMatrix(3,0) = (avU*avU + avV*avV)*0.5;
  _transMatrix(3,1) = avRho*avU;
  _transMatrix(3,2) = avRho*avV;
  _transMatrix(3,3) = 1./(_model->getGamma() - 1.);
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NavierStokes

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
