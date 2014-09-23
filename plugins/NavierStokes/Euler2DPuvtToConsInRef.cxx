#include "NavierStokes/NavierStokes.hh"
#include "Euler2DPuvtToConsInRef.hh"
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

Environment::ObjectProvider<Euler2DPuvtToConsInRef, VarSetMatrixTransformer,
NavierStokesModule, 1> euler2DPuvtToConsInRefProvider("Euler2DPuvtToConsInRef");

//////////////////////////////////////////////////////////////////////////////

Euler2DPuvtToConsInRef::Euler2DPuvtToConsInRef(Common::SafePtr<Framework::PhysicalModelImpl> model) :
  VarSetMatrixTransformer(model),
  _model(model->getConvectiveTerm().d_castTo<EulerTerm>())
{
}

//////////////////////////////////////////////////////////////////////////////

Euler2DPuvtToConsInRef::~Euler2DPuvtToConsInRef()
{
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DPuvtToConsInRef::setMatrixFromRef()
{
  cf_assert(_model.isNotNull());
  const RealVector& linearData = _model->getPhysicalData();

  const CFreal R = _model->getR();
  //unused//  const CFreal p = linearData[EulerTerm::P];
  const CFreal T = linearData[EulerTerm::T];
  const CFreal rho = linearData[EulerTerm::RHO];
  const CFreal rhoOvT = rho/T;
  const CFreal invRT = 1./(R*T);
  const CFreal u = linearData[EulerTerm::VX];
  const CFreal v = linearData[EulerTerm::VY];
  const CFreal V2 = linearData[EulerTerm::V]*linearData[EulerTerm::V];
  const CFreal gammaMinus1 = _model->getGamma() - 1.;

  _transMatrix(0,0) = invRT;
  _transMatrix(0,3) = -rhoOvT;

  _transMatrix(1,0) = u*invRT;
  _transMatrix(1,1) = rho;
  _transMatrix(1,3) = -rhoOvT*u;

  _transMatrix(2,0) = v*invRT;
  _transMatrix(2,2) = rho;
  _transMatrix(2,3) = -rhoOvT*v;

  _transMatrix(3,0) = 1./gammaMinus1 + 0.5*V2*invRT;
  _transMatrix(3,1) = rho*u;
  _transMatrix(3,2) = rho*v;
  _transMatrix(3,3) = -0.5*rhoOvT*V2;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NavierStokes

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
