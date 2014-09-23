#include "NavierStokes/NavierStokes.hh"
#include "Euler2DRhovtToConsInRef.hh"
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

Environment::ObjectProvider<Euler2DRhovtToConsInRef, VarSetMatrixTransformer, 
NavierStokesModule, 1> euler2DRhovtToConsInRefProvider("Euler2DRhovtToConsInRef");

//////////////////////////////////////////////////////////////////////////////

Euler2DRhovtToConsInRef::Euler2DRhovtToConsInRef(Common::SafePtr<Framework::PhysicalModelImpl> model) :
  VarSetMatrixTransformer(model),
  _model(model->getConvectiveTerm().d_castTo<EulerTerm>())
{
}

//////////////////////////////////////////////////////////////////////////////

Euler2DRhovtToConsInRef::~Euler2DRhovtToConsInRef()
{
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DRhovtToConsInRef::setMatrixFromRef()
{
  cf_assert(_model.isNotNull());
  const RealVector& linearData = _model->getPhysicalData();
  
  const CFreal R = _model->getR();
  const CFreal rho = linearData[EulerTerm::RHO];
  const CFreal T = linearData[EulerTerm::T];
  cf_assert(T > 0.);
  
  const CFreal u = linearData[EulerTerm::VX];
  const CFreal v = linearData[EulerTerm::VY];
  const CFreal V2 = linearData[EulerTerm::V]*linearData[EulerTerm::V];
  const CFreal ovGammaMinus1 = 1./(_model->getGamma() - 1.);
  
  _transMatrix(0,0) = 1.0;
  
  _transMatrix(1,0) = u;
  _transMatrix(1,1) = rho;
  
  _transMatrix(2,0) = v;
  _transMatrix(2,2) = rho;
  
  _transMatrix(3,0) = R*T*ovGammaMinus1 + 0.5*V2;
  _transMatrix(3,1) = rho*u;
  _transMatrix(3,2) = rho*v;
  _transMatrix(3,3) = rho*R*ovGammaMinus1;  
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NavierStokes

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
