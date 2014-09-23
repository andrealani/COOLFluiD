#include "NavierStokes/NavierStokes.hh"
#include "Euler2DRhovtToConsInRhovt.hh"
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

Environment::ObjectProvider<Euler2DRhovtToConsInRhovt, VarSetMatrixTransformer, 
			    NavierStokesModule, 1> 
euler2DRhovtToConsInRhovtProvider("Euler2DRhovtToConsInRhovt");

//////////////////////////////////////////////////////////////////////////////

Euler2DRhovtToConsInRhovt::Euler2DRhovtToConsInRhovt(Common::SafePtr<Framework::PhysicalModelImpl> model) :
  VarSetMatrixTransformer(model),
  _model(model->getConvectiveTerm().d_castTo<EulerTerm>())
{
}

//////////////////////////////////////////////////////////////////////////////

Euler2DRhovtToConsInRhovt::~Euler2DRhovtToConsInRhovt()
{
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DRhovtToConsInRhovt::setMatrix(const RealVector& state)
{
  const CFreal R = _model->getR();
  const CFreal rho = state[0];
  const CFreal T = state[3];
  const CFreal u = state[1];
  const CFreal v = state[2];
  const CFreal V2 = u*u + v*v;
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
