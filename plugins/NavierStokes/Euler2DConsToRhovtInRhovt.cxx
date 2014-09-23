#include "NavierStokes/NavierStokes.hh"
#include "Euler2DConsToRhovtInRhovt.hh"
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

Environment::ObjectProvider<Euler2DConsToRhovtInRhovt, VarSetMatrixTransformer, 
			    NavierStokesModule, 1> 
euler2DConsToRhovtInRhovtProvider("Euler2DConsToRhovtInRhovt");

//////////////////////////////////////////////////////////////////////////////

Euler2DConsToRhovtInRhovt::Euler2DConsToRhovtInRhovt(Common::SafePtr<Framework::PhysicalModelImpl> model) :
  VarSetMatrixTransformer(model),
  _model(model->getConvectiveTerm().d_castTo<EulerTerm>())
{
}

//////////////////////////////////////////////////////////////////////////////

Euler2DConsToRhovtInRhovt::~Euler2DConsToRhovtInRhovt()
{
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DConsToRhovtInRhovt::setMatrix(const RealVector& state)
{
  const CFreal R = _model->getR();
  const CFreal rho = state[0];
  const CFreal invRho = 1/rho;
  const CFreal u = state[1];
  const CFreal v = state[2];
  const CFreal T = state[3];
  const CFreal V2 = u*u + v*v;
  const CFreal gammaMinus1 = _model->getGamma() - 1.;
  const CFreal kg = gammaMinus1/(R*rho);
  
  _transMatrix(0,0) = 1.0;

  _transMatrix(1,0) = -u*invRho;
  _transMatrix(1,1) = invRho;
  
  _transMatrix(2,0) = -v*invRho;
  _transMatrix(2,2) = invRho;
  
  _transMatrix(3,0) = kg*(0.5*V2 - R*T/gammaMinus1);
  _transMatrix(3,1) = -kg*u;
  _transMatrix(3,2) = -kg*v;
  _transMatrix(3,3) = kg; 
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NavierStokes

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
