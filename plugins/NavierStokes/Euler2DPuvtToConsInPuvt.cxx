#include "NavierStokes/NavierStokes.hh"
#include "Euler2DPuvtToConsInPuvt.hh"
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

Environment::ObjectProvider<Euler2DPuvtToConsInPuvt, VarSetMatrixTransformer, NavierStokesModule, 1> euler2DPuvtToConsInPuvtProvider("Euler2DPuvtToConsInPuvt");

//////////////////////////////////////////////////////////////////////////////

Euler2DPuvtToConsInPuvt::Euler2DPuvtToConsInPuvt(Common::SafePtr<Framework::PhysicalModelImpl> model) :
  VarSetMatrixTransformer(model),
  _model(model->getConvectiveTerm().d_castTo<EulerTerm>())
{
}

//////////////////////////////////////////////////////////////////////////////

Euler2DPuvtToConsInPuvt::~Euler2DPuvtToConsInPuvt()
{
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DPuvtToConsInPuvt::setMatrix(const RealVector& state)
{
  const CFreal R = _model->getR();
  const CFreal p = _model->getPressureFromState(state[0]);
  const CFreal T = state[3];
  const CFreal rho = _model->getDensity(p,T);
  const CFreal rhoOvT = rho/T;
  const CFreal invRT = 1./(R*T);
  const CFreal u = state[1];
  const CFreal v = state[2];
  const CFreal V2 = u*u + v*v;
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

  // RealMatrix m(4,4,0.);
//   RealMatrix m2(4,4,0.);
//   const CFreal invRho = 1./rho;
//   m(0,0) = (u*u + v*v)*0.5*gammaMinus1;
//   m(0,1) = -gammaMinus1*u;
//   m(0,2) = -gammaMinus1*v;
//   m(0,3) = gammaMinus1;
  
//   m(1,0) = -u*invRho;
//   m(1,1) = invRho;
  
//   m(2,0) = -v*invRho;
//   m(2,2) = invRho;
  
//   m(3,0) = (m(0,0) - p*invRho)*invRho/R;
//   m(3,1) = m(0,1)*invRho/R;
//   m(3,2) = m(0,2)*invRho/R;
//   m(3,3) = m(0,3)*invRho/R;
  
//   m2 = m*_transMatrix;
  
//   CFout << "m2 = "<< "\n";
//   CFout << m2 << "\n" << "\n";
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NavierStokes

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
