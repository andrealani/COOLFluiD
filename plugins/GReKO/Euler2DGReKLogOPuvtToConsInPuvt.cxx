#include "GReKO/GReKO.hh"
#include "Euler2DGReKLogOPuvtToConsInPuvt.hh"
#include "NavierStokes/EulerPhysicalModel.hh"
#include "Framework/PhysicalModel.hh"
#include "Environment/ObjectProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace GReKO {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<Euler2DGReKLogOPuvtToConsInPuvt, VarSetMatrixTransformer,
GReKOModule, 1> Euler2DGReKLogOPuvtToConsInPuvtProvider("Euler2DGReKLogOPuvtToConsInPuvt");

//////////////////////////////////////////////////////////////////////////////

Euler2DGReKLogOPuvtToConsInPuvt::Euler2DGReKLogOPuvtToConsInPuvt(Common::SafePtr<Framework::PhysicalModelImpl> model) :
  VarSetMatrixTransformer(model),
  _model(model->getConvectiveTerm().d_castTo<EulerGReKLogOTerm>())
{
}

//////////////////////////////////////////////////////////////////////////////

Euler2DGReKLogOPuvtToConsInPuvt::~Euler2DGReKLogOPuvtToConsInPuvt()
{
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DGReKLogOPuvtToConsInPuvt::setMatrix(const RealVector& state)
{
  const CFreal R = _model->getR();
  const CFreal p = state[0];
  const CFreal T = state[3];
  const CFreal K = state[4];
  const CFreal Omega = state[5];
  const CFreal Ga = state[6];
  const CFreal Re = state[7];
  const CFreal rho = p/(R*T);
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
  _transMatrix(3,4) = rho;

  _transMatrix(4,0) = K*invRT;
  _transMatrix(4,3) = -rhoOvT*K;
  _transMatrix(4,4) = rho;

  _transMatrix(5,0) = Omega*invRT;
  _transMatrix(5,3) = -rhoOvT*Omega;
  _transMatrix(5,5) = rho;

  _transMatrix(6,0) = Ga*invRT;
  _transMatrix(6,3) = -rhoOvT*Ga;
  _transMatrix(6,6) = rho;
  
  _transMatrix(7,0) = Re*invRT;
  _transMatrix(7,3) = -rhoOvT*Re;
  _transMatrix(7,7) = rho;
  // RealMatrix m(4,4,0.);
  //   RealMatrix m2(4,4,0.);

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

    } // namespace GReKO

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
