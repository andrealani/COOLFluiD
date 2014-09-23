#include "KOmega/KOmega.hh"
#include "Euler3DKOmegaPuvtToConsInPuvt.hh"
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

    namespace KOmega {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<Euler3DKOmegaPuvtToConsInPuvt, VarSetMatrixTransformer,
KOmegaModule, 1> euler3DKOmegaPuvtToConsInPuvtProvider("Euler3DKOmegaPuvtToConsInPuvt");

//////////////////////////////////////////////////////////////////////////////

Euler3DKOmegaPuvtToConsInPuvt::Euler3DKOmegaPuvtToConsInPuvt(Common::SafePtr<Framework::PhysicalModelImpl> model) :
  VarSetMatrixTransformer(model),
  _model(model->getConvectiveTerm().d_castTo<EulerKOmegaTerm>())
{
}

//////////////////////////////////////////////////////////////////////////////

Euler3DKOmegaPuvtToConsInPuvt::~Euler3DKOmegaPuvtToConsInPuvt()
{
}

//////////////////////////////////////////////////////////////////////////////

void Euler3DKOmegaPuvtToConsInPuvt::setMatrix(const RealVector& state)
{
  const CFreal R = _model->getR();
  const CFreal p = state[0];
  const CFreal u = state[1];
  const CFreal v = state[2];
  const CFreal w = state[3];
  const CFreal T = state[4];
  const CFreal K = state[5];
  const CFreal Omega = state[6];
  const CFreal rho = p/(R*T);
  const CFreal rhoOvT = rho/T;
  const CFreal invRT = 1./(R*T);
  const CFreal V2 = u*u + v*v + w*w;
  const CFreal gammaMinus1 = _model->getGamma() - 1.;

  _transMatrix(0,0) = invRT;
  _transMatrix(0,4) = -rhoOvT;

  _transMatrix(1,0) = u*invRT;
  _transMatrix(1,1) = rho;
  _transMatrix(1,4) = -rhoOvT*u;

  _transMatrix(2,0) = v*invRT;
  _transMatrix(2,2) = rho;
  _transMatrix(2,4) = -rhoOvT*v;

  _transMatrix(3,0) = w*invRT;
  _transMatrix(3,3) = rho;
  _transMatrix(3,4) = -rhoOvT*w;

  _transMatrix(4,0) = 1./gammaMinus1 + 0.5*V2*invRT;
  _transMatrix(4,1) = rho*u;
  _transMatrix(4,2) = rho*v;
  _transMatrix(4,3) = rho*w;
  _transMatrix(4,4) = -0.5*rhoOvT*V2;
  _transMatrix(4,5) = rho;

  _transMatrix(5,0) = K*invRT;
  _transMatrix(5,4) = -rhoOvT*K;
  _transMatrix(5,5) = rho;

  _transMatrix(6,0) = Omega*invRT;
  _transMatrix(6,4) = -rhoOvT*Omega;
  _transMatrix(6,6) = rho;

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

    } // namespace KOmega

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
