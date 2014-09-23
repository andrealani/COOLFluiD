#include "NavierStokes/NavierStokes.hh"
#include "Euler3DPvtToConsInPvt.hh"
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

Environment::ObjectProvider<Euler3DPvtToConsInPvt, VarSetMatrixTransformer, NavierStokesModule, 1> euler3DPvtToConsInPvtProvider("Euler3DPvtToConsInPvt");

//////////////////////////////////////////////////////////////////////////////

Euler3DPvtToConsInPvt::Euler3DPvtToConsInPvt(Common::SafePtr<Framework::PhysicalModelImpl> model) :
  VarSetMatrixTransformer(model),
  _model(model->getConvectiveTerm().d_castTo<EulerTerm>())
{
}

//////////////////////////////////////////////////////////////////////////////

Euler3DPvtToConsInPvt::~Euler3DPvtToConsInPvt()
{
}

//////////////////////////////////////////////////////////////////////////////

void Euler3DPvtToConsInPvt::setMatrix(const RealVector& state)
{
  const CFreal R = _model->getR();
  const CFreal p = _model->getPressureFromState(state[0]);
  const CFreal T = state[4];
  const CFreal rho = _model->getDensity(p,T);
  const CFreal rhoOvT = rho/T;
  const CFreal invRT = 1./(R*T);
  const CFreal u = state[1];
  const CFreal v = state[2];
  const CFreal w = state[3];
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
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NavierStokes

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
