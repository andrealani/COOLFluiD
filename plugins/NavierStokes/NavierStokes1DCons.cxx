#include "NavierStokes/NavierStokes.hh"
#include "NavierStokes1DCons.hh"
#include "Environment/ObjectProvider.hh"
#include "NavierStokes/EulerTerm.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace NavierStokes {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<NavierStokes1DCons, DiffusiveVarSet,
	       NavierStokesModule, 2>
ns1DConsProvider("NavierStokes1DCons");

//////////////////////////////////////////////////////////////////////////////

NavierStokes1DCons::NavierStokes1DCons
(const std::string& name, Common::SafePtr<Framework::PhysicalModelImpl> model) :
  NavierStokes1DVarSet(name, model),
  _eulerModel(model->getConvectiveTerm().d_castTo<EulerTerm>())
{
  vector<std::string> names(3);
  names[0] = "rho";
  names[1] = "rhoU";
  names[2] = "rhoE";
  setVarNames(names);
}

//////////////////////////////////////////////////////////////////////////////

NavierStokes1DCons::~NavierStokes1DCons()
{
}

//////////////////////////////////////////////////////////////////////////////

void NavierStokes1DCons::setGradientVars(const vector<RealVector*>& states,
					 RealMatrix& values,
					 const CFuint stateSize)
{
  const CFreal R = _eulerModel->getR();
  const CFreal cv = R/(_eulerModel->getGamma() - 1.);
  const CFuint nbStates = stateSize;

  cf_assert(values.nbRows() == PhysicalModelStack::getActive()->getNbEq());

  for (CFuint i = 0; i < nbStates; ++i) {
    const RealVector& state = *states[i];
    values(1,i) = state[1]/state[0];
    values(2,i) = (state[2] - 0.5*state[0]*values(1,i)*values(1,i))/(state[0]*cv);
    values(0,i) = R*state[0]*values(2,i);
  }
}

//////////////////////////////////////////////////////////////////////////////

void NavierStokes1DCons::setGradientVarGradients(const vector<RealVector*>& states,
                                                 const vector< vector<RealVector*> >& stateGradients,
                                                 vector< vector<RealVector*> >& gradVarGradients,
                                                 const CFuint stateSize)
{
  // get heat capacity at constant volume
  const CFreal R = _eulerModel->getR();
  const CFreal gammaM1 = _eulerModel->getGamma() - 1.;
  const CFreal cvInv = gammaM1/R;

  for (CFuint j = 0; j < stateSize; ++j)
  {
    // dereference states
    const CFreal rho      = (*states[j])[0];
    const CFreal rhoU     = (*states[j])[1];
    const CFreal rhoE     = (*states[j])[2];

    // dereference state gradients
    const RealVector& gradRho  = *stateGradients[j][0];
    const RealVector& gradRhoU = *stateGradients[j][1];
    const RealVector& gradRhoE = *stateGradients[j][2];

    // helper variables
    const CFreal rhoInv   = 1.0/rho;
    const CFreal rhoInvSq = rhoInv*rhoInv;
    const CFreal rhoVelMagSq = rhoU*rhoU;
    const RealVector rhoVelTimesGradRhoVel = rhoU*gradRhoU;

    // pressure gradient
    *gradVarGradients[j][0] = gammaM1*(
                                       0.5*rhoVelMagSq*rhoInvSq*gradRho
                                      - rhoInv*rhoVelTimesGradRhoVel
                                      + gradRhoE
                                      );

    // velocity gradient
    *gradVarGradients[j][1] = rhoInv*gradRhoU - rhoU*rhoInvSq*gradRho;

    // temperature gradient
    *gradVarGradients[j][2] = cvInv*rhoInvSq*(
                                              (rhoVelMagSq*rhoInv - rhoE)*gradRho
                                              - rhoVelTimesGradRhoVel
                                              + rho*gradRhoE
                                             );
  }
}

//////////////////////////////////////////////////////////////////////////////

void NavierStokes1DCons::setStateGradients(const vector<RealVector*>& states,
                                           const vector< std::vector<RealVector*> >& gradVarGradients,
                                           vector< vector<RealVector*> >& stateGradients,
                                           const CFuint stateSize)
{
  // get heat capacity at constant volume
  const CFreal R = _eulerModel->getR();
  const CFreal RInv = 1.0/R;
  const CFreal gammaM1 = _eulerModel->getGamma() - 1.;
  const CFreal gammaM1Inv = 1.0/gammaM1;

  for (CFuint j = 0; j < stateSize; ++j)
  {
    // dereference states
    const CFreal rho      = (*states[j])[0];
    const CFreal rhoU     = (*states[j])[1];
    const CFreal rhoE     = (*states[j])[2];

    // dereference gradient variable gradients
    const RealVector& gradP = *gradVarGradients[j][0];
    const RealVector& gradU = *gradVarGradients[j][1];
    const RealVector& gradT = *gradVarGradients[j][2];

    // helper variables
    const CFreal rhoInv = 1.0/rho;
    const CFreal u = rhoU*rhoInv;
    const CFreal eKin = 0.5*u*u;
    const CFreal p = gammaM1*(rhoE - rho*eKin);
    const CFreal T = p*rhoInv*RInv;

    // mass density gradient
    RealVector& gradRho = *stateGradients[j][0];
    gradRho = rho*(gradP/p - gradT/T);

    // momentum gradient
    *stateGradients[j][1] = rho*gradU + u*gradRho;

    // total energy gradient
    *stateGradients[j][2] = gammaM1Inv*gradP + rho*u*gradU + eKin*gradRho;
  }
}

//////////////////////////////////////////////////////////////////////////////

CFreal NavierStokes1DCons::getDynViscosity(const RealVector& state,
                                           const vector<RealVector*>& gradients)
{
  // here it is assumed that state is in Cons variables
  const CFreal R = _eulerModel->getR();
  const CFreal cv = R/(_eulerModel->getGamma() - 1.);
  const CFreal u = state[1]/state[0];
  const CFreal T = (state[2] - 0.5*state[0]*u*u)/(state[0]*cv);
  const CFreal p = R*state[0]*T;
  const CFreal Tdim = _eulerModel->getTempRef()*T;
  const CFreal pdim = _eulerModel->getPressRef()*p;

  return getModel().getDynViscosityDim(pdim, Tdim)/
    (getModel().getReferencePhysicalData())[NSTerm::MU];
}

//////////////////////////////////////////////////////////////////////////////

CFreal NavierStokes1DCons::getDensity(const RealVector& state)
{
  return state[0];
}

//////////////////////////////////////////////////////////////////////////////

void NavierStokes1DCons::setGradientState(const RealVector& state)
{
  cf_assert(_gradState.size() == state.size());
  const CFreal R = _eulerModel->getR();
  const CFreal cv = R/(_eulerModel->getGamma() - 1.);

  // _gradState = [p u v T]
  _gradState[1] = state[1]/state[0];
  const CFreal V2 = _gradState[1]*_gradState[1];

  _gradState[2] = (state[2] - 0.5*state[0]*V2)/(state[0]*cv);
  _gradState[0] = R*state[0]*_gradState[2];
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NavierStokes

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
