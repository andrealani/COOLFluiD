#include "LES/LESModule.hh"
#include "LES3DCons.hh"
#include "Environment/ObjectProvider.hh"
#include "NavierStokes/EulerTerm.hh"
//#include "NavierStokes/NSTerm.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace LES {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<LES3DCons, DiffusiveVarSet,
	                    LESModule, 2>
les3DConsProvider("LES3DCons");

//////////////////////////////////////////////////////////////////////////////

LES3DCons::LES3DCons
(const std::string& name, Common::SafePtr<Framework::PhysicalModelImpl> model) :
  LES3DVarSet(name, model)
{
  vector<std::string> names(5);
  names[0] = "rho";
  names[1] = "rhoU";
  names[2] = "rhoV";
  names[3] = "rhoW";
  names[4] = "rhoE";
  setVarNames(names);
}

//////////////////////////////////////////////////////////////////////////////

LES3DCons::~LES3DCons()
{
}

//////////////////////////////////////////////////////////////////////////////

void LES3DCons::setGradientVars(const vector<RealVector*>& states,
                                RealMatrix& values,
                                const CFuint stateSize)
{
  const CFreal R = m_eulerModel->getR();
  const CFreal cv = R/(m_eulerModel->getGamma() - 1.);
  const CFuint nbStates = stateSize;

  cf_assert(values.nbRows() == PhysicalModelStack::getActive()->getNbEq());

  for (CFuint i = 0; i < nbStates; ++i)
  {
    const RealVector& state = *states[i];
    values(1,i) = state[1]/state[0];
    values(2,i) = state[2]/state[0];
    values(3,i) = state[3]/state[0];
    values(4,i) = (state[4] - 0.5*state[0]*(
                                            values(1,i)*values(1,i) +
                                            values(2,i)*values(2,i) +
                                            values(3,i)*values(3,i)
                                           )
                  )/(state[0]*cv);
    values(0,i) = R*state[0]*values(4,i);
  }
}

//////////////////////////////////////////////////////////////////////////////

void LES3DCons::setGradientVarGradients(const vector<RealVector*>& states,
                                        const vector< vector<RealVector*> >& stateGradients,
                                        vector< vector<RealVector*> >& gradVarGradients,
                                        const CFuint stateSize)
{
  // get heat capacity at constant volume
  const CFreal R = m_eulerModel->getR();
  const CFreal gammaM1 = m_eulerModel->getGamma() - 1.;
  const CFreal cvInv = gammaM1/R;

  for (CFuint j = 0; j < stateSize; ++j)
  {
    // dereference states
    const CFreal rho      = (*states[j])[0];
    const CFreal rhoU     = (*states[j])[1];
    const CFreal rhoV     = (*states[j])[2];
    const CFreal rhoW     = (*states[j])[3];
    const CFreal rhoE     = (*states[j])[4];

    // dereference state gradients
    const RealVector& gradRho  = *stateGradients[j][0];
    const RealVector& gradRhoU = *stateGradients[j][1];
    const RealVector& gradRhoV = *stateGradients[j][2];
    const RealVector& gradRhoW = *stateGradients[j][3];
    const RealVector& gradRhoE = *stateGradients[j][4];

    // helper variables
    const CFreal rhoInv   = 1.0/rho;
    const CFreal rhoInvSq = rhoInv*rhoInv;
    const CFreal rhoVelMagSq = rhoU*rhoU + rhoV*rhoV + rhoW*rhoW;
    const RealVector rhoVelTimesGradRhoVel = rhoU*gradRhoU + rhoV*gradRhoV + rhoW*gradRhoW;

    // pressure gradient
    *gradVarGradients[j][0] = gammaM1*(
                                       0.5*rhoVelMagSq*rhoInvSq*gradRho
                                       - rhoInv*rhoVelTimesGradRhoVel
                                       + gradRhoE
                                      );

    // velocity gradients
    *gradVarGradients[j][1] = rhoInv*gradRhoU - rhoU*rhoInvSq*gradRho;
    *gradVarGradients[j][2] = rhoInv*gradRhoV - rhoV*rhoInvSq*gradRho;
    *gradVarGradients[j][3] = rhoInv*gradRhoW - rhoW*rhoInvSq*gradRho;

    // temperature gradient
    *gradVarGradients[j][4] = cvInv*rhoInvSq*(
                                              (rhoVelMagSq*rhoInv - rhoE)*gradRho
                                              - rhoVelTimesGradRhoVel
                                              + rho*gradRhoE
                                             );
  }
}

//////////////////////////////////////////////////////////////////////////////

void LES3DCons::setStateGradients(const vector<RealVector*>& states,
                                           const vector< std::vector<RealVector*> >& gradVarGradients,
                                           vector< vector<RealVector*> >& stateGradients,
                                           const CFuint stateSize)
{
  // get heat capacity at constant volume
  const CFreal R = m_eulerModel->getR();
  const CFreal RInv = 1.0/R;
  const CFreal gammaM1 = m_eulerModel->getGamma() - 1.;
  const CFreal gammaM1Inv = 1.0/gammaM1;

  for (CFuint j = 0; j < stateSize; ++j)
  {
    // dereference states
    const CFreal rho      = (*states[j])[0];
    const CFreal rhoU     = (*states[j])[1];
    const CFreal rhoV     = (*states[j])[2];
    const CFreal rhoW     = (*states[j])[3];
    const CFreal rhoE     = (*states[j])[4];

    // dereference gradient variable gradients
    const RealVector& gradP = *gradVarGradients[j][0];
    const RealVector& gradU = *gradVarGradients[j][1];
    const RealVector& gradV = *gradVarGradients[j][2];
    const RealVector& gradW = *gradVarGradients[j][3];
    const RealVector& gradT = *gradVarGradients[j][4];

    // helper variables
    const CFreal rhoInv = 1.0/rho;
    const CFreal u = rhoU*rhoInv;
    const CFreal v = rhoV*rhoInv;
    const CFreal w = rhoW*rhoInv;
    const CFreal eKin = 0.5*(u*u + v*v + w*w);
    const CFreal p = gammaM1*(rhoE - rho*eKin);
    const CFreal T = p*rhoInv*RInv;

    // mass density gradient
    RealVector& gradRho = *stateGradients[j][0];
    gradRho = rho*(gradP/p - gradT/T);

    // momentum gradients
    *stateGradients[j][1] = rho*gradU + u*gradRho;
    *stateGradients[j][2] = rho*gradV + v*gradRho;
    *stateGradients[j][3] = rho*gradW + w*gradRho;

    // total energy gradient
    *stateGradients[j][4] = gammaM1Inv*gradP + rho*(u*gradU + v*gradV + w*gradW) + eKin*gradRho;
  }
}

//////////////////////////////////////////////////////////////////////////////

CFreal LES3DCons::getDynViscosity(const RealVector& state,
                                  const vector<RealVector*>& gradients)
{
  // here it is assumed that state is in Cons variables
  const CFreal R = m_eulerModel->getR();
  const CFreal cv = R/(m_eulerModel->getGamma() - 1.);
  const CFreal u = state[1]/state[0];
  const CFreal v = state[2]/state[0];
  const CFreal w = state[3]/state[0];
  const CFreal T = (state[4] - 0.5*state[0]*(u*u + v*v + w*w))/(state[0]*cv);
  const CFreal p = R*state[0]*T;
  const CFreal Tdim = m_eulerModel->getTempRef()*T;
  const CFreal pdim = m_eulerModel->getPressRef()*p;

  return getModel().getDynViscosityDim(pdim, Tdim)/
    (getModel().getReferencePhysicalData())[Physics::NavierStokes::NSTerm::MU];
}

//////////////////////////////////////////////////////////////////////////////

CFreal LES3DCons::getDensity(const RealVector& state)
{
  return state[0];
}

//////////////////////////////////////////////////////////////////////////////

void LES3DCons::setGradientState(const RealVector& state)
{
  cf_assert(_gradState.size() == state.size());
  const CFreal R = m_eulerModel->getR();
  const CFreal cv = R/(m_eulerModel->getGamma() - 1.);

  // _gradState = [p u v w T]
  _gradState[1] = state[1]/state[0];
  _gradState[2] = state[2]/state[0];
  _gradState[3] = state[3]/state[0];
  const CFreal V2 = _gradState[1]*_gradState[1] + _gradState[2]*_gradState[2] + _gradState[3]*_gradState[3];

  _gradState[4] = (state[4] - 0.5*state[0]*V2)/(state[0]*cv);
  _gradState[0] = R*state[0]*_gradState[4];
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace LES

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
