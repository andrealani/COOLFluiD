#include "LESvki/LESvki.hh"
#include "Clark2DCons.hh"
#include "Environment/ObjectProvider.hh"
#include "NavierStokes/EulerTerm.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace LESvki {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<Clark2DCons, DiffusiveVarSet,
	       LESvkiModule, 2>
Clark2DConsProvider("Clark2DCons");

//////////////////////////////////////////////////////////////////////////////

Clark2DCons::Clark2DCons
(const std::string& name, Common::SafePtr<Framework::PhysicalModelImpl> model) :
  Clark2DVarSet(name, model),
  _eulerModel(model->getConvectiveTerm().d_castTo<NavierStokes::EulerTerm>())
{
  vector<std::string> names(4);
  names[0] = "rho";
  names[1] = "rhoU";
  names[2] = "rhoV";
  names[3] = "rhoE";
  setVarNames(names);
}

//////////////////////////////////////////////////////////////////////////////

Clark2DCons::~Clark2DCons()
{
}

//////////////////////////////////////////////////////////////////////////////

void Clark2DCons::setGradientVars(const vector<RealVector*>& states,
					 RealMatrix& values,
					 const CFuint stateSize)
{
  cf_assert(values.nbRows() == PhysicalModelStack::getActive()->getNbEq());
  
  const CFreal R = _eulerModel->getR();
  const CFreal ovCv = (_eulerModel->getGamma() - 1.)/R;
  
  for (CFuint i = 0; i < stateSize; ++i) {
    const RealVector& state = *states[i];
    const CFreal ovRho = 1./state[0]; 
    values(1,i) = state[1]*ovRho;
    values(2,i) = state[2]*ovRho;
    
    const CFreal V2 = values(1,i)*values(1,i) + values(2,i)*values(2,i);
    values(3,i) = (state[3]*ovRho - 0.5*V2)*ovCv;
    values(0,i) = R*state[0]*values(3,i);
  }
}

//////////////////////////////////////////////////////////////////////////////

void Clark2DCons::setGradientVarGradients(const vector<RealVector*>& states,
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
    const CFreal rhoV     = (*states[j])[2];
    const CFreal rhoE     = (*states[j])[3];

    // dereference state gradients
    const RealVector& gradRho  = *stateGradients[j][0];
    const RealVector& gradRhoU = *stateGradients[j][1];
    const RealVector& gradRhoV = *stateGradients[j][2];
    const RealVector& gradRhoE = *stateGradients[j][3];

    // helper variables
    const CFreal rhoInv   = 1.0/rho;
    const CFreal rhoInvSq = rhoInv*rhoInv;
    const CFreal rhoVelMagSq = rhoU*rhoU + rhoV*rhoV;
    const RealVector rhoVelTimesGradRhoVel = rhoU*gradRhoU + rhoV*gradRhoV;

    // pressure gradient
    *gradVarGradients[j][0] = gammaM1*(
                                       0.5*rhoVelMagSq*rhoInvSq*gradRho
                                      - rhoInv*rhoVelTimesGradRhoVel
                                      + gradRhoE
                                      );

    // velocity gradients
    *gradVarGradients[j][1] = rhoInv*gradRhoU - rhoU*rhoInvSq*gradRho;
    *gradVarGradients[j][2] = rhoInv*gradRhoV - rhoV*rhoInvSq*gradRho;

    // temperature gradient
    *gradVarGradients[j][3] = cvInv*rhoInvSq*(
                                              (rhoVelMagSq*rhoInv - rhoE)*gradRho
                                              - rhoVelTimesGradRhoVel
                                              + rho*gradRhoE
                                             );
  }
}

//////////////////////////////////////////////////////////////////////////////

void Clark2DCons::setStateGradients(const vector<RealVector*>& states,
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
    const CFreal rhoV     = (*states[j])[2];
    const CFreal rhoE     = (*states[j])[3];

    // dereference gradient variable gradients
    const RealVector& gradP = *gradVarGradients[j][0];
    const RealVector& gradU = *gradVarGradients[j][1];
    const RealVector& gradV = *gradVarGradients[j][2];
    const RealVector& gradT = *gradVarGradients[j][3];

    // helper variables
    const CFreal rhoInv = 1.0/rho;
    const CFreal u = rhoU*rhoInv;
    const CFreal v = rhoV*rhoInv;
    const CFreal eKin = 0.5*(u*u + v*v);
    const CFreal p = gammaM1*(rhoE - rho*eKin);
    const CFreal T = p*rhoInv*RInv;

    // mass density gradient
    RealVector& gradRho = *stateGradients[j][0];
    gradRho = rho*(gradP/p - gradT/T);

    // momentum gradients
    *stateGradients[j][1] = rho*gradU + u*gradRho;
    *stateGradients[j][2] = rho*gradV + v*gradRho;

    // total energy gradient
    *stateGradients[j][3] = gammaM1Inv*gradP + rho*(u*gradU + v*gradV) + eKin*gradRho;
  }
}

//////////////////////////////////////////////////////////////////////////////

CFreal Clark2DCons::getDynViscosity(const RealVector& state,
                                           const vector<RealVector*>& gradients)
{
  // here it is assumed that state is in Cons variables
  const CFreal R = _eulerModel->getR();
  const CFreal cv = R/(_eulerModel->getGamma() - 1.);
  const CFreal u = state[1]/state[0];
  const CFreal v = state[2]/state[0];
  const CFreal T = (state[3] - 0.5*state[0]*(u*u + v*v))/(state[0]*cv);
  const CFreal p = R*state[0]*T;
  const CFreal Tdim = _eulerModel->getTempRef()*T;
  const CFreal pdim = _eulerModel->getPressRef()*p;

  return getModel().getDynViscosityDim(pdim, Tdim)/
    (getModel().getReferencePhysicalData())[NSTerm::MU];
}

//////////////////////////////////////////////////////////////////////////////

CFreal Clark2DCons::getDensity(const RealVector& state)
{
  return state[0];
}

//////////////////////////////////////////////////////////////////////////////

void Clark2DCons::setGradientState(const RealVector& state)
{
  cf_assert(_gradState.size() == state.size());
  const CFreal R = _eulerModel->getR();
  const CFreal cv = R/(_eulerModel->getGamma() - 1.);

  // _gradState = [p u v T]
  _gradState[1] = state[1]/state[0];
  _gradState[2] = state[2]/state[0];
  const CFreal V2 = _gradState[1]*_gradState[1] + _gradState[2]*_gradState[2];

  _gradState[3] = (state[3] - 0.5*state[0]*V2)/(state[0]*cv);
  _gradState[0] = R*state[0]*_gradState[3];
}


////////////////////////////////////////////////////////////////////////////////
CFreal Clark2DCons::getTurbDynViscosityFromGradientVars(const RealVector& state, const vector<RealVector*>& gradients)
{
 const RealVector& gradU = *gradients[1];
  const RealVector& gradV = *gradients[2];
  const CFreal Cs = 0.1; // Smagorinsky constant
  const CFreal SXX = gradU[XX];
  const CFreal SXY = 0.5*(gradU[YY]+gradV[XX]);
  const CFreal SYY = gradV[YY];

  const CFreal S = sqrt(2.*pow(SXX,2)+2.*pow(SYY,2)+pow(SXY,2)); 
  const CFreal mut = pow(Cs,2)*S; // Smagorinsky Subgrid viscosity
  return mut;
}
 
/////////////////////////////////////////////////////////////////////////////

    } // namespace NavierStokes

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
