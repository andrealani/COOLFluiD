#include "NavierStokes/NavierStokes.hh"
#include "NavierStokes2DRhovt.hh"
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

Environment::ObjectProvider<NavierStokes2DRhovt, DiffusiveVarSet, NavierStokesModule, 2>
 ns2DRhovtProvider("NavierStokes2DRhovt");

//////////////////////////////////////////////////////////////////////////////

NavierStokes2DRhovt::NavierStokes2DRhovt(const std::string& name, Common::SafePtr<Framework::PhysicalModelImpl> model) :
  NavierStokes2DVarSet(name, model),
  _eulerModel(model->getConvectiveTerm().d_castTo<EulerTerm>())
{
  vector<std::string> names(4);
  names[0] = "rho";
  names[1] = "u";
  names[2] = "v";
  names[3] = "T";
  setVarNames(names);
}

//////////////////////////////////////////////////////////////////////////////

NavierStokes2DRhovt::~NavierStokes2DRhovt()
{
}

//////////////////////////////////////////////////////////////////////////////

void NavierStokes2DRhovt::setGradientVars(const vector<RealVector*>& states,
					  RealMatrix& values,
					  const CFuint stateSize)
{
  const CFreal R = _eulerModel->getR();
  const CFuint nbStates = stateSize;

  cf_assert(values.nbRows() == PhysicalModelStack::getActive()->getNbEq());

  for (CFuint i = 0; i < nbStates; ++i) {
    const RealVector& state = *states[i];
    values(1,i) = state[1];
    values(2,i) = state[2];
    values(3,i) = state[3];
    values(0,i) = R*state[0]*values(3,i);
  }
}

//////////////////////////////////////////////////////////////////////////////

void NavierStokes2DRhovt::setGradientVarGradients(const vector<RealVector*>& states,
                                                  const vector< vector<RealVector*> >& stateGradients,
                                                  vector< vector<RealVector*> >& gradVarGradients,
                                                  const CFuint stateSize)
{
  throw Common::NotImplementedException(FromHere(),"setGradientVarGradients");
}

//////////////////////////////////////////////////////////////////////////////

void NavierStokes2DRhovt::setStateGradients(const vector<RealVector*>& states,
                                            const vector< vector<RealVector*> >& gradVarGradients,
                                            vector< vector<RealVector*> >& stateGradients,
                                            const CFuint stateSize)
{
  throw Common::NotImplementedException(FromHere(),"setStateGradients");
}

//////////////////////////////////////////////////////////////////////////////

CFreal NavierStokes2DRhovt::getDynViscosity(const RealVector& state,
                                            const vector<RealVector*>& gradients)
{
  const CFreal Tdim = _eulerModel->getTempRef()*state[3];
  const CFreal pdim = _eulerModel->getPressRef()*state[0]*_eulerModel->getR()*Tdim;
  return getModel().getDynViscosityDim(pdim, Tdim)/
    (getModel().getReferencePhysicalData())[NSTerm::MU];
}

//////////////////////////////////////////////////////////////////////////////

CFreal NavierStokes2DRhovt::getDensity(const RealVector& state)
{
  return state[0];
}

//////////////////////////////////////////////////////////////////////////////

void NavierStokes2DRhovt::setGradientState(const RealVector& state)
{
  cf_assert(_gradState.size() == state.size());
  _gradState[0] = state[0]*_eulerModel->getR()*state[3];
  _gradState[1] = state[1];
  _gradState[2] = state[2];
  _gradState[3] = state[3];
}

//////////////////////////////////////////////////////////////////////////////

} // namespace NavierStokes

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
