#include "NavierStokes/NavierStokes.hh"
#include "NavierStokes1DRhovt.hh"
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

Environment::ObjectProvider<NavierStokes1DRhovt, DiffusiveVarSet, NavierStokesModule, 2>
 ns1DRhovtProvider("NavierStokes1DRhovt");

//////////////////////////////////////////////////////////////////////////////

NavierStokes1DRhovt::NavierStokes1DRhovt(const std::string& name, Common::SafePtr<Framework::PhysicalModelImpl> model) :
  NavierStokes1DVarSet(name, model),
  _eulerModel(model->getConvectiveTerm().d_castTo<EulerTerm>())
{
  vector<std::string> names(4);
  names[0] = "rho";
  names[1] = "v";
  names[2] = "T";
  setVarNames(names);
}

//////////////////////////////////////////////////////////////////////////////

NavierStokes1DRhovt::~NavierStokes1DRhovt()
{
}

//////////////////////////////////////////////////////////////////////////////

void NavierStokes1DRhovt::setGradientVars(const vector<RealVector*>& states,
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
    values(0,i) = R*state[0]*values(2,i);
  }
}

//////////////////////////////////////////////////////////////////////////////

void NavierStokes1DRhovt::setGradientVarGradients(const vector<RealVector*>& states,
                                                  const vector< vector<RealVector*> >& stateGradients,
                                                  vector< vector<RealVector*> >& gradVarGradients,
                                                  const CFuint stateSize)
{
  throw Common::NotImplementedException(FromHere(),"setGradientVarGradients");
}

//////////////////////////////////////////////////////////////////////////////

void NavierStokes1DRhovt::setStateGradients(const vector<RealVector*>& states,
                                            const vector< vector<RealVector*> >& gradVarGradients,
                                            vector< vector<RealVector*> >& stateGradients,
                                            const CFuint stateSize)
{
  throw Common::NotImplementedException(FromHere(),"setStateGradients");
}

//////////////////////////////////////////////////////////////////////////////

CFreal NavierStokes1DRhovt::getDynViscosity(const RealVector& state,
                                            const vector<RealVector*>& gradients)
{
  const CFreal Tdim = _eulerModel->getTempRef()*state[2];
  const CFreal pdim = _eulerModel->getPressRef()*state[0]*_eulerModel->getR()*Tdim;
  return getModel().getDynViscosityDim(pdim, Tdim)/
    (getModel().getReferencePhysicalData())[NSTerm::MU];
}

//////////////////////////////////////////////////////////////////////////////

CFreal NavierStokes1DRhovt::getDensity(const RealVector& state)
{
  return state[0];
}

//////////////////////////////////////////////////////////////////////////////

void NavierStokes1DRhovt::setGradientState(const RealVector& state)
{
  cf_assert(_gradState.size() == state.size());
  _gradState[0] = state[0]*_eulerModel->getR()*state[2];
  _gradState[1] = state[1];
  _gradState[2] = state[2];
}

//////////////////////////////////////////////////////////////////////////////

} // namespace NavierStokes

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
