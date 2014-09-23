#include "LTE.hh"
#include "NavierStokes3DPvtLTE.hh"
#include "Environment/ObjectProvider.hh"
#include "Framework/PhysicalChemicalLibrary.hh"
#include "NavierStokes/EulerTerm.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace LTE {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<NavierStokes3DPvtLTE, DiffusiveVarSet, LTEModule, 2> ns3DPvtLTEProvider("NavierStokes3DPvtLTE");

//////////////////////////////////////////////////////////////////////////////

NavierStokes3DPvtLTE::NavierStokes3DPvtLTE(const std::string& name, Common::SafePtr<Framework::PhysicalModelImpl> model) :
  NavierStokes3DVarSet(name, model),
  _library(CFNULL),
  _eulerModel(model->getConvectiveTerm().d_castTo<EulerTerm>()),
  _tempX()
{
  vector<std::string> names(5);
  names[0] = (!_eulerModel->isIncompressible()) ? "p" : "dp";
  names[1] = "u";
  names[2] = "v";
  names[3] = "w";
  names[4] = "T";
  setVarNames(names);
}

//////////////////////////////////////////////////////////////////////////////

NavierStokes3DPvtLTE::~NavierStokes3DPvtLTE()
{
}

//////////////////////////////////////////////////////////////////////////////

void NavierStokes3DPvtLTE::setGradientVars(const vector<RealVector*>& states,
					   RealMatrix& values,
					   const CFuint stateSize)
{
  const CFuint nbValues = values.nbRows();
  for (CFuint i = 0; i < nbValues; ++i) {
    for (CFuint j = 0; j < stateSize; ++j) {
      values(i,j) = (*states[j])[i];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void NavierStokes3DPvtLTE::setGradientVarGradients(const vector<RealVector*>& states,
                                                   const vector< vector<RealVector*> >& stateGradients,
                                                   vector< vector<RealVector*> >& gradVarGradients,
                                                   const CFuint stateSize)
{
  cf_assert(stateGradients.size() > 0);
  const CFuint nbValues = stateGradients[0].size();
  for (CFuint i = 0; i < nbValues; ++i)
  {
    for (CFuint j = 0; j < stateSize; ++j)
    {
      *gradVarGradients[j][i] = *stateGradients[j][i];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void NavierStokes3DPvtLTE::setStateGradients(const vector<RealVector*>& states,
                                             const vector< vector<RealVector*> >& gradVarGradients,
                                             vector< vector<RealVector*> >& stateGradients,
                                             const CFuint stateSize)
{
  cf_assert(stateGradients.size() > 0);
  const CFuint nbValues = stateGradients[0].size();
  for (CFuint i = 0; i < nbValues; ++i)
  {
    for (CFuint j = 0; j < stateSize; ++j)
    {
      *stateGradients[j][i] = *gradVarGradients[j][i];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

CFreal NavierStokes3DPvtLTE::getDynViscosity(const RealVector& state, const vector<RealVector*>& gradients)
{
  CFreal Tdim = _eulerModel->getTempRef()*state[4];
  CFreal pdim = _eulerModel->getPressureFromState(state[0])*
    (_eulerModel->getReferencePhysicalData())[EulerTerm::P];
 return _library->eta(Tdim, pdim, CFNULL) /
   (getModel().getReferencePhysicalData())[NSTerm::MU];
}

//////////////////////////////////////////////////////////////////////////////

CFreal NavierStokes3DPvtLTE::getDensity(const RealVector& state)
{
  CFreal Tdim = _eulerModel->getTempRef()*state[4];
  CFreal pdim = _eulerModel->getPressureFromState(state[0])*
    (_eulerModel->getReferencePhysicalData())[EulerTerm::P];
  const CFreal rhoRef = (_eulerModel->getReferencePhysicalData())[EulerTerm::RHO];
  return _library->density(Tdim,pdim)/rhoRef;
}

//////////////////////////////////////////////////////////////////////////////

CFreal NavierStokes3DPvtLTE::getThermConductivity(const RealVector& state,
						  const CFreal& dynViscosity)
{
  CFreal Tdim = _eulerModel->getTempRef()*state[4];
  CFreal pdim = _eulerModel->getPressureFromState(state[0])*
    (_eulerModel->getReferencePhysicalData())[EulerTerm::P];
  return _library->lambdaEQ(Tdim,pdim) /
    (getModel().getReferencePhysicalData())[NSTerm::LAMBDA];
}

//////////////////////////////////////////////////////////////////////////////

void NavierStokes3DPvtLTE::setup()
{
  //  fluct split has to call setup()
  NavierStokes3DVarSet::setup();

  _library = PhysicalModelStack::getActive()->getImplementor()->
    getPhysicalPropertyLibrary<PhysicalChemicalLibrary>();
  cf_assert(_library.isNotNull());

  _tempX.resize(_library->getNbSpecies());
}

//////////////////////////////////////////////////////////////////////////////

void NavierStokes3DPvtLTE::setComposition(const RealVector& state,
					  const bool isPerturb,
					  const CFuint iVar)
{
  // this is to avoid useless expensive re-computations
  useBackUpValues(false);
  setBackUpValues(false);

  if (isPerturb && (iVar == 1 || iVar == 2 || iVar == 3)) {
    useBackUpValues(true);
  }
  else if (isPerturb && (iVar == 0 || iVar == 4)) {
    _library->resetComposition(_tempX);
  }
  else if (!isPerturb) {
    CFreal Tdim = _eulerModel->getTempRef()*state[4];
    CFreal pdim = _eulerModel->getPressureFromState(state[0])*
    (_eulerModel->getReferencePhysicalData())[EulerTerm::P];
    _library->setComposition(Tdim,pdim,&_tempX);
    // set and store the back up values only if an unperturbed flux
    // is computed
    setBackUpValues(true);
  }
}

//////////////////////////////////////////////////////////////////////////////

void NavierStokes3DPvtLTE::setGradientState(const RealVector& state)
{
  cf_assert(_gradState.size() == state.size());
  _gradState = state;
}

//////////////////////////////////////////////////////////////////////////////

} // namespace LTE

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
