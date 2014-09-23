#include "LES/LESModule.hh"
#include "LES2DPuvt.hh"
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

Environment::ObjectProvider<LES2DPuvt, DiffusiveVarSet,
	                    LESModule, 2>
les2DPuvtProvider("LES2DPuvt");

//////////////////////////////////////////////////////////////////////////////

LES2DPuvt::LES2DPuvt
(const std::string& name, Common::SafePtr<Framework::PhysicalModelImpl> model) :
  LES2DVarSet(name, model),
  m_eulerModel(model->getConvectiveTerm().d_castTo<Physics::NavierStokes::EulerTerm>())
{
  vector<std::string> names(4);
  names[0] = "p";
  names[1] = "u";
  names[2] = "v";
  names[3] = "T";
  setVarNames(names);
}

//////////////////////////////////////////////////////////////////////////////

LES2DPuvt::~LES2DPuvt()
{
}

//////////////////////////////////////////////////////////////////////////////

void LES2DPuvt::setGradientVars(const vector<RealVector*>& states,
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

void LES2DPuvt::setGradientVarGradients(const vector<RealVector*>& states,
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

void LES2DPuvt::setStateGradients(const vector<RealVector*>& states,
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

CFreal LES2DPuvt::getDynViscosity(const RealVector& state,
                                  const vector<RealVector*>& gradients)
{
  const CFreal Tdim = m_eulerModel->getTempRef()*state[3];
  const CFreal pdim = m_eulerModel->getPressRef()*state[0];

  return getModel().getDynViscosityDim(pdim, Tdim)/
      (getModel().getReferencePhysicalData())[Physics::NavierStokes::NSTerm::MU];
}

//////////////////////////////////////////////////////////////////////////////

CFreal LES2DPuvt::getDensity(const RealVector& state)
{
  return state[0]/(m_eulerModel->getR()*state[3]);
}

//////////////////////////////////////////////////////////////////////////////

void LES2DPuvt::setGradientState(const RealVector& state)
{
  cf_assert(_gradState.size() == state.size());
  _gradState = state;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace LES

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
