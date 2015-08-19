#include "Poisson/Poisson.hh"
#include "Poisson/PoissonDiff3DCons.hh"
#include "Common/StringOps.hh"
#include "Environment/ObjectProvider.hh"
#include "Poisson/PoissonConvTerm.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace Poisson {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<PoissonDiff3DCons, DiffusiveVarSet, PoissonModule, 2>
poissonDiff3DConsProvider("Poisson3DCons");

//////////////////////////////////////////////////////////////////////////////

PoissonDiff3DCons::PoissonDiff3DCons(const std::string& name, Common::SafePtr<Framework::PhysicalModelImpl> model) :
  PoissonDiffVarSet(name, model),
  m_eulerModel(model->getConvectiveTerm().d_castTo<PTERM>())
{
  
  const CFuint totalNbEqs = 1;
  vector<std::string> names(totalNbEqs);
  
  names[0] = "Phi";
  
  setVarNames(names);
}

//////////////////////////////////////////////////////////////////////////////

PoissonDiff3DCons::~PoissonDiff3DCons()
{
}

//////////////////////////////////////////////////////////////////////////////

void PoissonDiff3DCons::setGradientVars(const vector<RealVector*>& states,
					  RealMatrix& values,
					  const CFuint stateSize)

{
  using namespace std;
  using namespace COOLFluiD::Framework;

  cf_assert(values.nbRows() == PhysicalModelStack::getActive()->getNbEq());
  
  const CFuint nbValues = values.nbRows();
  for (CFuint i = 0; i < nbValues; ++i) {
    for (CFuint j = 0; j < stateSize; ++j) {
      values(i,j) = (*states[j])[i];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void PoissonDiff3DCons::setGradientVarGradients(const vector<RealVector*>& states,
                                                  const vector< vector<RealVector*> >& stateGradients,
                                                  vector< vector<RealVector*> >& gradVarGradients,
                                                  const CFuint stateSize)
{
  using namespace std;
  using namespace COOLFluiD::Framework;
  
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

void PoissonDiff3DCons::setStateGradients(const vector<RealVector*>& states,
                                            const vector< vector<RealVector*> >& gradVarGradients,
                                            vector< vector<RealVector*> >& stateGradients,
                                            const CFuint stateSize)
{
  using namespace std;
  using namespace COOLFluiD::Framework;
  
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

void PoissonDiff3DCons::setGradientState(const RealVector& state)
{
  cf_assert(_gradState.size() == state.size());
  _gradState = state;

}

//////////////////////////////////////////////////////////////////////////////

    } // namespace Poisson

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
