#include "NavierStokes/NavierStokes.hh"
#include "NavierStokes1DPvt.hh"
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

Environment::ObjectProvider<NavierStokes1DPvt, DiffusiveVarSet, NavierStokesModule, 2>
 ns1DPvtProvider("NavierStokes1DPvt");

//////////////////////////////////////////////////////////////////////////////

NavierStokes1DPvt::NavierStokes1DPvt(const std::string& name, Common::SafePtr<Framework::PhysicalModelImpl> model) :
  NavierStokes1DVarSet(name, model),
  _eulerModel(model->getConvectiveTerm().d_castTo<EulerTerm>())
{
  vector<std::string> names(3);
  names[0] = (!_eulerModel->isIncompressible()) ? "p" : "dp";
  names[1] = "v";
  names[2] = "T";
  setVarNames(names);
}
      
//////////////////////////////////////////////////////////////////////////////

NavierStokes1DPvt::~NavierStokes1DPvt()
{
}

//////////////////////////////////////////////////////////////////////////////

void NavierStokes1DPvt::setGradientVars(const std::vector<RealVector*>& states,
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

void NavierStokes1DPvt::setGradientVarGradients(const vector<RealVector*>& states,
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

void NavierStokes1DPvt::setStateGradients(const vector<RealVector*>& states,
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

CFreal NavierStokes1DPvt::getDynViscosity(const RealVector& state, const vector<RealVector*>& gradients)
{
  const CFreal Tdim = _eulerModel->getTempRef()*state[2];
  const CFreal pdim = _eulerModel->getPressRef()*_eulerModel->getPressureFromState(state[0]);
  return getModel().getDynViscosityDim(pdim, Tdim)/
    (getModel().getReferencePhysicalData())[NSTerm::MU];
}

//////////////////////////////////////////////////////////////////////////////

CFreal NavierStokes1DPvt::getDensity(const RealVector& state)
{
 const CFreal p = _eulerModel->getPressureFromState(state[0]);
 return _eulerModel->getDensity(p, state[2]);
}

//////////////////////////////////////////////////////////////////////////////

void NavierStokes1DPvt::setGradientState(const RealVector& state)
{
  cf_assert(_gradState.size() == state.size());
  _gradState = state;
}

//////////////////////////////////////////////////////////////////////////////

void NavierStokes1DPvt::computeFluxJacobian(const RealVector& state,
					     const RealVector& gradientJacob,
					     const RealVector& normal,
					     const CFreal& radius,
					     RealMatrix& fluxJacob)
{
  const CFreal nx = normal[XX];
  
  RealVector& nsData = getModel().getPhysicalData();
  const CFreal mu = nsData[NSTerm::MU];
  const CFreal tc = nsData[NSTerm::LAMBDA];
  
  const CFreal twoThird = 2./3.;
  const CFreal tauxx_uL = twoThird*mu*gradientJacob[XX]; 
  const CFreal qx_T = -tc*gradientJacob[XX];
  const CFreal avU = state[1];
  const CFreal coeffTau = getModel().getCoeffTau();
  const CFreal coeffQ = getModel().getCoeffQ();
  
  fluxJacob(1,1) = coeffTau*tauxx_uL*nx;
  
  fluxJacob(2,1) = fluxJacob(1,1)*avU;
  fluxJacob(2,2) = -coeffQ*qx_T*nx;

}

//////////////////////////////////////////////////////////////////////////////

} // namespace NavierStokes

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
