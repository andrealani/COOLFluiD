#include "NavierStokes/NavierStokes.hh"
#include "NavierStokes2DPuvt.hh"
#include "Environment/ObjectProvider.hh"
#include "NavierStokes/EulerTerm.hh"
#include "Framework/PhysicalChemicalLibrary.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace NavierStokes {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<NavierStokes2DPuvt, DiffusiveVarSet, NavierStokesModule, 2>
 ns2DPuvtProvider("NavierStokes2DPuvt");

//////////////////////////////////////////////////////////////////////////////

NavierStokes2DPuvt::NavierStokes2DPuvt(const std::string& name, Common::SafePtr<Framework::PhysicalModelImpl> model) :
  NavierStokes2DVarSet(name, model),
  _eulerModel(model->getConvectiveTerm().d_castTo<EulerTerm>()),
  _tempX()
{
  vector<std::string> names(4);
  names[0] = (!_eulerModel->isIncompressible()) ? "p" : "dp";
  names[1] = "u";
  names[2] = "v";
  names[3] = "T";
  setVarNames(names);
}

//////////////////////////////////////////////////////////////////////////////

NavierStokes2DPuvt::~NavierStokes2DPuvt()
{
}

//////////////////////////////////////////////////////////////////////////////

void NavierStokes2DPuvt::setGradientVars(const std::vector<RealVector*>& states,
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

void NavierStokes2DPuvt::setGradientVarGradients(const vector<RealVector*>& states,
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

void NavierStokes2DPuvt::setStateGradients(const vector<RealVector*>& states,
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

CFreal NavierStokes2DPuvt::getDynViscosity(const RealVector& state, const vector<RealVector*>& gradients)
{
  const CFreal Tdim = _eulerModel->getTempRef()*state[3];
  const CFreal pdim = _eulerModel->getPressRef()*_eulerModel->getPressureFromState(state[0]);
  return getModel().getDynViscosityDim(pdim, Tdim)/
    (getModel().getReferencePhysicalData())[NSTerm::MU];
}

//////////////////////////////////////////////////////////////////////////////

CFreal NavierStokes2DPuvt::getDensity(const RealVector& state)
{
 const CFreal p = _eulerModel->getPressureFromState(state[0]);
 return _eulerModel->getDensity(p, state[3]);
}

//////////////////////////////////////////////////////////////////////////////

void NavierStokes2DPuvt::setGradientState(const RealVector& state)
{
  cf_assert(_gradState.size() == state.size());
  _gradState = state;
}

//////////////////////////////////////////////////////////////////////////////

void NavierStokes2DPuvt::computeFluxJacobian(const RealVector& state,
					     const RealVector& gradientJacob,
					     const RealVector& normal,
					     const CFreal& radius,
					     RealMatrix& fluxJacob)
{
  const CFreal nx = normal[XX];
  const CFreal ny = normal[YY];
  
  RealVector& nsData = getModel().getPhysicalData();
  const CFreal mu = nsData[NSTerm::MU];
  const CFreal tc = nsData[NSTerm::LAMBDA];
  
  const CFreal fourThird = 4./3.;
  const CFreal twoThird = 2./3.;

  const CFreal tauxx_uL = twoThird*mu*gradientJacob[XX]; 
  const CFreal tauxx_vL = -fourThird*mu*gradientJacob[YY]; 
  const CFreal tauxy_uL = mu*gradientJacob[YY]; 
  const CFreal tauxy_vL = mu*gradientJacob[XX];
  const CFreal tauyy_uL = -fourThird*mu*gradientJacob[XX];
  const CFreal tauyy_vL = twoThird*mu*gradientJacob[YY];
  const CFreal qx_T = -tc*gradientJacob[XX];
  const CFreal qy_T = -tc*gradientJacob[YY];
  const CFreal avU = state[1];
  const CFreal avV = state[2];
  const CFreal coeffTau = getModel().getCoeffTau();
  const CFreal coeffQ = getModel().getCoeffQ();
  
  fluxJacob(1,1) = coeffTau*(tauxx_uL*nx + tauxy_uL*ny);
  fluxJacob(1,2) = coeffTau*(tauxx_vL*nx + tauxy_vL*ny);
  
  fluxJacob(2,1) = coeffTau*(tauxy_uL*nx + tauyy_uL*ny);
  fluxJacob(2,2) = coeffTau*(tauxy_vL*nx + tauyy_vL*ny);
  
  fluxJacob(3,1) = fluxJacob(1,1)*avU + fluxJacob(2,1)*avV;
  fluxJacob(3,2) = fluxJacob(1,2)*avU + fluxJacob(2,2)*avV;
  fluxJacob(3,3) = -coeffQ*(qx_T*nx + qy_T*ny);
}

//////////////////////////////////////////////////////////////////////////////

void NavierStokes2DPuvt::setComposition(const RealVector& state,
					const bool isPerturb,
					const CFuint iVar)
{
  // this is to avoid useless expensive re-computations
  useBackUpValues(false);
  setBackUpValues(false);
  
  SafePtr<PhysicalChemicalLibrary> library = 
    PhysicalModelStack::getActive()->getImplementor()->
    getPhysicalPropertyLibrary<PhysicalChemicalLibrary>();
  if (!library.isNotNull()) {
    if (isPerturb && iVar != _TID) {
      useBackUpValues(true);
    }
    else if (!isPerturb) {
      setBackUpValues(true); 
    }
  }
  else {
    _tempX.resize(library->getNbSpecies());
    if (isPerturb && (iVar == 1 || iVar == 2)) {
      useBackUpValues(true);
    }
    else if (isPerturb && (iVar == 0 || iVar == 3)) {
      /// @todo test if it is possible to use it like this
      //      	library->resetComposition(_tempX);
      
      CFreal Tdim   = _eulerModel->getTempRef()*state[3];
      CFreal pdim   =  _eulerModel->getPressureFromState(state[0])*
	(_eulerModel->getReferencePhysicalData())[EulerTerm::P];
      library->setComposition(Tdim,pdim);
    }
    else if (!isPerturb) {
      CFreal Tdim   = _eulerModel->getTempRef()*state[3];
      CFreal pdim   =  _eulerModel->getPressureFromState(state[0])*
	(_eulerModel->getReferencePhysicalData())[EulerTerm::P];
      library->setComposition(Tdim,pdim,&_tempX);
      // set and store the back up values only if an unperturbed flux
      // is computed
      setBackUpValues(true); 
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

} // namespace NavierStokes

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
