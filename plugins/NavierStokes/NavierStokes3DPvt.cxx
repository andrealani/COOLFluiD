#include "NavierStokes/NavierStokes.hh"
#include "NavierStokes3DPvt.hh"
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

Environment::ObjectProvider<NavierStokes3DPvt, DiffusiveVarSet, NavierStokesModule, 2> 
ns3DPvtProvider("NavierStokes3DPvt");

Environment::ObjectProvider<NavierStokes3DPvt, DiffusiveVarSet, NavierStokesModule, 2> 
ns3DRotationPvtProvider("NavierStokes3DRotationPvt");

//////////////////////////////////////////////////////////////////////////////

NavierStokes3DPvt::NavierStokes3DPvt(const std::string& name, Common::SafePtr<Framework::PhysicalModelImpl> model) :
  NavierStokes3DVarSet(name, model),
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

NavierStokes3DPvt::~NavierStokes3DPvt()
{
}

//////////////////////////////////////////////////////////////////////////////

void NavierStokes3DPvt::setGradientVars(const vector<RealVector*>& states,
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

void NavierStokes3DPvt::setGradientVarGradients(const vector<RealVector*>& states,
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

void NavierStokes3DPvt::setStateGradients(const vector<RealVector*>& states,
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

CFreal NavierStokes3DPvt::getDynViscosity(const RealVector& state, const vector<RealVector*>& gradients)
{
  const CFreal Tdim = _eulerModel->getTempRef()*state[4];
  const CFreal pdim = _eulerModel->getPressRef()*_eulerModel->getPressureFromState(state[0]);
  return getModel().getDynViscosityDim(pdim, Tdim)/
    (getModel().getReferencePhysicalData())[NSTerm::MU];
}

//////////////////////////////////////////////////////////////////////////////

CFreal NavierStokes3DPvt::getDensity(const RealVector& state)
{
 const CFreal p = _eulerModel->getPressureFromState(state[0]);
 return _eulerModel->getDensity(p, state[4]);
}

//////////////////////////////////////////////////////////////////////////////

void NavierStokes3DPvt::setGradientState(const RealVector& state)
{
  cf_assert(_gradState.size() == state.size());
  _gradState = state;
}

//////////////////////////////////////////////////////////////////////////////

void NavierStokes3DPvt::computeFluxJacobian(const RealVector& state,
					    const RealVector& gradientJacob,
					    const RealVector& normal,
					    const CFreal& radius,
					    RealMatrix& fluxJacob)
{
  const CFreal nx = normal[XX];
  const CFreal ny = normal[YY];
  const CFreal nz = normal[ZZ];
  
  RealVector& nsData = getModel().getPhysicalData();
  const CFreal mu = nsData[NSTerm::MU];
  const CFreal tc = nsData[NSTerm::LAMBDA];
  
  const CFreal fourThird = 4./3.;
  const CFreal twoThird = 2./3.;

  const CFreal tauxx_uL = twoThird*mu*gradientJacob[XX]; 
  const CFreal tauxx_vL = -fourThird*mu*gradientJacob[YY];
  const CFreal tauxx_wL = -fourThird*mu*gradientJacob[ZZ];
  
  const CFreal tauyy_uL = -fourThird*mu*gradientJacob[XX];
  const CFreal tauyy_vL = twoThird*mu*gradientJacob[YY];
  const CFreal tauyy_wL = -fourThird*mu*gradientJacob[ZZ];
  
  const CFreal tauzz_uL = -fourThird*mu*gradientJacob[XX];
  const CFreal tauzz_vL = -fourThird*mu*gradientJacob[YY];
  const CFreal tauzz_wL = twoThird*mu*gradientJacob[ZZ];
  
  const CFreal tauxy_uL = mu*gradientJacob[YY]; 
  const CFreal tauxy_vL = mu*gradientJacob[XX];
  
  const CFreal tauxz_uL = mu*gradientJacob[ZZ]; 
  const CFreal tauxz_wL = mu*gradientJacob[XX];
  
  const CFreal tauyz_vL = mu*gradientJacob[ZZ]; 
  const CFreal tauyz_wL = mu*gradientJacob[YY];
  
  const CFreal qx_T = -tc*gradientJacob[XX];
  const CFreal qy_T = -tc*gradientJacob[YY];
  const CFreal qz_T = -tc*gradientJacob[ZZ];
  
  const CFreal avU = state[1];
  const CFreal avV = state[2];
  const CFreal avW = state[3];
  
  const CFreal coeffTau = getModel().getCoeffTau();
  const CFreal coeffQ = getModel().getCoeffQ();
  
  fluxJacob(1,1) = coeffTau*(tauxx_uL*nx + tauxy_uL*ny + tauxz_uL*nz);
  fluxJacob(1,2) = coeffTau*(tauxx_vL*nx + tauxy_vL*ny);
  fluxJacob(1,3) = coeffTau*(tauxx_wL*nx + tauxz_wL*nz);
  
  fluxJacob(2,1) = coeffTau*(tauxy_uL*nx + tauyy_uL*ny);
  fluxJacob(2,2) = coeffTau*(tauxy_vL*nx + tauyy_vL*ny + tauyz_vL*nz);
  fluxJacob(2,3) = coeffTau*(tauyy_wL*ny + tauyz_wL*nz);
  
  fluxJacob(3,1) = coeffTau*(tauxz_uL*nx + tauzz_uL*nz);
  fluxJacob(3,2) = coeffTau*(tauyz_vL*ny + tauzz_vL*nz);
  fluxJacob(3,3) = coeffTau*(tauxz_wL*nx + tauyz_wL*ny + tauzz_wL*nz);
  
  fluxJacob(4,1) = fluxJacob(1,1)*avU + fluxJacob(2,1)*avV + fluxJacob(3,1)*avW;
  fluxJacob(4,2) = fluxJacob(1,2)*avU + fluxJacob(2,2)*avV + fluxJacob(3,2)*avW;
  fluxJacob(4,3) = fluxJacob(1,3)*avU + fluxJacob(2,3)*avV + fluxJacob(3,3)*avW;
  fluxJacob(4,4) = -coeffQ*(qx_T*nx + qy_T*ny + qz_T*nz);
}

//////////////////////////////////////////////////////////////////////////////

void NavierStokes3DPvt::setComposition(const RealVector& state,
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
    if (isPerturb && (iVar == 1 || iVar == 2 || iVar == 3)) {
      useBackUpValues(true);
    }
    else if (isPerturb && (iVar == 0 || iVar == 4)) {
      library->resetComposition(_tempX);
    }
    else if (!isPerturb) {
      CFreal Tdim = _eulerModel->getTempRef()*state[4];
      CFreal pdim = _eulerModel->getPressureFromState(state[0])*
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
