#include "LTE.hh"
#include "NavierStokes3DLTEDemixPvt.hh"
#include "Environment/ObjectProvider.hh"
#include "Framework/PhysicalChemicalLibrary.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace LTE {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<NavierStokes3DLTEDemixPvt, DiffusiveVarSet,
	       LTEModule, 2>
ns3DLTEDemixPvtProvider("NavierStokes3DLTEDemixPvt");

//////////////////////////////////////////////////////////////////////////////

NavierStokes3DLTEDemixPvt::NavierStokes3DLTEDemixPvt
(const std::string& name, SafePtr<PhysicalModelImpl> model) :
  NavierStokesLTEDemixVarSet(name, model),
  _library(CFNULL),
  _eulerModel(model->getConvectiveTerm().d_castTo<EulerLTEDemixTerm>()),
  _tempX(),
  _ye(),
  _lambdaEL(),
  _eldifcoef(),
  _eltdifcoef()
{
  const CFuint nbElements = _eulerModel->getNbScalarVars(0);

  vector<std::string> names(5 + nbElements);
  names[0] = (!_eulerModel->isIncompressible()) ? "p" : "dp";
  names[1] = "u";
  names[2] = "v";
  names[3] = "w";
  names[4] = "T";
  
  // Names for the elemental fractions
  for (CFuint ie = 0; ie < nbElements; ++ie) {
    names[5 + ie] = "ye" + StringOps::to_str(ie);
  }

  setVarNames(names);
}

//////////////////////////////////////////////////////////////////////////////

NavierStokes3DLTEDemixPvt::~NavierStokes3DLTEDemixPvt()
{
}

//////////////////////////////////////////////////////////////////////////////

void NavierStokes3DLTEDemixPvt::setGradientVars(const vector<RealVector*>& states,
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

void NavierStokes3DLTEDemixPvt::setGradientVarGradients(const vector<RealVector*>& states,
                                                        const vector< vector<RealVector*> >& stateGradients,
                                                        vector< vector<RealVector*> >& gradVarGradients,
                                                        const CFuint stateSize)
{
  throw Common::NotImplementedException(FromHere(),"setGradientVarGradients");
}

//////////////////////////////////////////////////////////////////////////////

void NavierStokes3DLTEDemixPvt::setStateGradients(const vector<RealVector*>& states,
                                                  const vector< vector<RealVector*> >& gradVarGradients,
                                                  vector< vector<RealVector*> >& stateGradients,
                                                  const CFuint stateSize)
{
  throw Common::NotImplementedException(FromHere(),"setStateGradients");
}

//////////////////////////////////////////////////////////////////////////////

CFreal NavierStokes3DLTEDemixPvt::getDynViscosity(const RealVector& state, const vector<RealVector*>& gradients)
{
  CFreal Tdim = _eulerModel->getTempRef()*state[4];
  CFreal pdim = _eulerModel->getPressureFromState(state[0])*
    (_eulerModel->getReferencePhysicalData())[EulerTerm::P];
  return _library->eta(Tdim, pdim, CFNULL) /
    (getModel().getReferencePhysicalData())[NSTerm::MU];
}

//////////////////////////////////////////////////////////////////////////////

CFreal NavierStokes3DLTEDemixPvt::getDensity(const RealVector& state)
{
  CFreal Tdim = _eulerModel->getTempRef()*state[3];
  CFreal pdim = _eulerModel->getPressureFromState(state[0])*
    (_eulerModel->getReferencePhysicalData())[EulerTerm::P];
  const CFreal rhoRef = (_eulerModel->getReferencePhysicalData())[EulerTerm::RHO];
  return (PhysicalModelStack::getActive()->isAdimensional()) ?
    _library->density(Tdim,pdim)/rhoRef : _library->density(Tdim,pdim);
}

//////////////////////////////////////////////////////////////////////////////

void NavierStokes3DLTEDemixPvt::setup()
{
  //  fluct split has to call setup()
  NavierStokesLTEDemixVarSet::setup();

  _library = PhysicalModelStack::getActive()->getImplementor()->
    getPhysicalPropertyLibrary<PhysicalChemicalLibrary>();
  cf_assert(_library.isNotNull());

  const CFuint nbElements = _eulerModel->getNbScalarVars(0);

  _tempX.resize(_library->getNbSpecies());
  _ye.resize(nbElements);
  _lambdaEL.resize(nbElements);
  _eldifcoef.resize(nbElements, nbElements);
  _eltdifcoef.resize(nbElements);
}

//////////////////////////////////////////////////////////////////////////////

void NavierStokes3DLTEDemixPvt::setComposition(const RealVector& state,
                                               const bool isPerturb,
                                               const CFuint iVar)
{
  // this is to avoid useless expensive re-computations
  useBackUpValues(false);
  setBackUpValues(false);

  if (isPerturb && (iVar == 1 || iVar == 2 || iVar == 3)) {
    useBackUpValues(true);
  }
  else if (isPerturb && (iVar == 0 || iVar >= 4)) {
    _library->resetComposition(_tempX);
  }
  else if (!isPerturb) {
    CFreal Tdim = _eulerModel->getTempRef()*state[4];
    CFreal pdim = _eulerModel->getPressureFromState(state[0])*
      (_eulerModel->getReferencePhysicalData())[EulerTerm::P];
    
    // Set the elemental fractions
    const CFuint nbElements = _eulerModel->getNbScalarVars(0);
    
    for (CFuint ie = 0; ie < nbElements; ++ie) {
      _ye[ie] = state[5 + ie];
    }
    
    _library->setElemFractions(_ye);
    _library->setComposition(Tdim,pdim,&_tempX);
    // set and store the back up values only if an unperturbed flux
    // is computed
    setBackUpValues(true);
  }
}

//////////////////////////////////////////////////////////////////////////////

RealVector& NavierStokes3DLTEDemixPvt::getFlux(const RealVector& values,
                                               const vector<RealVector*>& gradients,
                                               const RealVector& normal,
                                               const CFreal& radius)
{
  const RealVector& gradU = *gradients[1];
  const RealVector& gradV = *gradients[2];
  const RealVector& gradW = *gradients[3];
  const RealVector& gradT = *gradients[4];
  const CFreal twoThirdDivV = 2./3.*(gradU[XX] + gradV[YY] + gradW[ZZ]);

  const CFreal avU = values[1];
  const CFreal avV = values[2];
  const CFreal avW = values[3];

  RealVector& nsData = getModel().getPhysicalData();
  
  CFdouble avTdim = _eulerModel->getTempRef()*values[4];
  CFdouble avpdim = _eulerModel->getPressureFromState(values[0])*
    (_eulerModel->getReferencePhysicalData())[EulerTerm::P];
  
  CFdouble lambdaR = 0.0; // thermal reactive conductivity
  CFdouble lambdaD = 0.0; // thermal demixing conductivity

  // get transport properties
  _library->getTransportCoefs(avTdim,
			      avpdim,
			      lambdaR,
			      lambdaD,
			      _lambdaEL,
			      _eldifcoef,
			      _eltdifcoef);

  // adimensional dynamical viscosity
  if (_useBackUpValues) {
    nsData[NSTerm::MU] = _dynViscCoeff;
    nsData[NSTerm::LAMBDA] = _thermCondCoeff;
  }
  else {
    // adimensional dynamical viscosity
    nsData[NSTerm::MU] = getDynViscosity(values, gradients);

    // adimensional thermal conductivity
    nsData[NSTerm::LAMBDA] =
      (_library->lambdaEQ(avTdim, avpdim) + lambdaR + lambdaD)/
      (getModel().getReferencePhysicalData())[NSTerm::LAMBDA];

    if (_setBackUpValues) {
      _dynViscCoeff = nsData[NSTerm::MU];
      _thermCondCoeff = nsData[NSTerm::LAMBDA];
    }
  }

  const CFreal coeffTauMu = getModel().getCoeffTau()*nsData[NSTerm::MU];
  const CFreal tauXX = coeffTauMu*(2.*gradU[XX] - twoThirdDivV);
  const CFreal tauXY = coeffTauMu*(gradU[YY] + gradV[XX]);
  const CFreal tauYY = coeffTauMu*(2.*gradV[YY] - twoThirdDivV);
  const CFreal tauXZ = coeffTauMu*(gradU[ZZ] + gradW[XX]);
  const CFreal tauYZ = coeffTauMu*(gradW[YY] + gradV[ZZ]);
  const CFreal tauZZ = coeffTauMu*(2.*gradW[ZZ] - twoThirdDivV);
  const CFreal nx = normal[XX];
  const CFreal ny = normal[YY];
  const CFreal nz = normal[ZZ];

  // Elemental heat transfer coefficient term
  const CFuint nbElements = _eulerModel->getNbScalarVars(0);
  CFdouble lambdaelGradYex = 0.0;
  CFdouble lambdaelGradYey = 0.0;
  CFdouble lambdaelGradYez = 0.0;
  for (CFuint ie = 0; ie < nbElements; ++ie) {
    lambdaelGradYex = lambdaelGradYex + _lambdaEL[ie] *
      (*gradients[5 + ie])[XX];
    lambdaelGradYey = lambdaelGradYey + _lambdaEL[ie] *
      (*gradients[5 + ie])[YY];
    lambdaelGradYez = lambdaelGradYez + _lambdaEL[ie] *
      (*gradients[5 + ie])[ZZ];
  }

  const CFreal qFlux = -getModel().getCoeffQ()*nsData[NSTerm::LAMBDA]*
    (gradT[XX]*nx + gradT[YY]*ny + gradT[ZZ]*nz) -
    (lambdaelGradYex*nx + lambdaelGradYey*ny + lambdaelGradYez*nz);

  _flux[1] = tauXX*nx + tauXY*ny + tauXZ*nz;
  _flux[2] = tauXY*nx + tauYY*ny + tauYZ*nz;
  _flux[3] = tauXZ*nx + tauYZ*ny + tauZZ*nz;
  _flux[4] = (tauXX*avU + tauXY*avV + tauXZ*avW)*nx +
    (tauXY*avU + tauYY*avV + tauYZ*avW)*ny +
    (tauXZ*avU + tauYZ*avV + tauZZ*avW)*nz - qFlux;

  // Fluxes for the element continuities
  for (CFuint ie = 0; ie < nbElements; ++ie) {
    // Elemental multicomponent diffusion coefficient terms
    CFdouble rhoDGradYex = 0.0;
    CFdouble rhoDGradYey = 0.0;
    CFdouble rhoDGradYez = 0.0;
    for (CFuint je = 0; je < nbElements; ++je) {
      rhoDGradYex = rhoDGradYex + _eldifcoef(ie, je) * (*gradients[5 + je])[XX];
      rhoDGradYey = rhoDGradYey + _eldifcoef(ie, je) * (*gradients[5 + je])[YY];
      rhoDGradYez = rhoDGradYez + _eldifcoef(ie, je) * (*gradients[5 + je])[ZZ];
    }

    _flux[5 + ie] = (_eltdifcoef[ie] * gradT[XX] + rhoDGradYex) * nx +
                    (_eltdifcoef[ie] * gradT[YY] + rhoDGradYey) * ny +
                    (_eltdifcoef[ie] * gradT[ZZ] + rhoDGradYez) * nz;
  }

  return _flux;
}

//////////////////////////////////////////////////////////////////////////////

RealMatrix& NavierStokes3DLTEDemixPvt::getFlux(const RealVector& values,
                                               const vector<RealVector*>& gradients,
                                               const CFreal& radius)
{
  const RealVector& gradU = *gradients[1];
  const RealVector& gradV = *gradients[2];
  const RealVector& gradW = *gradients[3];
  const RealVector& gradT = *gradients[4];
  const CFreal twoThirdDivV = 2./3.*(gradU[XX] + gradV[YY] + gradW[ZZ]);

  const CFreal avU = values[1];
  const CFreal avV = values[2];
  const CFreal avW = values[3];

  RealVector& nsData = getModel().getPhysicalData();

  CFdouble avTdim = _eulerModel->getTempRef()*values[4];
  CFdouble avpdim = _eulerModel->getPressureFromState(values[0])*
    (_eulerModel->getReferencePhysicalData())[EulerTerm::P];
  
  CFdouble lambdaR = 0.0; // thermal reactive conductivity
  CFdouble lambdaD = 0.0; // thermal demixing conductivity

  // get transport properties
  _library->getTransportCoefs(avTdim,
                              avpdim,
                              lambdaR,
                              lambdaD,
                              _lambdaEL,
                              _eldifcoef,
                              _eltdifcoef);

  // adimensional dynamical viscosity
  if (_useBackUpValues) {
    nsData[NSTerm::MU] = _dynViscCoeff;
    nsData[NSTerm::LAMBDA] = _thermCondCoeff;
  }
  else {
    // adimensional dynamical viscosity
    nsData[NSTerm::MU] = getDynViscosity(values, gradients);

    // adimensional thermal conductivity
    nsData[NSTerm::LAMBDA] =
      (_library->lambdaEQ(avTdim, avpdim) + lambdaR + lambdaD)/
      (getModel().getReferencePhysicalData())[NSTerm::LAMBDA];

    if (_setBackUpValues) {
      _dynViscCoeff = nsData[NSTerm::MU];
      _thermCondCoeff = nsData[NSTerm::LAMBDA];
    }
  }

  const CFreal coeffTauMu = getModel().getCoeffTau()*nsData[NSTerm::MU];
  const CFreal tauXX = coeffTauMu*(2.*gradU[XX] - twoThirdDivV);
  const CFreal tauXY = coeffTauMu*(gradU[YY] + gradV[XX]);
  const CFreal tauYY = coeffTauMu*(2.*gradV[YY] - twoThirdDivV);
  const CFreal tauXZ = coeffTauMu*(gradU[ZZ] + gradW[XX]);
  const CFreal tauYZ = coeffTauMu*(gradW[YY] + gradV[ZZ]);
  const CFreal tauZZ = coeffTauMu*(2.*gradW[ZZ] - twoThirdDivV);

  // Elemental heat transfer coefficient term
  const CFuint nbElements = _eulerModel->getNbScalarVars(0);
  RealVector lambdaelGradYe(0.0,3);
  for (CFuint ie = 0; ie < nbElements; ++ie) {
    lambdaelGradYe += _lambdaEL[ie]*(*gradients[5 + ie]);
  }

  const RealVector qFlux = -getModel().getCoeffQ()*nsData[NSTerm::LAMBDA]*gradT - lambdaelGradYe;

  _fluxVec(1,XX) = tauXX;
  _fluxVec(1,YY) = tauXY;
  _fluxVec(1,ZZ) = tauXZ;

  _fluxVec(2,XX) = tauXY;
  _fluxVec(2,YY) = tauYY;
  _fluxVec(2,ZZ) = tauYZ;

  _fluxVec(3,XX) = tauXZ;
  _fluxVec(3,YY) = tauYZ;
  _fluxVec(3,ZZ) = tauZZ;

  _fluxVec(4,XX) = tauXX*avU + tauXY*avV + tauXZ*avW - qFlux[XX];
  _fluxVec(4,YY) = tauXY*avU + tauYY*avV + tauYZ*avW - qFlux[YY];
  _fluxVec(4,ZZ) = tauXZ*avU + tauYZ*avV + tauZZ*avW - qFlux[ZZ];

  // Fluxes for the element continuities
  for (CFuint ie = 0; ie < nbElements; ++ie) {
    // Elemental multicomponent diffusion coefficient terms
    RealVector rhoDGradYe(0.0,3);
    for (CFuint je = 0; je < nbElements; ++je) {
      rhoDGradYe += _eldifcoef(ie, je)*(*gradients[5 + je]);
    }
    rhoDGradYe += _eltdifcoef[ie]*gradT;
    _fluxVec.setRow(rhoDGradYe,5 + ie);
  }

  return _fluxVec;
}

//////////////////////////////////////////////////////////////////////////////

void NavierStokes3DLTEDemixPvt::setGradientState(const RealVector& state)
{
  cf_assert(_gradState.size() == state.size());
  _gradState = state;
}

//////////////////////////////////////////////////////////////////////////////

CFreal NavierStokes3DLTEDemixPvt::getHeatFlux(const RealVector& state,
					      const std::vector<RealVector*>& gradients,
					      const RealVector& normal)
{
  const RealVector& gradT = *gradients[4];
  CFdouble avTdim = _eulerModel->getTempRef()*state[4];
  CFreal avpdim   = _eulerModel->getPressureFromState(state[0])*
    (_eulerModel->getReferencePhysicalData())[EulerTerm::P];
  CFdouble lambdaR = 0.0; // thermal reactive conductivity
  CFdouble lambdaD = 0.0; // thermal demixing conductivity
  
  // get transport properties
  _library->getTransportCoefs(avTdim,
			      avpdim,
			      lambdaR,
			      lambdaD,
			      _lambdaEL,
			      _eldifcoef,
			      _eltdifcoef);

  // adimensional thermal conductivity
  const CFreal lambda = (_library->lambdaEQ(avTdim, avpdim) + lambdaR + lambdaD)/
    (getModel().getReferencePhysicalData())[NSTerm::LAMBDA];
  
  // Elemental heat transfer coefficient term
  CFdouble lambdaelGradYex = 0.0;
  CFdouble lambdaelGradYey = 0.0;
  CFdouble lambdaelGradYez = 0.0;
  const CFuint nbElements = _eulerModel->getNbScalarVars(0);
  for (CFuint ie = 0; ie < nbElements; ++ie) {
    lambdaelGradYex += _lambdaEL[ie] * (*gradients[5 + ie])[XX];
    lambdaelGradYey += _lambdaEL[ie] * (*gradients[5 + ie])[YY];
    lambdaelGradYez += _lambdaEL[ie] * (*gradients[5 + ie])[ZZ];
  }

  return (-getModel().getCoeffQ()*lambda*(gradT[XX]*normal[XX] +
					  gradT[YY]*normal[YY] +
					  gradT[ZZ]*normal[ZZ])
	  - (lambdaelGradYex*normal[XX] + lambdaelGradYey*normal[YY]
	     + lambdaelGradYez*normal[ZZ]));
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace LTE

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
