#include "DiffMFMHD2DVarSet.hh"
#include "EulerMFMHDTerm.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace MultiFluidMHD {

//////////////////////////////////////////////////////////////////////////////

void DiffMFMHD2DVarSet::setup()
{
  DiffMFMHDVarSet::setup();
}

//////////////////////////////////////////////////////////////////////////////

RealMatrix& DiffMFMHD2DVarSet::getFlux(const RealVector& state,
				       const vector<RealVector*>& gradients,
				       const CFreal& radius)
{
  setGradientState(state);
  // dummy normal here
  computeTransportProperties(state, gradients, _normal);
  computeStressTensor(state, gradients, radius);

  const CFuint endEM = 8;
  for (CFuint i = 0; i < endEM; ++i){
    _fluxVec(i,XX) = _fluxVec(i,YY) = 0.;
  }
  
  const CFuint nbSpecies = getModel().getNbSpecies();    
  for (CFuint i = 0; i < nbSpecies; ++i) {
    const RealVector& gradT = *gradients[m_TID[i]];
    const CFreal avU = _gradState[m_uID[i]];
    const CFreal avV = _gradState[m_vID[i]];
    const CFreal tauXX = m_tau[i](XX, XX);
    const CFreal tauXY = m_tau[i](XX, YY);
    const CFreal tauYY = m_tau[i](YY, YY);

    if(getModel().isBraginskii()) {
      computeHeatFluxBraginskii(gradients, i);
    }
    else {
      computeHeatFluxScalar(gradients, i);
    }
    
    _fluxVec(m_uID[i],XX) = tauXX;
    _fluxVec(m_uID[i],YY) = tauXY;

    _fluxVec(m_vID[i],XX) = tauXY;
    _fluxVec(m_vID[i],YY) = tauYY;

    _fluxVec(m_TID[i],XX) = tauXX*avU + tauXY*avV - _qFluxVect[i][XX];
    _fluxVec(m_TID[i],YY) = tauXY*avU + tauYY*avV - _qFluxVect[i][YY];
  }
  
  return _fluxVec;
}

//////////////////////////////////////////////////////////////////////////////

void DiffMFMHD2DVarSet::getAxiSourceTerm(const RealVector& physicalData,
					 const RealVector& state,
					 const vector<RealVector*>& gradients,
					 const CFreal& radius,
					 RealVector& source)
{
  throw Common::NotImplementedException (FromHere(),"DiffMFMHD2DVarSet::getAxiSourceTerm()");
}

//////////////////////////////////////////////////////////////////////////////

RealVector& DiffMFMHD2DVarSet::getFlux(const RealVector& state,
				       const vector<RealVector*>& gradients,
				       const RealVector& normal,
				       const CFreal& radius)
{  
  setGradientState(state);
  computeTransportProperties(state, gradients, normal);
  computeStressTensor(state, gradients, radius);
  
  const CFuint endEM = 8;
  for (CFuint i = 0; i < endEM; ++i) {_flux[i] = 0;}

  const CFuint nbSpecies = getModel().getNbSpecies();  
  for (CFuint i = 0; i < nbSpecies; ++i) {
    _flux[m_uID[i]] = m_tau[i](XX,XX)*normal[XX] + m_tau[i](XX,YY)*normal[YY];
    _flux[m_vID[i]] = m_tau[i](YY,XX)*normal[XX] + m_tau[i](YY,YY)*normal[YY];
    
    if (getModel().isBraginskii()) {
      computeHeatFluxBraginskii(gradients, i);
    }
    else {
      computeHeatFluxScalar(gradients, i);
    }
    _qFlux[i] = MathFunctions::innerProd(_qFluxVect[i], normal);
    
    _flux[m_TID[i]] =
      (m_tau[i](XX,XX)*_gradState[m_uID[i]] + m_tau[i](XX,YY)*_gradState[m_vID[i]])*normal[XX] +
      (m_tau[i](XX,YY)*_gradState[m_uID[i]] + m_tau[i](YY,YY)*_gradState[m_vID[i]])*normal[YY] - _qFlux[i];
  }
  
  return _flux;
}

//////////////////////////////////////////////////////////////////////////////

void DiffMFMHD2DVarSet::computeHeatFluxBraginskii(const vector<RealVector*>& gradients,
						  const CFuint i)
{
  const CFuint nbSpecies = getModel().getNbSpecies();  
  if (nbSpecies == 2){ // Model with plasma + neutrals
    const RealVector& diffMFMHDData = getModel().getPhysicalData();
    const CFuint endVisc = 2;
    const RealVector& gradTi = *gradients[m_TID[i]];

    if(i == 0) { //plasma
      const RealVector& gradTi = *gradients[m_TID[0]];
      const CFreal Tx = gradTi[XX];
      const CFreal Ty = gradTi[YY];
      
      // AL: double check this
      _qFluxVect[0][XX] = -((diffMFMHDData[endVisc + 0] + diffMFMHDData[endVisc + 3])*Tx + (diffMFMHDData[endVisc + 1] + diffMFMHDData[endVisc + 4])*Ty);
      _qFluxVect[0][YY] = -((diffMFMHDData[endVisc + 1] + diffMFMHDData[endVisc + 4])*Tx + (diffMFMHDData[endVisc + 2] + diffMFMHDData[endVisc + 5])*Ty);
    }
    else { //neutrals
      _qFluxVect[i] = -diffMFMHDData[endVisc + 7]*gradTi;
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

} // namespace MultiFluidMHD

} // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
