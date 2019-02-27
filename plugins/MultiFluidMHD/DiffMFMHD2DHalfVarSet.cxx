#include "DiffMFMHD2DHalfVarSet.hh"
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

void DiffMFMHD2DHalfVarSet::setup()
{
  DiffMFMHDVarSet::setup();
}

//////////////////////////////////////////////////////////////////////////////
      
RealMatrix& DiffMFMHD2DHalfVarSet::getFlux(const RealVector& state,
					   const vector<RealVector*>& gradients,
					   const CFreal& radius)
{
  setGradientState(state);

  // dummy normal here
  computeTransportProperties(state, gradients, _normal);
  computeStressTensor(state, gradients, radius);
  
  const CFuint endEM = 8;
  for (CFuint i = 0; i < endEM; ++i){
    _fluxVec(i,XX) = _fluxVec(i,YY) = _fluxVec(i,ZZ) = 0.;
  }
  
  const CFuint nbSpecies = getModel().getNbSpecies();    
  for (CFuint i = 0; i < nbSpecies; ++i) {  
    const RealVector& gradT = *gradients[m_TID[i]];
    const CFreal avU = _gradState[m_uID[i]];
    const CFreal avV = _gradState[m_vID[i]];
    const CFreal avW = _gradState[m_vID[i]];
    const CFreal tauXX = m_tau[i](XX, XX);
    const CFreal tauXY = m_tau[i](XX, YY);
    const CFreal tauXZ = m_tau[i](XX, ZZ); 
    const CFreal tauYY = m_tau[i](YY, YY);
    const CFreal tauYZ = m_tau[i](YY, ZZ);
    const CFreal tauZZ = m_tau[i](ZZ, ZZ);
    
    if(getModel().isBraginskii()) {
      computeHeatFluxBraginskii(gradients, i);
    }
    else {
      computeHeatFluxScalar(gradients, i);
    }
    
    _fluxVec(m_uID[i],XX) = tauXX;
    _fluxVec(m_uID[i],YY) = tauXY;
    _fluxVec(m_uID[i],ZZ) = tauXZ; //YM: this term should be zero in 2.5D
    
    _fluxVec(m_vID[i],XX) = tauXY;
    _fluxVec(m_vID[i],YY) = tauYY;
    _fluxVec(m_vID[i],ZZ) = tauYZ; //YM: this term should be zero in 2.5D
    
    _fluxVec(m_wID[i],XX) = tauXZ;
    _fluxVec(m_wID[i],YY) = tauYZ;
    _fluxVec(m_wID[i],ZZ) = tauZZ; //YM: this term should be zero in 2.5D	  
    
    _fluxVec(m_TID[i],XX) = tauXX*avU + tauXY*avV + tauXZ*avW - _qFluxVect[i][XX];
    _fluxVec(m_TID[i],YY) = tauXY*avU + tauYY*avV + tauYZ*avW - _qFluxVect[i][YY];
    _fluxVec(m_TID[i],ZZ) = tauXZ*avU + tauYZ*avV + tauZZ*avW - _qFluxVect[i][ZZ]; //YM: this term should be zero in 2.5D
  }

  return _fluxVec;
}

//////////////////////////////////////////////////////////////////////////////

void DiffMFMHD2DHalfVarSet::getAxiSourceTerm(const RealVector& physicalData,
					    const RealVector& state,
					    const vector<RealVector*>& gradients,
					    const CFreal& radius,
					    RealVector& source)
{
  throw Common::NotImplementedException (FromHere(),"DiffMFMHD2DHalfVarSet::getAxiSourceTerm()");
}

//////////////////////////////////////////////////////////////////////////////

RealVector& DiffMFMHD2DHalfVarSet::getFlux(const RealVector& state,
					  const vector<RealVector*>& gradients,
					  const RealVector& normal,
					  const CFreal& radius)
{  
  setGradientState(state);
  computeTransportProperties(state, gradients, normal);
  computeStressTensor(state, gradients, radius);

  const CFuint endEM = 8;
  for (CFuint i = 0; i < endEM; ++i) {_flux[i] = 0.;}
  
  const CFreal normalzz = 0.; // the computational domain is 2D
  const CFuint nbSpecies = getModel().getNbSpecies();  
  for (CFuint i = 0; i < nbSpecies; ++i) {
    _flux[m_uID[i]] =
      m_tau[i](XX,XX)*normal[XX] + m_tau[i](XX,YY)*normal[YY] + m_tau[i](XX,ZZ)*normalzz;
    _flux[m_vID[i]] =
      m_tau[i](YY,XX)*normal[XX] + m_tau[i](YY,YY)*normal[YY] + m_tau[i](YY,ZZ)*normalzz;
    _flux[m_wID[i]] =
      m_tau[i](ZZ,XX)*normal[XX] + m_tau[i](YY,ZZ)*normal[YY] + m_tau[i](ZZ,ZZ)*normalzz;
    
    if (getModel().isBraginskii()) {
      computeHeatFluxBraginskii(gradients, i);
    }
    else {
      computeHeatFluxScalar(gradients, i);
    }
    _qFlux[i] = MathFunctions::innerProd(_qFluxVect[i], normal);
    
    _flux[m_TID[i]] =
      (m_tau[i](XX,XX)*_gradState[m_uID[i]] + m_tau[i](XX,YY)*_gradState[m_vID[i]] + m_tau[i](XX,ZZ)*_gradState[m_wID[i]])*normal[XX] + 
      (m_tau[i](XX,YY)*_gradState[m_uID[i]] + m_tau[i](YY,YY)*_gradState[m_vID[i]] + m_tau[i](YY,ZZ)*_gradState[m_wID[i]])*normal[YY] +
      (m_tau[i](XX,ZZ)*_gradState[m_uID[i]] + m_tau[i](YY,ZZ)*_gradState[m_vID[i]] + m_tau[i](ZZ,ZZ)*_gradState[m_wID[i]])*normalzz - _qFlux[i];
  }
  
  return _flux;
}
      
//////////////////////////////////////////////////////////////////////////////

void DiffMFMHD2DHalfVarSet::computeHeatFluxBraginskii(const vector<RealVector*>& gradients,
						      const CFuint i)
{
  const CFuint nbSpecies = getModel().getNbSpecies();  
  if (nbSpecies == 2){ // Model with plasma + neutrals
    const RealVector& diffMFMHDData = getModel().getPhysicalData();
    const RealVector& gradTi = *gradients[m_TID[i]];
    const CFuint endVisc = 2;
    
    if(i == 0){			//plasma
      const RealVector& gradTi = *gradients[m_TID[0]];
      const CFreal Tx = gradTi[XX];
      const CFreal Ty = gradTi[YY];	
      
      _qFluxVect[0][XX] = -((diffMFMHDData[endVisc + 0] + diffMFMHDData[endVisc + 3])*Tx + (diffMFMHDData[endVisc + 1] + diffMFMHDData[endVisc + 4])*Ty);
      _qFluxVect[0][YY] = -((diffMFMHDData[endVisc + 1] + diffMFMHDData[endVisc + 4])*Tx + (diffMFMHDData[endVisc + 2] + diffMFMHDData[endVisc + 5])*Ty);
    }
    else{				//neutrals
      _qFluxVect[i] = -diffMFMHDData[endVisc + 7]*gradTi;
    }
  }
}
    
//////////////////////////////////////////////////////////////////////////////

    } // namespace MultiFluidMHD

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
