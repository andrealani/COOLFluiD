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
  const CFuint nbSpecies = getModel().getNbSpecies();
  const CFuint dim = PhysicalModelStack::getActive()->getDim();

  
  _qFluxVect.resize(nbSpecies);
  _qFlux.resize(nbSpecies);
  
  for (CFuint i = 0; i < nbSpecies; ++i) {
    _qFluxVect[i].resize(dim);
  }  
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
  const CFuint nbSpecies = getModel().getNbSpecies();    

  RealVector& DiffMFMHDData = getModel().getPhysicalData();
  
  const EquationSubSysDescriptor& eqSS = PhysicalModelStack::getActive()->getEquationSubSysDescriptor();
  const CFuint totalNbEqs = PhysicalModelStack::getActive()->getNbEq();	// Maxwell's Eqs.+ Multifluid NavierStokes Eqs.
  const CFuint nbEqs = eqSS.getNbEqsSS();  				
  const CFuint iEqSS = eqSS.getEqSS();  
  
  if (nbEqs == totalNbEqs || iEqSS == 0) {
    const CFuint endEM = 8;
    for (CFuint i = 0; i < endEM; ++i){
     _fluxVec(i,XX) = 0;
     _fluxVec(i,YY) = 0;
    }
  }
  
  if (nbEqs == totalNbEqs || iEqSS == 1) {
    //Model with Braginskii
    if(getModel().isBraginskii() == true) {
      if(nbSpecies == 2){
	for (CFuint i = 0; i < nbSpecies; ++i) {  
	  
	  const RealVector& gradT = *gradients[m_TID[i]];
	  const CFreal avU = _gradState[m_uID[i]];
	  const CFreal avV = _gradState[m_vID[i]];
	  const CFreal tauXX = m_tau[i](XX, XX);
	  const CFreal tauXY = m_tau[i](XX, YY);
	  const CFreal tauYY = m_tau[i](YY, YY);
	  const CFuint endVisc = 2;
	  if(i == 0){			//plasma
	    const RealVector& gradTi = *gradients[m_TID[0]];
	    const CFreal Tx = gradTi[XX];
	    const CFreal Ty = gradTi[YY];	
	    
	    _qFluxVect[0][XX] = -((DiffMFMHDData[endVisc + 0] + DiffMFMHDData[endVisc + 3])*Tx + (DiffMFMHDData[endVisc + 1] + DiffMFMHDData[endVisc + 4])*Ty);
	    _qFluxVect[0][YY] = -((DiffMFMHDData[endVisc + 1] + DiffMFMHDData[endVisc + 4])*Tx + (DiffMFMHDData[endVisc + 2] + DiffMFMHDData[endVisc + 5])*Ty);
	  }
	  else{				//neutrals
	    _qFluxVect[i] = -DiffMFMHDData[endVisc + 7]*gradT;
	  }
	
	  _fluxVec(m_uID[i],XX) = tauXX;
	  _fluxVec(m_uID[i],YY) = tauXY;
	
	  _fluxVec(m_vID[i],XX) = tauXY;
	  _fluxVec(m_vID[i],YY) = tauYY;
	
	  _fluxVec(m_TID[i],XX) = tauXX*avU + tauXY*avV - _qFluxVect[i][XX];
	  _fluxVec(m_TID[i],YY) = tauXY*avU + tauYY*avV - _qFluxVect[i][YY];
	}  
      }
    }
    else{		// Default case with constant viscosity
      for (CFuint i = 0; i < nbSpecies; ++i) {  
	
	const RealVector& gradT = *gradients[m_TID[i]];
	const CFreal avU = _gradState[m_uID[i]];
	const CFreal avV = _gradState[m_vID[i]];
	const CFreal tauXX = m_tau[i](XX, XX);
	const CFreal tauXY = m_tau[i](XX, YY);
	const CFreal tauYY = m_tau[i](YY, YY);
	_qFluxVect[i] = -DiffMFMHDData[nbSpecies + i]*gradT;
      
	_fluxVec(m_uID[i],XX) = tauXX;
	_fluxVec(m_uID[i],YY) = tauXY;
      
	_fluxVec(m_vID[i],XX) = tauXY;
	_fluxVec(m_vID[i],YY) = tauYY;
      
	_fluxVec(m_TID[i],XX) = tauXX*avU + tauXY*avV - _qFluxVect[i][XX];
	_fluxVec(m_TID[i],YY) = tauXY*avU + tauYY*avV - _qFluxVect[i][YY];
      }
    }
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
  const CFuint nbSpecies = getModel().getNbSpecies();  
  
  const EquationSubSysDescriptor& eqSS = PhysicalModelStack::getActive()->getEquationSubSysDescriptor();
  const CFuint totalNbEqs = PhysicalModelStack::getActive()->getNbEq();	// Maxwell's Eqs.+ Multifluid NavierStokes Eqs.
  const CFuint nbEqs = eqSS.getNbEqsSS();  				
  const CFuint iEqSS = eqSS.getEqSS();  
  if (nbEqs == totalNbEqs || iEqSS == 0) {
    const CFuint endEM = 8;
    for (CFuint i = 0; i < endEM; ++i){
     _flux[i] = 0;
    }
  }
  if (nbEqs == totalNbEqs || iEqSS == 1) {
    //Model with Braginskii
    if(getModel().isBraginskii() == true) {
      if(nbSpecies == 2){		// Model with plasma + neutrals
	for (CFuint i = 0; i < nbSpecies; ++i) {
	  const RealVector& DiffMFMHDData = getModel().getPhysicalData();
	  
	  //const CFuint dim = PhysicalModelStack::getActive()->getDim();

	  _flux[m_uID[i]] = m_tau[i](XX,XX)*normal[XX] + m_tau[i](XX,YY)*normal[YY];
	  _flux[m_vID[i]] = m_tau[i](YY,XX)*normal[XX] + m_tau[i](YY,YY)*normal[YY];
	  const CFuint endVisc = 2;  
	  
	  if(i == 0){
	    const RealVector& gradTi = *gradients[m_TID[0]];
	    const CFreal Tx = gradTi[XX];
	    const CFreal Ty = gradTi[YY];	
	    
	    CFreal qiX = -((DiffMFMHDData[endVisc + 0] + DiffMFMHDData[endVisc + 3])*Tx + (DiffMFMHDData[endVisc + 1] + DiffMFMHDData[endVisc + 4])*Ty);
	    CFreal qiY = -((DiffMFMHDData[endVisc + 1] + DiffMFMHDData[endVisc + 4])*Tx + (DiffMFMHDData[endVisc + 2] + DiffMFMHDData[endVisc + 5])*Ty);
	    _qFlux[0] = qiX*normal[XX] + qiY*normal[YY];
	  }
	  else{		//neutrals
	    _qFlux[i] = -DiffMFMHDData[endVisc + 7]*MathFunctions::innerProd(*gradients[m_TID[i]], normal);	
	  }  
	  // AL: old implementation, slower 
	  _flux[m_TID[i]] = (m_tau[i](XX,XX)*_gradState[m_uID[i]] + m_tau[i](XX,YY)*_gradState[m_vID[i]])*normal[XX] + 
	    (m_tau[i](XX,YY)*_gradState[m_uID[i]] + m_tau[i](YY,YY)*_gradState[m_vID[i]])*normal[YY] - _qFlux[i];
	    
	}
      }
    }
    else{		// Default case with constant viscosity and consuctivity
      for (CFuint i = 0; i < nbSpecies; ++i) {
	const RealVector& DiffMFMHDData = getModel().getPhysicalData();
	const CFreal lambda = DiffMFMHDData[nbSpecies + i];		//WATCH OUT: only reads one value of the thermal conductivity
	//cout <<"DiffMFMHD2DVarSet::getFlux--> species = " << i <<"\n";
	//cout <<"DiffMFMHD2DVarSet::getFlux--> lambda  = " << lambda <<"\n";
	
	_qFlux[i] = -lambda*MathFunctions::innerProd(*gradients[m_TID[i]], normal);/*-getModel().getCoeffQ()*DiffMFMHDData[DiffMFMHDTerm::LAMBDA]*
				  MathFunctions::innerProd(*gradients[m_TID[i]], normal);*/ // to be changed
      
	// const CFuint dim = PhysicalModelStack::getActive()->getDim();
	///OLD IMPLEMENTATION
	// _flux.slice(m_uID[i],dim) = m_tau[i]*normal;
	///New one (giving same result)
	_flux[m_uID[i]] = m_tau[i](XX,XX)*normal[XX] + m_tau[i](XX,YY)*normal[YY];
	_flux[m_vID[i]] = m_tau[i](YY,XX)*normal[XX] + m_tau[i](YY,YY)*normal[YY];
	
	//_flux[m_TID[i]] = MathFunctions::innerProd(_flux.slice(m_uID[i],dim), _gradState.slice(m_uID[i],dim)) - qFlux;
      
	//AL: old implementation, slower 
	_flux[m_TID[i]] = (m_tau[i](XX,XX)*_gradState[m_uID[i]] + m_tau[i](XX,YY)*_gradState[m_vID[i]])*normal[XX] + 
	  (m_tau[i](XX,YY)*_gradState[m_uID[i]] + m_tau[i](YY,YY)*_gradState[m_vID[i]])*normal[YY] - _qFlux[i];
      }
    }
  }
  return _flux;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace MultiFluidMHD

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
