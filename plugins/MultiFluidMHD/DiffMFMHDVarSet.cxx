#include "DiffMFMHDVarSet.hh"
#include "EulerMFMHDTerm.hh"
#include "Framework/PhysicalConsts.hh"
#include "MathTools/MathConsts.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace MultiFluidMHD {

//////////////////////////////////////////////////////////////////////////////
       
DiffMFMHDVarSet::DiffMFMHDVarSet(const std::string& name, 
				 SafePtr<PhysicalModelImpl> model) :
  DiffusiveVarSet(name, model),
  m_model(model->getDiffusiveTerm().d_castTo<DiffMFMHDTerm>()),
  m_eulerModel(model->getConvectiveTerm().d_castTo<PTERM>()),
  _twoThird(2./3.),
  _useBackUpValues(false),
  _setBackUpValues(false),
  _wallDistance(0.),
  _qFluxVect(),
  _qFlux(),
  _gradState(),
  _normal(),
  _heatFlux(),
  m_dynViscCoeff(),
  m_thermCondCoeff(),
  m_tau(),
  m_uID(),
  m_vID(),
  m_wID(),
  m_TID(),
  m_lambda()
{ 
}
  
//////////////////////////////////////////////////////////////////////////////
      
DiffMFMHDVarSet::~DiffMFMHDVarSet()
{
}
    
//////////////////////////////////////////////////////////////////////////////
   
void DiffMFMHDVarSet::setup()
{
  DiffusiveVarSet::setup();
  
  const CFuint nbSpecies = getModel().getNbSpecies();
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  
  _qFluxVect.resize(nbSpecies);
  for (CFuint i = 0; i < nbSpecies; ++i) {
    _qFluxVect[i].resize(dim);
  }  
  _qFlux.resize(nbSpecies);
  
  _gradState.resize(PhysicalModelStack::getActive()->getNbEq());
  _normal.resize(dim);
  const CFuint endEM = 8;
  
  _heatFlux.resize(nbSpecies);
  m_tau.resize(nbSpecies);
  m_uID.resize(nbSpecies);
  m_vID.resize(nbSpecies);
  m_wID.resize(nbSpecies);
  m_TID.resize(nbSpecies);
  m_dynViscCoeff.resize(nbSpecies);
  m_currDynViscosity.resize(nbSpecies);
  m_lambda.resize(nbSpecies);
  if(getModel().isBraginskii() == true) {
    if(nbSpecies == 2){
      m_thermCondCoeff.resize(8);
    }
  }
  else{		// Default case with constant viscosity and thermal conductivity
    m_thermCondCoeff.resize(nbSpecies);
  }
  
  
  for (CFuint i = 0; i < nbSpecies; ++i) {
    m_uID[i] = endEM + nbSpecies + i*dim;
    m_vID[i] = endEM + nbSpecies + i*dim + 1;
    m_TID[i] = endEM + (dim + 1)*nbSpecies + i;
   
    if (dim == DIM_3D) {
      m_uID[i] = endEM + nbSpecies + i*dim;
      m_vID[i] = endEM + nbSpecies + i*dim + 1;
      m_wID[i] = endEM + nbSpecies + i*dim + 2;
      m_TID[i] = endEM + nbSpecies + nbSpecies*dim + i;
    }
    m_tau[i].resize(dim,dim);
  }
}

//////////////////////////////////////////////////////////////////////////////

void DiffMFMHDVarSet::computeTransportProperties(const RealVector& state,
						 const std::vector<RealVector*>& gradients,
						 const RealVector& normal)
{  
  RealVector& diffMFMHDData = getModel().getPhysicalData();
  const CFuint nbSpecies = getModel().getNbSpecies();
    
  if (_useBackUpValues || _freezeDiffCoeff) {
    // here you just store precomputed values of the transport properties, this is typically
    // called if the user asks to freeze transport properties when computing the numerical
    // jacobian of the diffusive fluxes
    for (CFuint i = 0; i < nbSpecies; ++i){
      diffMFMHDData[i] = m_dynViscCoeff[i];
      diffMFMHDData[nbSpecies+i] = m_thermCondCoeff[i];
    }
  }
  else {
    // this is the default case, transport properties are computed on the fly
    const RealVector& mu = getDynViscosityVec(state, gradients);
    
    for (CFuint i = 0; i < nbSpecies; ++i) {
      diffMFMHDData[i] = mu[i]; //reads the dimensional value imposed in the options but only takes the max
      diffMFMHDData[nbSpecies+i] = getThermConductivity(state, mu[i]);	//reads the dimensional value imposed in the options but only takes the max
    }
        
    if (_setBackUpValues) {
      for (CFuint i = 0; i < nbSpecies; ++i) {
        m_dynViscCoeff[i] = diffMFMHDData[i];
        m_thermCondCoeff[i] = diffMFMHDData[nbSpecies+i];
	
	CFLog(DEBUG_MAX, "DiffMFMHDVarSet::computeTransportProperties() [" << i <<  "] => "
	      << m_dynViscCoeff[i] << ", "  << m_thermCondCoeff[i] << "\n");
      }
    }
  }
  
  /// Compute solar properties for one fluid
  if (getModel().isSolarTransport1F()) {
    computeSolarProperties1F(state, gradients, diffMFMHDData);
  }
  
  /// Compute Braginskii properties
  if (getModel().isBraginskii()) {
    computeBraginskiiProperties(state, gradients, diffMFMHDData);
  }
  
  /// Increase the value of the viscosity in the extended domain
  if(getModel().isExtendedDomain()){
    computePropertiesInExtendedDomain(state, gradients, diffMFMHDData);
  }
}

//////////////////////////////////////////////////////////////////////////////

void DiffMFMHDVarSet::computeSolarProperties1F(const RealVector& state,
					       const std::vector<RealVector*>& gradients,
					       RealVector& diffMFMHDData)
{
  CFLog(DEBUG_MAX, "DiffMFMHDVarSet::computeSolarTransport1F()\n");
  
  // Peter: here is your part to compute and store the thermal conductivities to be used later
  //        inside DiffMFMHD3DVarSet::computeHeatFluxSolar1F()
  
  // diffMFMHDData[0] = ...; 
  // diffMFMHDData[1] = ...; 
}
      
//////////////////////////////////////////////////////////////////////////////

void DiffMFMHDVarSet::computeBraginskiiProperties(const RealVector& state,
						  const std::vector<RealVector*>& gradients,
						  RealVector& diffMFMHDData)
{
  const CFuint nbSpecies = getModel().getNbSpecies();
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  
  // Implementation of Braginskii heat Flux in 2D
  // this overwrites the values previously stored 
  if (nbSpecies == 2){
    if (_useBackUpValues || _freezeDiffCoeff) {
      for (CFuint i = 0; i < nbSpecies; ++i){
	diffMFMHDData[i] = m_dynViscCoeff[i];
      }
      // AL: here to be commented properly
      for(CFuint i = 0; i < 8; i++){
	diffMFMHDData[nbSpecies + i] = m_thermCondCoeff[i];
      }
    } //End of freezeDiffCoeff
    else {		//this is the default case
      if(dim == DIM_2D){
	const RealVector& mu = getDynViscosityVec(state, gradients);
	if(getModel().isVariableCoeff()) {
	  computeVariableCoeffs(state);
	  diffMFMHDData[0] = _ionVisc;
	  diffMFMHDData[1] = _neutralVisc;
	}
	else { // Constant viscosities
	  for (CFuint i = 0; i < nbSpecies; ++i) { // Reading viscosities
	    diffMFMHDData[i] = mu[i]; //reads the dimensional value imposed in the options
	  }
	}
	
	computeBraginskiiThermConduct(state);             // Compute the Braginskii coeffs.
	const CFreal kappaParallel = _kappaParallel;
	const CFreal kappaPerpendicular = _kappaPerpendicular;
	const CFreal betaWedge = _betaWedge;
	
	//unitary vector in B direction
	const CFreal Bx = state[0];
	const CFreal By = state[1];
	const CFreal B2 = Bx*Bx + By*By;
	const CFreal B = std::sqrt(B2);
	
	CFreal bx = 0;
	CFreal by = 0;
	if (std::abs(B) > 1e-15){//controlling that one doesn't divide by zero
	  bx = Bx/B;
	  by = By/B;
	}
	
	//matrix
	const CFreal bxx = bx*bx;
	const CFreal bxy = bx*by;
	const CFreal byy = by*by;
	const CFuint endVisc = 2;
	diffMFMHDData[endVisc + 0] = kappaParallel*bxx;
	diffMFMHDData[endVisc + 1] = kappaParallel*bxy;
	diffMFMHDData[endVisc + 2] = kappaParallel*byy;
	
	diffMFMHDData[endVisc + 3] = kappaPerpendicular*(1-bxx);
	diffMFMHDData[endVisc + 4] = -kappaPerpendicular*bxy;
	diffMFMHDData[endVisc + 5] = kappaPerpendicular*(1-byy);
	
	diffMFMHDData[endVisc + 6] = betaWedge;
	
	if(getModel().isVariableCoeff()) {  
	  diffMFMHDData[endVisc + 7] = _neutralThermCond; //neutrals' thermal conduction
	}
	else {
	  diffMFMHDData[endVisc + 7] = getThermConductivity(state, mu[1]); //neutrals' thermal conduction
	}
      }
    }
    if (_setBackUpValues) {			//To be implenemted
      for (CFuint i = 0; i < nbSpecies; ++i){
	m_dynViscCoeff[i] = diffMFMHDData[i];
      }
      for(CFuint i = 0; i < 8; i++){
	m_thermCondCoeff[i] = diffMFMHDData[nbSpecies + i];
      }
    }
  } // End of if nbSpecies==2
}

//////////////////////////////////////////////////////////////////////////////

void DiffMFMHDVarSet::computePropertiesInExtendedDomain(const RealVector& state,
							const std::vector<RealVector*>& gradients,
							RealVector& diffMFMHDData)
{
  cf_assert(PhysicalModelStack::getActive()->getDim() == DIM_2D);
  
  const RealVector& mu_0 = getDynViscosityVec(state, gradients);
  const RealVector& mu_f = getModel().getIncreasedDynViscosityDim();
  const CFreal y_boundary = getModel().getTopHeight();
  const CFreal y_0 = getModel().getDampingHeight();
  const CFreal Delta_y = std::abs(y_0 - y_boundary)/5;   // It's taken 5 since tanh(5) is almost 1
  const CFreal y_coord = _faceCoord[YY];
  const CFuint nbSpecies = getModel().getNbSpecies();
  for (CFuint i = 0; i < nbSpecies; ++i) {
    if(y_coord >= y_boundary){ //If the coordinate of the state is in the extended domain, it increases the viscosity
      //reads the dimensional value imposed in the options
      diffMFMHDData[i] =  0.5*mu_f[i]*std::tanh((y_coord - y_0)/Delta_y) + mu_0[i] - 0.5*mu_f[i]*std::tanh((y_boundary - y_0)/Delta_y); 
    }
    else{//If not, it uses the same viscosity
      diffMFMHDData[i] = mu_0[i];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////
      
void DiffMFMHDVarSet::computeStressTensor(const RealVector& state,
					  const vector<RealVector*>& gradients,
					  const CFreal& radius) 
{
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  const CFuint nbSpecies = getModel().getNbSpecies();   
  
  for (CFuint i = 0; i < nbSpecies; ++i) {
    CFreal divTerm = 0.0;
    if (dim == DIM_2D && radius > MathConsts::CFrealEps()) {
      // if the face is a boundary face, the radius could be 0
      // check against eps instead of 0. for safety
      divTerm = _gradState[m_vID[i]]/radius;
    }
    else if (dim == DIM_3D) {
      const RealVector& gradW = *gradients[m_wID[i]];
      divTerm = gradW[ZZ];
    }
    
    const RealVector& gradU = *gradients[m_uID[i]];
    const RealVector& gradV = *gradients[m_vID[i]];
    const CFreal twoThirdDivV = _twoThird*(gradU[XX] + gradV[YY] + divTerm);
    const RealVector& diffMFMHDData = getModel().getPhysicalData();
    const CFreal coeffTauMu = diffMFMHDData[i];    
    CFLog(DEBUG_MAX, "DiffMFMHDVarSet::computeStressTensor() => coeffTauMu = " << coeffTauMu << "\n");
    
    m_tau[i](XX,XX) = coeffTauMu*(2.*gradU[XX] - twoThirdDivV);
    m_tau[i](XX,YY) = m_tau[i](YY,XX) = coeffTauMu*(gradU[YY] + gradV[XX]);
    m_tau[i](YY,YY) = coeffTauMu*(2.*gradV[YY] - twoThirdDivV);
    
    if (dim == DIM_3D) {
      const RealVector& gradW = *gradients[m_wID[i]];
      m_tau[i](XX,ZZ) = m_tau[i](ZZ,XX) = coeffTauMu*(gradU[ZZ] + gradW[XX]);
      m_tau[i](YY,ZZ) = m_tau[i](ZZ,YY) = coeffTauMu*(gradV[ZZ] + gradW[YY]);
      m_tau[i](ZZ,ZZ) = coeffTauMu*(2.*gradW[ZZ] - twoThirdDivV);
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

RealVector& DiffMFMHDVarSet::getHeatFlux(const RealVector& state,
					 const vector<RealVector*>& gradients,
					 const RealVector& normal)
{
  const CFuint nbSpecies = getModel().getNbSpecies(); 
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  const CFreal mu = getDynViscosity(state, gradients);
  m_lambda = getThermConductivityVec(state, mu);
  
  //Model with Braginskii 
  if(getModel().isBraginskii() == true) {
    if (nbSpecies == 2){
      if(dim == DIM_2D){
	const RealVector& diffMFMHDData = getModel().getPhysicalData();
	const RealVector& gradTi = *gradients[m_TID[0]];
	const CFreal Tx = gradTi[XX];
	const CFreal Ty = gradTi[YY];
	
	const CFuint endVisc = 2;
	
	const CFreal TBx = ((diffMFMHDData[endVisc + 0] + diffMFMHDData[endVisc + 3])*Tx + (diffMFMHDData[endVisc + 1] + diffMFMHDData[endVisc + 4])*Ty);
	const CFreal TBy = ((diffMFMHDData[endVisc + 1] + diffMFMHDData[endVisc + 4])*Tx + (diffMFMHDData[endVisc + 2] + diffMFMHDData[endVisc + 5])*Ty);
	
	_heatFlux[0] = -TBx*normal[XX] - TBy*normal[YY];
	
	const RealVector& gradT = *gradients[m_TID[1]];
	_heatFlux[1] = (-m_lambda[1]*(MathFunctions::innerProd(gradT,normal)));
      }
      if(dim == DIM_3D){
        for (CFuint i = 0; i < nbSpecies; ++i) {
          std::cout<<"DiffMFMHDVarSet::getHeatFlux ==> Braginskii 3D not implemented";
        }
      }
    }
  }
  else {    
    for (CFuint i = 0; i < nbSpecies; ++i) {  
      const RealVector& gradT = *gradients[m_TID[i]];
      _heatFlux[i] = (-m_lambda[i]*(MathFunctions::innerProd(gradT,normal)));
    }
  } 
  return _heatFlux;
}

//////////////////////////////////////////////////////////////////////////////

void DiffMFMHDVarSet::computeBraginskiiThermConduct(const RealVector& state)
{
  //Implemented for a 2D case, to be extended to 3D in the future
  //Physical magnitudes in S.I. units
  CFreal Bx = state[0];
  CFreal By = state[1];
  CFreal Bz = state[2]; 
  CFreal Bnorm = std::sqrt(Bx*Bx + By*By + Bz*Bz);	//[T]
  _kappaParallel = 0;
  _kappaPerpendicular = 0;
  _betaWedge = 0; 
  
  const CFuint endEM = 8;
  const CFuint TiID = endEM + 6;
  const CFuint rhoiID = endEM;
  CFreal T_i     = state[TiID];				//[K]
  CFreal rho_i   = state[rhoiID];			//[kg/m3]
  CFreal m_i     = m_eulerModel->getMolecularMass3();	//[kg]
  CFreal m_e     = m_eulerModel->getMolecularMass1();	//[kg]
  CFreal n_i     = rho_i/m_i;				//[m-3]
  CFreal c       = m_eulerModel->getLightSpeed();	//[m/s]
  CFreal eCharge = Framework::PhysicalConsts::ElectronCharge();	//[C]
  CFreal kBoltz = 1.380648813e-16*11604.5052;				//[erg/eV]  
  
  //Conversion to cgs system
  Bnorm *= 1e4;						//[gauss]
  T_i *= 1/11604.5052;					//[erg]
  n_i *= 1e-6;						//[cm-3]
  c *= 100;						//[cm/s]
  eCharge *= c/10;					//[Fr]
  m_i *= 1000;						//[gr]
  m_e *= 1000;						//[gr]
  
  CFreal n_e = n_i; 
  CFreal T_e = T_i;
  
  //General parameters
  const CFreal lambda = 23 - 0.5*std::log(n_i) + 1.5*std::log(T_i); /*23.4 - 1.15*std::log(n_i) + 3.45*std::log(T_i);*/
  const CFreal pi = MathTools::MathConsts::CFrealPi();
  //ion's parameters
  const CFreal mu = 1; //m_i/m_p ion-proton mass ratio. Here ions are considered protons
  const CFreal tau_i = 2.09e7*std::pow(T_i, 3/2)/(lambda*n_i)*mu;// NRL Plasma Formulary /*3*std::sqrt(m_i)*std::pow(T_i, 3/2)/(4*std::sqrt(pi)*lambda*std::pow(eCharge, 4)*n_i);*/
  const CFreal omega_i = (eCharge*Bnorm)/(m_i*c);
  const CFreal chi_i = omega_i*tau_i;
  const CFreal Delta_i = std::pow(chi_i, 4) + 2.7*std::pow(chi_i, 2) + 0.677; 
  //electron's parameters
  const CFreal tau_e = 3.44e5*std::pow(T_i, 3/2)/(lambda*n_i)*mu;// NRL Plasma Formulary 3*std::sqrt(m_e)*std::pow(T_e, 3/2)/(4*std::sqrt(2*pi)*lambda*std::pow(eCharge, 4)*n_e);
  const CFreal omega_e = (eCharge*Bnorm)/(m_e*c);
  const CFreal chi_e = omega_e*tau_e;
  const CFreal Delta_e = std::pow(chi_e, 4) + 14.79*std::pow(chi_e, 2) + 3.7703;  
  
  _kappaParallel = (3.906*tau_i/m_i + 3.1616*tau_e/m_e)*n_i*kBoltz*kBoltz*T_i; 
  _kappaPerpendicular = (((2*std::pow(chi_i, 2)+ 2.645)/Delta_i)*(tau_i/m_i) + ((4.664*std::pow(chi_e, 2)+ 11.92)/Delta_e)*(tau_e/m_e))*n_i*kBoltz*kBoltz*T_i;
  _betaWedge = (chi_e*(3/2*std::pow(chi_e, 2) + 3.053)/Delta_e)*n_e*kBoltz*kBoltz*T_e;
  
  //Conversion into S.I. units
  _kappaParallel *= 1/(1e5);//sure about this
  _kappaPerpendicular *= 1/(1e5);//sure about this
  _betaWedge *= 1/(1e5); //Not sure about this
  
  //   
//   std::cout<<"\n";
//   std::cout<<"Ion's Properties in cgs \n";
//   std::cout<<"lambda  = " << lambda <<"\n";
//   std::cout<<"\t log(n_i)  = " << std::log(n_i) <<"\n";
//   std::cout<<"\t log(T_i)  = " << std::log(T_i) <<"\n";
//   std::cout<<"tau_i   = " << tau_i <<"\n";
//   std::cout<<"omega_i = " << omega_i <<"\n";
//   std::cout<<"chi_i   = " << chi_i <<"\n";
//   std::cout<<"Delta_i = " << Delta_i <<"\n";
//   
//   std::cout<<"\n";
//   std::cout<<"Electron's Properties in cgs \n";
//   std::cout<<"tau_e   = " << tau_e <<"\n";
//   std::cout<<"omega_e = " << omega_e <<"\n";
//   std::cout<<"chi_e   = " << chi_e <<"\n";
//   std::cout<<"Delta_e = " << Delta_e <<"\n";
 
  //std::cout<<"kappaParallel[W/mK] = "<< _kappaParallel <<"\n";
  //std::cout <<"DiffMFMHDVarSet::computeBraginskiiThermConduct => _kappaParallel      =" << _kappaParallel <<"\n";
  //std::cout <<"DiffMFMHDVarSet::computeBraginskiiThermConduct => _kappaPerpendicular =" << _kappaPerpendicular <<"\n";
  //std::cout <<"DiffMFMHDVarSet::computeBraginskiiThermConduct => _betaWedge          =" << _betaWedge <<"\n";
  //This function should be debugged
  //   _kappaParallel = 0.;
  //   _kappaPerpendicular = 0.;
  //   _betaWedge = 0.;
  //   bool stop = true;
  //   cf_assert(stop == false);
}

////////////////////////////////////////////////////////////////////////////// 

void DiffMFMHDVarSet::computeVariableCoeffs(const RealVector& state){

  const CFuint endEM  = 8;
  const CFuint rhoiID = endEM;
  const CFuint rhonID = endEM + 1;
  const CFuint TiID   = endEM + 6;  
  const CFuint TnID   = endEM + 7;
  const CFreal rho_i   = state[rhoiID];                       //[kg/m3]
  const CFreal rho_n   = state[rhonID];                       //[kg/m3]
  const CFreal T_i     = state[TiID];				//[K]
  const CFreal T_n     = state[TnID];                         //[K]
  const CFreal m_i     = m_eulerModel->getMolecularMass3();   //[kg]
  const CFreal m_n     = m_eulerModel->getMolecularMass2();   //[kg]
  const CFreal m_e     = m_eulerModel->getMolecularMass1();   //[kg]
  const CFreal n_i     = rho_i/m_i;                           //[m-3]
  const CFreal n_n     = rho_n/m_n;                           //[m-3]
  const CFreal c       = m_eulerModel->getLightSpeed();       //[m/s]
  const CFreal eCharge = Framework::PhysicalConsts::ElectronCharge(); //[C]
  const CFreal kBoltz  = Framework::PhysicalConsts::Boltzmann(); //[J/K]
  const CFreal epsilon = Framework::PhysicalConsts::VacuumPermittivity();
  
  const CFreal Epsilon_nn = 7.73e-19;
  const CFreal pi         = MathConsts::CFrealPi();
  const CFreal nu_nn      = n_n*Epsilon_nn*std::sqrt(16*kBoltz*T_n/(pi*m_n));
  _neutralVisc = n_n*kBoltz*T_n/nu_nn;  
  _neutralThermCond = 4*n_n*kBoltz*kBoltz*T_n/(m_n*nu_nn);

  const CFreal lambda_C   = 10;
  const CFreal r_Debye    = eCharge*eCharge/(4*pi*epsilon*kBoltz*T_i);
  const CFreal Epsilon_ii = lambda_C*pi*r_Debye*r_Debye;
  const CFreal nu_ii      = 4./3.*n_i*Epsilon_ii*std::sqrt(2.*kBoltz*T_i/(pi*m_i));
  _ionVisc = n_i*kBoltz*T_i/nu_ii; 

  //cout << "DiffMFMHDVarSet::computeVariableCoeffs \n";
  //cout << "_neutralVisc = " << _neutralVisc << "\n";
  //cout << "_neutralThermCond = " << _neutralThermCond << "\n";
  //cout << "_ionVisc = " << _ionVisc << "\n";
}

//////////////////////////////////////////////////////////////////////////////

void DiffMFMHDVarSet::computeHeatFluxScalar(const vector<RealVector*>& gradients,
					    const CFuint i)
{
  const RealVector& diffMFMHDData = getModel().getPhysicalData();
  const CFuint nbSpecies = getModel().getNbSpecies();
  const CFreal lambda = diffMFMHDData[nbSpecies + i]; // single value
  _qFluxVect[i] = -lambda*(*gradients[m_TID[i]]);
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NavierStokes

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
