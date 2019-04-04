#include "MutationppI/MutationLibraryppDebug.hh"
#include "MutationppI/Mutationpp.hh"
#include "Common/CFLog.hh"
#include "Common/Stopwatch.hh"
#include "Environment/ObjectProvider.hh"
#include "Common/StringOps.hh"
#include <fstream>

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::Mutation2OLD;

using Mutation::Thermodynamics::X_TO_Y;
using Mutation::Thermodynamics::Y_TO_X;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace Mutationpp {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<MutationLibraryppDebug,
			    PhysicalPropertyLibrary,
			    MutationppModule,
			    1>
mutationppDebugLibraryProvider("MutationppDebug");

//////////////////////////////////////////////////////////////////////////////

void MutationLibraryppDebug::defineConfigOptions(Config::OptionList& options)
{
  // the following is already defined in Mutation2OLD
  // options.addConfigOption< std::string >("mixtureName","Name of the mixture."); 
  options.addConfigOption< std::string >
    ("StateModelName","Name of the state model (e.g. \"Equil\", \"ChemNonEq1T\", \"ChemNonEq1TTv\").");
  options.addConfigOption< CFreal >("MinRhoi","Minimum partial density."); 
  options.addConfigOption< CFreal >("MinT","Minimum temperature."); 
  options.addConfigOption< CFdouble >("deltaT","Delta temperature.");
  options.addConfigOption< bool >("useLookUpTable","Flag telling if to use the look up tables.");
	options.addConfigOption< bool >("pLogScale","Flag to activate the logarithmic scale for the pressure in the lookup Table.");
  options.addConfigOption< std::vector<std::string> >("lookUpVars","Name of the vars to store in the table.");
  options.addConfigOption< CFdouble >("Tmax","Maximum temperature in the table.");
  options.addConfigOption< CFdouble >("Tmin","Minimum temperature in the table.");
  options.addConfigOption< CFdouble >("Pmax","Maximum pressure in the table.");
  options.addConfigOption< CFdouble >("Pmin","Minimum pressure in the table.");
  options.addConfigOption< CFdouble >("deltaP","Delta pressure.");
}
      
//////////////////////////////////////////////////////////////////////////////

MutationLibraryppDebug::MutationLibraryppDebug(const std::string& name) :
  MutationLibrary2OLD(name),
  m_gasMixture(CFNULL),
  m_gasMixtureEquil(CFNULL),
  m_smType(),
  _nbTvibLocal(0),
  m_y(),
  m_x(),
  m_yn(),
  m_xn(),
  m_molarmassp(),
  m_df(),
  m_rhoivBkp(),
  m_rhoiv(),
  m_ht(),
  m_hr(),
  m_hf(),  
  m_Tstate(),
  _nameToIdxVar(), //@modif_LkT
  _lookUpTables()
{
  addConfigOptionsTo(this);
  
  // the following is already definied in Mutation2OLD
  //_mixtureName = "";
  //setParameter("mixtureName",&_mixtureName);
  
  _stateModelName = "Equil"; // equilibrium by default
  setParameter("StateModelName",&_stateModelName);
  
  _minRhoi = 0.;
  setParameter("MinRhoi",&_minRhoi);
  
  _minT = 0.;
  setParameter("MinT",&_minT);
  
  // change default
  m_shiftHO = true;
  _electrEnergyID = 0;

  _useLookUpTable = false;
  setParameter("useLookUpTable",&_useLookUpTable);

  _lkpVarNames = vector<std::string>();
  setParameter("lookUpVars",&_lkpVarNames);

  _Tmin = 100.0;
  setParameter("Tmin",&_Tmin);

  _Tmax = 2000.0;
  setParameter("Tmax",&_Tmax);

  _deltaT = 10.0;
  setParameter("deltaT",&_deltaT);
  
  _pLogScale = false;
  setParameter("pLogScale",&_pLogScale); 

  _pmin = 10000.0;
  setParameter("Pmin",&_pmin);

  _pmax = 1000000.0;
  setParameter("Pmax",&_pmax);

  _deltaP = 1000.0;
  setParameter("deltaP",&_deltaP);

}

//////////////////////////////////////////////////////////////////////////////

MutationLibraryppDebug::~MutationLibraryppDebug()
{
}

//////////////////////////////////////////////////////////////////////////////

void MutationLibraryppDebug::setup()
{
  CFLog(VERBOSE, "MutationLibraryppDebug::setup() => start\n"); 
  
  MutationLibrary2OLD::setup();
  
  Mutation::MixtureOptions mo(_mixtureName);
  mo.setStateModel(_stateModelName);
  m_gasMixture.reset(new Mutation::Mixture(mo));

  if (_stateModelName != "Equil") {
    Mutation::MixtureOptions moEquil(_mixtureName);
    moEquil.setStateModel("Equil");
    m_gasMixtureEquil = new Mutation::Mixture(moEquil);
  }
  else {
    m_gasMixtureEquil = m_gasMixture.get();
  }
  
  _NS = m_gasMixture->nSpecies();
  m_vecH0.resize(_NS, 0.); 
  m_y.resize(_NS);
  m_x.resize(_NS);
  m_yn.resize(m_gasMixture->nElements());
  m_xn.resize(m_gasMixture->nElements());
  m_charge.resize(_NS);
  m_df.resize(_NS);
  m_molarmassp.resize(_NS); 
  for (CFint i = 0; i < _NS; ++i) {
    m_molarmassp[i] = m_gasMixture->speciesMw(i);
  }
  m_rhoivBkp.resize(_NS);
  m_rhoiv.resize(_NS);
  m_ht.resize(_NS);
  m_hr.resize(_NS);
  m_hf.resize(_NS);  
 
  if (_stateModelName == "Equil") m_smType = LTE;
  if (_stateModelName == "ChemNonEq1T") m_smType = CNEQ;
  if (_stateModelName == "ChemNonEqTTv") m_smType = TCNEQ;
  
  // AL: check this for LTE
  // AL: this needs to be fixed and double checked for the case T-Te
  _nbTvibLocal = m_gasMixture->nEnergyEqns()-1;

  _hasElectrons = m_gasMixture->hasElectrons();  
  
  m_Tstate.resize(_nbTvibLocal+1);
  
  CFLog(VERBOSE, "MutationLibraryppDebug::setup() => _nbTvibLocal = " << _nbTvibLocal << "\n");
  
  // Setup the charge array
  for (CFint i = 0; i < _NS; ++i) {
    m_charge[i] = m_gasMixture->speciesCharge(i);
  }
  
  _Rgas = Mutation::RU;
  cf_assert(_Rgas > 0.);
  
  CFLog(VERBOSE, "MutationLibraryppDebug::setup() => m_shiftHO [" << m_shiftHO << "]\n");
  
  if (m_shiftHO) { 
    /// computation of the H formation at (close to) 0K
    const CFreal P0 = 10.; // pressure can be whatever in this case
    CFreal T0 = 10.; // the equilibrium solver doesn't work with T<10K
    m_gasMixtureEquil->setState(&P0, &T0, 1);
    // composition at (close to) 0K
    RealVector y0(_NS, const_cast<CFreal*>(m_gasMixtureEquil->Y()));
    T0 = 0.00000000001; // set T0 closer to 0K
    m_gasMixture->setState(&y0[0], &T0, 1); // suppose rho=1 (rho value is irrelevant at 0K)
    m_H0 = m_gasMixture->mixtureHMass(T0); 
    cf_always_assert(std::abs(m_H0) > 0.);
    CFLog(VERBOSE, "MutationLibraryppDebug::setup() => y0 = [ " << y0 << "], H0 = [" << m_H0 << "]\n");
    
    // species formation enthalpy at T=0K (this needs to be rechecked for TCNEQ)
    m_gasMixture->speciesHOverRT(&m_vecH0[0], CFNULL, CFNULL, CFNULL, CFNULL);   
    const CFreal RT0 = _Rgas*T0;
    for (CFuint i = 0; i < _NS; ++i) {
      cf_assert(m_molarmassp[i] > 0.);
      const CFreal k = RT0/m_molarmassp[i];
      m_vecH0[i] *= k;
    }
    CFLog(VERBOSE, "MutationLibraryppDebug::setup() => m_vecH0 = [ " << m_vecH0 << "]\n");
    RealVector yH0(m_vecH0*y0);
    CFLog(VERBOSE, "MutationLibraryppDebug::setup() => " << yH0.sum() << " == " << m_H0 << "\n");
  }

  //@modif_LkT
	CFLog(NOTICE, ">> In setup: _lkpVarNames.size() = " <<_lkpVarNames.size() << "\n");

	//////////////////////////////
	// Activate the Look up Table
  if (_lkpVarNames.size() > 0) {
	  CFLog(NOTICE, ">> In setup::Inside the loop IF(_lkpVarNames.size() > 0): " << "\n");
		CFLog(NOTICE, ">> -----> BEFORE the call of setLookUpTables()... " << "\n");
  	setLookUpTables(); //during the call of this function "_useLookUpTable" is still at false
		CFLog(NOTICE, ">> -----> AFTER the call of setLookUpTables() " << "\n");
    _useLookUpTable = true;
		CFLog(NOTICE, ">> -----> Boolean _useLookUpTable = " <<_useLookUpTable<< "\n");
		CFLog(NOTICE, ">> In setup::END loop IF" << "\n");
  }

  CFLog(VERBOSE, "MutationLibraryppDebug::setup() => end\n"); 
}
      
//////////////////////////////////////////////////////////////////////////////

void MutationLibraryppDebug::unsetup()
{
  CFLog(VERBOSE, "MutationLibraryppDebug::unsetup() => start\n"); 
  
  if (_stateModelName != "Equil") {
    delete m_gasMixtureEquil;
  }
  
  if(isSetup()) {
    MutationLibrary2OLD::unsetup();
  }
  
  CFLog(VERBOSE, "MutationLibraryppDebug::unsetup() => start\n"); 
}
      
//////////////////////////////////////////////////////////////////////////////

CFdouble MutationLibraryppDebug::lambdaNEQ(CFdouble& temperature,
				      CFdouble& pressure)
{
  CFreal kNEQ = m_gasMixture->frozenThermalConductivity();
  // RESET_TO_ZERO(kNEQ);
  CFLog(DEBUG_MAX, "Mutation::lambdaNEQ() => k = " << kNEQ << "\n");
  return kNEQ;
}

//////////////////////////////////////////////////////////////////////////////

void MutationLibraryppDebug::lambdaVibNEQ(CFreal& temperature,
				     RealVector& tVec,
				     CFdouble& pressure,
				     CFreal& lambdaTrRo,
				     RealVector& lambdaInt)
{
  RealVector lambdaTRV(_nbTvibLocal+1);
  m_gasMixture->frozenThermalConductivityVector(&lambdaTRV[0]);
  
  lambdaTrRo = lambdaTRV[0];
  for (CFuint i = 0; i < _nbTvibLocal; ++i) {
    lambdaInt[i] = lambdaTRV[i+1];
  }
  
  CFLog(DEBUG_MAX, "Mutation::lambdaVibNEQ() => " << lambdaTrRo 
	<< " " << lambdaInt << "\n");
}
      
//////////////////////////////////////////////////////////////////////////////

CFdouble MutationLibraryppDebug::sigma(CFdouble& temp, //electrical conductivity
				  CFdouble& pressure,
				  CFreal* tVec)
{
  if (temp < 100.) {temp = 100.;}
  //AL: here make sure that if Te is present, electricConductivity() uses that one 
  
  // we are assuming here that setState() has been called before!
  // this way we don't care about given pressure and temperature
  return m_gasMixture->electricConductivity();
}

//////////////////////////////////////////////////////////////////////////////

void MutationLibraryppDebug::gammaAndSoundSpeed(CFdouble& temp,
					   CFdouble& pressure,
					   CFdouble& rho,
					   CFdouble& gamma,
					   CFdouble& soundSpeed)
{
  gamma = m_gasMixture->mixtureEquilibriumGamma();
  soundSpeed = m_gasMixture->equilibriumSoundSpeed();
  
  CFLog(DEBUG_MAX, "Mutation::gammaAndSoundSpeed() => [t,p,rho] = [" 
	<< temp << ", "  << pressure << ", " << rho <<  "] => gamma = " 
	<< gamma << ", soundSpeed = " << soundSpeed << "\n");
}
      
//////////////////////////////////////////////////////////////////////////////

void MutationLibraryppDebug::frozenGammaAndSoundSpeed(CFdouble& temp,
						 CFdouble& pressure,
						 CFdouble& rho,
						 CFdouble& gamma,
						 CFdouble& soundSpeed,
						 RealVector* tVec)
{
  gamma = m_gasMixture->mixtureFrozenGamma();
  // soundSpeed = m_gasMixture->frozenSoundSpeed();
  soundSpeed = std::sqrt(gamma*pressure/rho); 
  
  CFLog(DEBUG_MAX, "Mutation::gammaAndSoundSpeed() => [t,p,rho] = [" 
	<< temp << ", "  << pressure << ", " << rho <<  "] => gamma = " 
	<< gamma << ", soundSpeed = " << soundSpeed << "\n");
}
      
//////////////////////////////////////////////////////////////////////////////
      
CFdouble MutationLibraryppDebug::soundSpeed(CFdouble& temp, CFdouble& pressure)
{
	if(!_useLookUpTable) {
  	m_gasMixture->setState(&pressure, &temp, 1);
  	return m_gasMixture->equilibriumSoundSpeed();
	}
	else {
		return _lookUpTables.get(temp, pressure, _nameToIdxVar.find("a"));
		}
}

//////////////////////////////////////////////////////////////////////////////
      
void MutationLibraryppDebug::setComposition(CFdouble& temp,
				       CFdouble& pressure,
				       RealVector* x)
{
	if(!_useLookUpTable){
		
		if (temp < 100.) {temp = 100.;}
		// ___________________________________________________________ 
		// Test for knowing the number of time this function is called
		//std::cout<<MutationLibrarypp::k<<std::endl;
		//MutationLibrarypp::k = MutationLibrarypp::k+1;

		 m_gasMixtureEquil->setState(&pressure, &temp, 1);

		const double* xm = m_gasMixtureEquil->X();
		
		if (x != CFNULL) {
		  for(CFint i = 0; i < _NS; ++i) {
		    (*x)[i] = xm[i];
		  }
		}
		
		m_gasMixtureEquil->convert<X_TO_Y>(xm, &m_y[0]);
	}
}
      
//////////////////////////////////////////////////////////////////////////////

void MutationLibraryppDebug::setDensityEnthalpyEnergy(CFdouble& temp,
						 CFdouble& pressure,
						 RealVector& dhe)
{
	//setComposition must be called before
  if(!_useLookUpTable){
  	CFLog(DEBUG_MAX, "Mutation::setDensityEnthalpyEnergy() => P = " 
		<< pressure << ", T = " << temp << "\n");
  
  	dhe[0] = m_gasMixture->density();
  	dhe[1] = m_gasMixture->mixtureHMass() - m_H0;
  	dhe[2] = dhe[1]-pressure/dhe[0];
  
  	CFLog(DEBUG_MAX, "Mutation::setDensityEnthalpyEnergy() => " << dhe << ", " <<  m_y << "\n");
	}
	else{
    dhe[0] = _lookUpTables.get(temp, pressure, _nameToIdxVar.find("d"));
    dhe[1] = _lookUpTables.get(temp, pressure, _nameToIdxVar.find("h"));
    dhe[2] = _lookUpTables.get(temp, pressure, _nameToIdxVar.find("e"));
  }
}
      
//////////////////////////////////////////////////////////////////////////////

void MutationLibraryppDebug::setDensityEnthalpyEnergy(CFdouble& temp,
						 RealVector& tVec,
						 CFdouble& pressure,
						 RealVector& dhe,
						 bool storeExtraData)
{
  CFLog(DEBUG_MAX, "Mutation::setDensityEnthalpyEnergy() => P = " 
	<< pressure << ", T = " << temp << ", Tv " << tVec << " \n");
  
  dhe[0] = m_gasMixture->density();
  dhe[1] = m_gasMixture->mixtureHMass() - m_H0;
  dhe[2] = dhe[1]-pressure/dhe[0];
 
  CFLog(DEBUG_MAX, "Mutation::setDensityEnthalpyEnergy() => " << dhe << ", " <<  m_y << "\n");
}
      
//////////////////////////////////////////////////////////////////////////////

CFdouble MutationLibraryppDebug::density(CFdouble& temp,
				    CFdouble& pressure,
				    CFreal* tVec)
{
  if (m_smType == LTE && _useLookUpTable) { 
    return _lookUpTables.get(temp, pressure, _nameToIdxVar.find("d"));
  }
  m_gasMixture->setState(&pressure, &temp, 1);
  return m_gasMixture->density();
}

//////////////////////////////////////////////////////////////////////////////

CFdouble MutationLibraryppDebug::pressure(CFdouble& rho,
				     CFdouble& temp,
				     CFreal* tVec)
{
  CFLog(DEBUG_MAX, "Mutation::pressure() => rho = " << rho << ", T = " << temp 
	<< ", y = " << m_y << "\n");
  
  // const CFreal p = m_gasMixture->pressure(temp, rho, &m_y[0]);
  const CFreal p = m_gasMixture->P();
  if (p <= 0.) {
    CFLog(DEBUG_MAX, "Mutation::pressure() => p = " << p << " with rho = " << rho 
	  << ", T = " << temp << ", y = " << m_y << "\n");
  }
  cf_assert(p>0.);
  
  //CFLog(DEBUG_MAX, "Mutation::pressure() => " << p << "\n");
  
  return p;
}

//////////////////////////////////////////////////////////////////////////////
      
CFdouble MutationLibraryppDebug::electronPressure(CFreal rhoE,
					     CFreal tempE)
{
  return rhoE*tempE*_Rgas/m_molarmassp[0];
}
      
//////////////////////////////////////////////////////////////////////////////

CFdouble MutationLibraryppDebug::energy(CFdouble& temp,
				   CFdouble& pressure)
  
{
	if(!_useLookUpTable){
  	m_gasMixture->setState(&pressure, &temp, 1);
  	return m_gasMixture->mixtureEnergyMass()- m_H0;
	}
  else {
		return _lookUpTables.get(temp, pressure, _nameToIdxVar.find("e"));
	}
}
      
//////////////////////////////////////////////////////////////////////////////

CFdouble MutationLibraryppDebug::enthalpy(CFdouble& temp,
				     CFdouble& pressure)
{
	if (!_useLookUpTable) {
  	m_gasMixture->setState(&pressure, &temp, 1);
  	return m_gasMixture->mixtureHMass() - m_H0;
	}
  else{
		 return _lookUpTables.get(temp, pressure, _nameToIdxVar.find("h"));
	}

}
      
//////////////////////////////////////////////////////////////////////////////

 void MutationLibraryppDebug::setElemFractions(const RealVector& yn)
 {
  // ICP
  //
   for (CFint ic = 0; ic < _NC; ++ic) {
     m_yn[ic] = yn[ic];

     if (!(m_yn[ic] >= 0.0 && m_yn[ic] <= 1.0)) {
       // cout << "Yn[ic] = " << Yn[ic] << endl;
       // abort();
     }

     assert(m_yn[ic] >= 0.0);
     assert(m_yn[ic] <= 1.0);
   }
  
   //m_gasMixture->convert<YE_TO_XE>(&m_yn[0], &m_xn[0]);
   throw NotImplementedException(FromHere(),"MutationLibraryppDebug::setElemFractions()");
 }

//////////////////////////////////////////////////////////////////////////////

 void MutationLibraryppDebug::setElementXFromSpeciesY(const RealVector& ys)
 {
   for (CFint ic = 0; ic < _NC; ++ic) {
     m_y[ic] = ys[ic];
   }
   //m_gasMixture->convert<Y_TO_XE>(&m_y[0], &m_xn[0]);
   throw NotImplementedException(FromHere(),"MutationLibraryppDebug::setElementXFromSpeciesY()");
 }

// //////////////////////////////////////////////////////////////////////////////

 void MutationLibraryppDebug::setElectronFraction(RealVector& ys)
 {
   // charge neutrality: xEl = sum(xIon)
   CFdouble yEl = 0.0;
   for (CFint is = 0; is < _NS; ++is) {
     if (m_charge[is] > 0) {
       yEl += ys[is] / m_molarmassp[is];
     }
   }
   
   yEl *= m_molarmassp[0]; // 1st species: electron
   ys[0] = std::max(0., yEl); // overwrite electron mass fraction (keep them positive)
   
   CFLog(DEBUG_MAX, "Mutation::setElectronFraction() => " << ys[0] << "\n");
 }

//////////////////////////////////////////////////////////////////////////////

void MutationLibraryppDebug::setSpeciesFractions(const RealVector& ys)
{
  MutationLibrary2OLD::setSpeciesFractions(ys);

  if (presenceElectron()) {
    setElectronFraction(const_cast<RealVector&>(ys));
  }
  
  for (CFint is = 0; is < _NS; ++is) {
    m_y[is] = ys[is];
    
    if (m_y[is] < 0.0) m_y[is] = 0.0;
    cf_assert(m_y[is] < 1.1);
  }
    
  CFLog(DEBUG_MAX, "Mutation::setSpeciesFractions() => " << m_y << "\n");
}
      
//////////////////////////////////////////////////////////////////////////////

void MutationLibraryppDebug::getSpeciesMolarFractions
(const RealVector& ys, RealVector& xs)
{
  RealVector& yss = const_cast<RealVector&>(ys);
  m_gasMixture->convert<Y_TO_X>(&yss[0], &xs[0]);
  CFLog(DEBUG_MAX, "Mutation::getSpeciesMolarFractions() => " << xs << "\n");
}
      
//////////////////////////////////////////////////////////////////////////////

 void MutationLibraryppDebug::getSpeciesMassFractions
 (const RealVector& xs, RealVector& ys)
 {
   RealVector& xss = const_cast<RealVector&>(xs);
   m_gasMixture->convert<X_TO_Y>(&xss[0], &ys[0]);
   CFLog(DEBUG_MAX, "Mutation::getSpeciesMolarFractions() => " << ys << "\n");
 }

//////////////////////////////////////////////////////////////////////////////

 void MutationLibraryppDebug::getSpeciesMassFractions(RealVector& ys)
 {
   for (CFint is = 0; is < _NS; ++is) {
     ys[is] = m_y[is];
   }
   CFLog(DEBUG_MAX, "Mutation::getSpeciesMolarFractions() => " << ys << "\n");
 }

//////////////////////////////////////////////////////////////////////////////

void MutationLibraryppDebug::getTransportCoefs(CFdouble& temp,
					  CFdouble& pressure,
					  CFdouble& lambda,
					  CFdouble& lambdacor,
					  RealVector& lambdael,
					  RealMatrix& eldifcoef,
					  RealVector& eltdifcoef)
{
 throw NotImplementedException(FromHere(),"MutationLibraryppDebug::getTransportCoefs()");
}

//////////////////////////////////////////////////////////////////////////////

void MutationLibraryppDebug::getMassProductionTerm(CFdouble& temperature,
					      RealVector& tVec,
					      CFdouble& pressure,
					      CFdouble& rho,
					      const RealVector& ys,
					      bool flagJac,
					      RealVector& omega,
					      RealMatrix& jacobian)
{
  // CFLog(DEBUG_MAX, "MutationLibraryppDebug::getMassProductionTerm() => P = " 
  // << pressure << ", T = " << temperature << "\n");
  
  // we assume setState() already called before
  if (!_freezeChemistry) {
    m_gasMixture->netProductionRates(&omega[0]);
  } 
  else {
    omega = 0.;
  }
    
  CFLog(DEBUG_MAX, "Mutation::getMassProductionTerm() => omega = " << omega << "\n\n");
  //EXIT_AT(1000);
  // TODO: Need to figure out how to do the Jacobian later
}
      
//////////////////////////////////////////////////////////////////////////////

void MutationLibraryppDebug::getSource(CFdouble& temperature,
				 RealVector& tVec,
				 CFdouble& pressure,
				 CFdouble& rho,
				 const RealVector& ys,
				 bool flagJac,
				 RealVector& omega,
				 RealVector& omegav,
				 CFdouble& omegaRad,
				 RealMatrix& jacobian)   
  
{
  // we assume setState() already called before
  m_gasMixture->netProductionRates(&omega[0]);
  
  m_gasMixture->energyTransferSource(&omegav[0]);
}
      
//////////////////////////////////////////////////////////////////////////////

void MutationLibraryppDebug::getRhoUdiff(CFdouble& temperature,
				    CFdouble& pressure,
				    RealVector& normConcGradients,
				    RealVector& normTempGradients,
				    CFreal* tVec,
				    RealVector& rhoUdiff,
				    bool fast)
{  
  // Set driving forces as gradients of molar fractions
  CFreal MMass = m_gasMixture->mixtureMw();
  CFreal normMMassGradient = 0.0;
  for (CFint is = 0; is < _NS; ++is) {
    normMMassGradient += normConcGradients[is] / m_molarmassp[is];
  }
  normMMassGradient *= -MMass*MMass;
  
  for (CFint is = 0; is < _NS; ++is) {
    m_df[is] = (MMass*normConcGradients[is] + m_y[is]*normMMassGradient) / m_molarmassp[is];
  }
  
  CFreal E = 0.0;
  m_gasMixture->stefanMaxwell(&m_df[0], &rhoUdiff[0], E);
  
  // CFLog(DEBUG_MAX, "Mutation::rhoUdiff() =>  rhoUdiff = " << rhoUdiff << "\n");
  
  const CFreal density = m_gasMixture->density();
  
  // CFLog(DEBUG_MAX, "Mutation::rhoUdiff() =>  rho = " << density << "\n");
  
  for (CFint is = 0; is < _NS; ++is) {
    rhoUdiff[is] *= m_y[is]*density;
  }
  
  // RESET_TO_ZERO(rhoUdiff); 
  CFLog(DEBUG_MAX, "Mutation::rhoUdiff() => rhoUdiff = " << rhoUdiff << "\n");
  // EXIT_AT(1000);
}
      
//////////////////////////////////////////////////////////////////////////////

void MutationLibraryppDebug::getDij_fick(RealVector& dx,
				    CFdouble& pressure,
				    CFdouble& temperature,
				    CFreal& Diff_coeff,
				    RealMatrix& Dij,
				    RealVector& rhoUdiff)
{
  throw NotImplementedException(FromHere(),"MutationLibraryppDebug::getDij_fick()");
}
      
//////////////////////////////////////////////////////////////////////////////

void MutationLibraryppDebug::getGammaN(CFreal& m_GN)
{
  throw NotImplementedException(FromHere(),"MutationLibraryppDebug::getGammaN()");
}
//////////////////////////////////////////////////////////////////////////////

void MutationLibraryppDebug::getGammaO(CFreal& m_GO)
{
  throw NotImplementedException(FromHere(),"MutationLibraryppDebug::getGammaO()");
}
				  
//////////////////////////////////////////////////////////////////////////////

 void MutationLibraryppDebug::setSpeciesMolarFractions(const RealVector& xs)
 {
   throw NotImplementedException(FromHere(),"MutationLibraryppDebug::setSpeciesMolarFractions()");
 }
      
//////////////////////////////////////////////////////////////////////////////

void MutationLibraryppDebug::getSpeciesTotEnthalpies(CFdouble& temp,
						RealVector& tVec,
						CFdouble& pressure,
						RealVector& hsTot,
						RealVector* hsVib,
						RealVector* hsEl)
{
  // CFLog(DEBUG_MAX, "Mutation::getSpeciesTotEnthalpies()\n");
  
  CFreal* hv = (hsVib != CFNULL) ? &(*hsVib)[0] : CFNULL; 
  CFreal* he = (hsEl  != CFNULL) ? &(*hsEl)[0] : CFNULL;
  CFreal* ht = (hsVib != CFNULL) ? &m_ht[0] : CFNULL;
  CFreal* hr = (hsVib != CFNULL) ? &m_hr[0] : CFNULL;
  CFreal* hf = (hsVib != CFNULL) ? &m_hf[0] : CFNULL;

  m_gasMixture->speciesHOverRT(&hsTot[0], ht, hr, hv, he, hf); 
  
  const CFreal RT = _Rgas*temp;
  for (CFuint i = 0; i < _NS; ++i) {
    const CFreal k0 = RT/m_molarmassp[i];
    hsTot[i] *= k0;
    if (hsVib != CFNULL) {(*hsVib)[i] *= k0;} 
  //  if (presenceElectron()) {
     if (hsEl != CFNULL) {(*hsEl)[i]  *= k0;}
  // }
  }
  
  hsTot -= m_vecH0; // shift on the formation enthalpy (considering H(T=0K)=0)
  CFLog(DEBUG_MAX, "Mutation::getSpeciesTotEnthalpies() => m_vecH0 = " << m_vecH0 << "\n");
  CFLog(DEBUG_MAX, "Mutation::getSpeciesTotEnthalpies() => hsTot = " << hsTot << "\n");
  
  //EXIT_AT(5);
}
      
//////////////////////////////////////////////////////////////////////////////

//@modif_LkT
void MutationLibraryppDebug::setLookUpTables()
{
  vector<ComputeQuantity> varComputeVec;
  varComputeVec.reserve(_lkpVarNames.size());

  // store the pointers to member functions for
  // computing physical quantities
  for (CFuint i = 0; i < _lkpVarNames.size(); ++i) {

    if (_lkpVarNames[i] == "e") {
			CFLog(NOTICE, ">> In SetLookUpTables():_lkpVarNames["<<i<<"] == 'e' \n");
      varComputeVec.push_back(&MutationLibraryppDebug::energy);
    }
    else if (_lkpVarNames[i] == "h") {
			CFLog(NOTICE, ">> In SetLookUpTables():_lkpVarNames["<<i<<"] == 'h' \n");
      varComputeVec.push_back(&MutationLibraryppDebug::enthalpy);
    }
    else if (_lkpVarNames[i] == "a") {
			CFLog(NOTICE, ">> In SetLookUpTables():_lkpVarNames["<<i<<"] == 'a'  \n");
      varComputeVec.push_back(&MutationLibraryppDebug::soundSpeed);
    }
    else if (_lkpVarNames[i] == "d") {
			CFLog(NOTICE, ">> In SetLookUpTables():_lkpVarNames["<<i<<"] == 'd'  \n");
      varComputeVec.push_back(&MutationLibraryppDebug::density);
    }
    else {
      throw Common::NoSuchValueException (FromHere(), "Variable name not found");
    }
  }

  // set the look up tables
  CFLog(NOTICE, ">> In SetLookUpTables(): -------> BEFORE the call of setTables() \n");
  setTables(varComputeVec);
	CFLog(NOTICE, ">> In SetLookUpTables(): -------> AFTER the call of setTables() \n");
}
////////////////////////////////////////////////////////////////////////////////

void MutationLibraryppDebug::setTables(vector<ComputeQuantity>& varComputeVec)
{
  Common::Stopwatch<Common::WallTime> stp;
  stp.start(); 
  const CFuint nbLookUpVars = _lkpVarNames.size();
  
  //_____________________________________________________________________________
  // create mapping variable names with corresponding idxs
	//https://github.com/andrealani/COOLFluiD/blob/master/src/Common/CFMap.hh
	CFLog(NOTICE, ">> In setTables(): Creating and printing the mapping CFmap... \n");
  _nameToIdxVar.reserve(nbLookUpVars);
  for (CFuint iVar = 0; iVar < nbLookUpVars; ++iVar) {
    _nameToIdxVar.insert(_lkpVarNames[iVar], iVar);
	  
  }	
	_nameToIdxVar.print();
  _nameToIdxVar.sortKeys();
	//a => 0
	//e => 1
	//h => 2
	//d => 3
	CFLog(NOTICE, ">> In setTables(): After sorting the mapping CFmap... \n");
	_nameToIdxVar.print();
	//a => 0
	//d => 3
	//e => 1
	//h => 2


  //_____________________________________________________________________________
  // compute the arrays containing all the keys for T and P
	// Based on the value set in CFcase file

	//_________________________________
	// Linear scale for the temperature
	//---------------------------------
	const CFuint nbT = static_cast<CFuint> ((_Tmax - _Tmin)/_deltaT) +1 ;
		vector<CFdouble> vT(nbT+2);

	CFLog(NOTICE, ">> 	           __________________________________\n");
	CFLog(NOTICE, ">> In setTables(): Linear scale for the temperature\n");
	CFLog(NOTICE, ">> In setTables(): Tmin           = "<<_Tmin<<"\n");
	CFLog(NOTICE, ">> In setTables(): Tmax           = "<<_Tmax<<"\n");
	CFLog(NOTICE, ">> In setTables(): deltaT         = "<<_deltaT<<"\n");
	CFLog(NOTICE, ">> In setTables(): nbT = "<<nbT<<"\n");
  for(CFuint i = 0; i <= nbT ; ++i) {
    vT[i] = _Tmin + i*_deltaT;
  }
		CFLog(NOTICE, ">> In setTables(): vT[0] = "<<vT[0]<<",  vT[1] = "<<vT[1]<<"... "<<"vT["<<nbT-1<<"] = "<<vT[nbT-1]<<"\n");

	//_________________________________
	// Linear scale for the pressure
	//---------------------------------

	const CFuint nbP_lin = static_cast<CFuint> ((_pmax - _pmin)/_deltaP)+1; 
	vector<CFdouble> vP_lin(nbP_lin+1);
	for(CFuint i = 0; i <= nbP_lin  ; ++i) {
    vP_lin[i] = _pmin + i*_deltaP;
  }

	//__________________________________
	//Logarithmic scale for the pressure
	//----------------------------------
  //Assumption : Pressure >= 1
	// step 0: Writing P in a Scientific expression form -> a.10^b 
	int b_min = log10(_pmin);
	CFdouble a_min = _pmin*pow(10,-b_min);
	int b_max = log10(_pmax);
	CFdouble a_max = _pmax*pow(10,-b_max);
	int pixels = (_pmax-_pmin)/_deltaP;
	CFdouble NbPowerOfTen = 0;

	/* 
	 * STEP 1: Number of 10th-power between _pmin and _pmax
	 */
	if(a_max == 1)  NbPowerOfTen = b_max - b_min + 1;
	else            NbPowerOfTen = b_max - b_min + 2;
	/* 
	 * STEP 2: Number of logarithmic units (bloc of size 10 => StartBloc[i+1]/StartBloc[i] = 10)
	 */
 	vector<CFdouble> StartBloc;                 
	for(CFuint i=0; i<NbPowerOfTen - 1;++i){
		StartBloc.push_back(pow(10,b_min+i));
	}	
	CFuint nbBlocs = NbPowerOfTen -1;
	/* 
	 * STEP 3: EquiDistribution of the pixels in each blocs 
	 */ 
	CFuint NbPixPerBloc = CFuint (pixels/nbBlocs); 
	/* 
	 * STEP 4: Set the boundaries of each blocs in the vector	
	 */
		CFuint nbP_log = nbBlocs*NbPixPerBloc + 1;
		vector<CFdouble>  vPbis(nbP_log);
		CFuint PixelPmax;
		for(CFuint i=0 ; i<nbBlocs ; ++i){
			vPbis[i*NbPixPerBloc] = StartBloc[i];
		}
		if (a_max ==1){
				 PixelPmax = nbBlocs*NbPixPerBloc;
				 vPbis[PixelPmax] = _pmax;
		}
		else vPbis[nbBlocs*NbPixPerBloc] = pow(10,b_max+1);
	/* 
	 * STEP 5: set pmin, pmax and the values 2, 3, 4, ...,9 between each bloc of ratio 10 in the vector
	 *					_pmin - 100 - 200 - ... - 1000 - 1100 - ... - 10 000 - 10 100 - ... - 100 000 - 100 100 - ... - _pmax
	 */
	vector<CFdouble> Res(9);
	CFuint PixelPmin = int(log10(a_min)*NbPixPerBloc);
	PixelPmax = int(log10(a_max)*NbPixPerBloc + (nbBlocs-1)*NbPixPerBloc);
	vPbis[PixelPmin] = _pmin;
	if(a_max>1) vPbis[PixelPmax] = _pmax;
	for(CFuint i=0 ; i<nbBlocs ; ++i){
			for(CFuint j=1 ; j<10 ; ++j){
					vPbis[int(NbPixPerBloc*log10(j)) + i*NbPixPerBloc] = j*StartBloc[i];
			}   
	}	
	/* 
	 * STEP 6: set the other values between the ones set previously (step 5)
	 */ 
	for(CFdouble pix=0 ; pix<NbPixPerBloc ; ++pix){	 	
	  for(CFuint k=0 ; k<nbBlocs ; ++k){
	    if(vPbis[pix+k*NbPixPerBloc] ==0){
	      vPbis[pix+k*NbPixPerBloc] = pow(10,double(pix/NbPixPerBloc))	*StartBloc[k];	
	    }
	  }		
	}
	/* 
	 * STEP 7: Resize the vector vPbis by removing the values in excess () in the boundaries
	 *        (10^bmin ...) - pmin - ........- pmax - (... 10^bmax+1)
	 */
	CFuint k = 0;
	nbP_log = PixelPmax-PixelPmin+1;
	vector<CFdouble> vP_log(nbP_log+1);
	 for(CFuint pix=PixelPmin ; pix<=PixelPmax ; ++pix){
	   vP_log[k]=vPbis[pix];
	   k=k+1;
	 }
	///... END STEPS

	// _________________________________________________________________
	// CHOICE (linear or logarithmic) OF THE PRESSURE SCALE BY THE USER:
	CFuint nbP = (1-_pLogScale)*nbP_lin + _pLogScale*nbP_log;
	vector<CFdouble>  vP( nbP +1 );
	if(_pLogScale){
		CFLog(NOTICE, ">> 	           __________________________________\n");
		CFLog(NOTICE, ">> In setTables(): Logarithmic scale for the pressure\n");
		CFLog(NOTICE, ">> In setTables(): pmin   = "<<a_min<<"x 10^"<<b_min<<"\n");
		CFLog(NOTICE, ">> In setTables(): pmax   = "<<a_max<<"x 10^"<<b_max<<"\n");
		CFLog(NOTICE, ">> In setTables(): deltaP = "<<_deltaP<<"\n");
		CFLog(NOTICE, ">> In setTables(): nbP_log before truncation = "<< pixels<<"\n");
		CFLog(NOTICE, ">> In setTables(): NbPowerOfTen = "<<NbPowerOfTen<<"\n");
		CFLog(NOTICE, ">> In setTables(): nbBlocs = "<<nbBlocs<<"\n");
		CFLog(NOTICE, ">> In setTables(): NbPixPerBloc = "<<NbPixPerBloc<<"\n");
	  CFLog(NOTICE, ">> In setTables(): nbP_log after removing the excess values = "<< nbP_log<<"\n");
		vP = vP_log;
		CFLog(NOTICE, ">> In setTables(): vP[0] = "<<vP[0]<<",  vP[1] = "<<vP[1]<<",  vP[2] = "<<vP[2]<<"... "<<"vP["<<nbP-1<<"] = "<<vP[nbP-1]<<"\n");

	}
	else{
		CFLog(NOTICE, ">> In setTables(): nbP_lin = "<<nbP_lin<<"\n");
		CFLog(NOTICE, ">> 	           __________________________________\n");
		CFLog(NOTICE, ">> In setTables(): Linear scale for the pressure\n");
		CFLog(NOTICE, ">> In setTables(): Pmin           = "<<_pmin<<"\n");
		CFLog(NOTICE, ">> In setTables(): Pmax           = "<<_pmax<<"\n");
		CFLog(NOTICE, ">> In setTables(): deltaP         = "<<_deltaP<<"\n");
		CFLog(NOTICE, ">> In setTables(): nb subinterval = "<<nbP_lin<<"\n");
		vP = vP_lin;
		CFLog(NOTICE, ">> In setTables(): vP[0] = "<<vP[0]<<",  vP[1] = "<<vP[1]<<",  vP[2] = "<<vP[2]<<"... "<<"vP["<<nbP-1<<"] = "<<vP[nbP-1]<<"\n");
	}

	//_____________________________________________________________________________
  // INITIALIZE THE TABLE AND INSERT THE VALUES
	//https://github.com/andrealani/COOLFluiD/blob/master/src/Common/LookupTable2D.ci
  _lookUpTables.initialize(vT, vP, nbLookUpVars);
		CFLog(NOTICE, ">> 	           ____________________________________________\n");
	CFLog(NOTICE, ">> In setTables(): inserting the values in the _lookUpTables (Takes time...) \n");
	CFdouble result;
	for(CFuint i = 0; i < nbT; ++i) {
    for(CFuint j = 0; j < nbP; ++j) {

      // set the composition at first (Here _useLookUpTable is still at false)
      setComposition(vT[i],vP[j], CFNULL); 
      for (CFuint iVar  = 0; iVar < nbLookUpVars; ++iVar) {
				
       result = (this->*varComputeVec[iVar])(vT[i],vP[j]); 
			 // Meaning of this command at each iter:
			 // result := varComputeVec[0](vT[i],vP[j]) := SoundSpeed(vT[i],vP[j])
			 // result := varComputeVec[1](vT[i],vP[j]) := Energy(vT[i],vP[j])
			 // result := varComputeVec[2](vT[i],vP[j]) := Enthalpy(vT[i],vP[j])
			 // result := varComputeVec[3](vT[i],vP[j]) := Density(vT[i],vP[j])


       _lookUpTables.insert(vT[i],vP[j], _nameToIdxVar.find(_lkpVarNames[iVar]),result);
      }
    }
  }
	//_____________________________________________________________________________
  // WRITE IN A FILE THE VALUES INSERTED IN THE LKT
	/*
  ofstream fichier("bliblablou.dat", ios::out | ios::trunc);  
	if(fichier){
		CFLog(NOTICE, "\n>> In setTables(): printing values computed in a .dat file:\n\n");

		//for(CFuint i = 0; i < nbT; ++i) {
			for(CFuint j = 0; j < nbP; ++j) {
				fichier << vP[j]<< "  "<<_lookUpTables.get(vP[j],vT[0], _nameToIdxVar.find("h"))<<endl;
			}
		//}

		fichier.close();
	}
	else cerr << "CAN NOT OPEN THE  !" << endl;
	*/
}
////////////////////////////////////////////////////////////////////////////////

void MutationLibraryppDebug::getSourceTermVT(CFdouble& temperature,
					RealVector& tVec,
					CFdouble& pressure,
					CFdouble& rho,
					RealVector& omegav,
					CFdouble& omegaRad)
{
  // AL: I guess that even the case T-Te should be treated with the same function
  m_gasMixture->energyTransferSource(&omegav[0]);
  
  CFLog(DEBUG_MAX, "Mutation::getSourceTermVT() => omegav = " << omegav << "\n");
}
      
//////////////////////////////////////////////////////////////////////////////

void MutationLibraryppDebug::getMolarMasses(RealVector& mm)
{
  assert(mm.size() == static_cast<CFuint>(_NS));

  // check the units
  for (CFint i = 0; i < _NS; ++i) {
    mm[i] = m_molarmassp[i];
  }
  
  CFLog(DEBUG_MAX, "Mutation::getMolarMasses() => " << mm << "\n");
}

//////////////////////////////////////////////////////////////////////////////
      
void MutationLibraryppDebug::transportCoeffNEQ(CFreal& temperature, 
					  CFdouble& pressure,
					  CFreal* tVec, 
					  RealVector& normConcGradients,
					  RealVector& normTempGradients,
					  CFreal& eta,
					  CFreal& lambdaTrRo, 
					  RealVector& lambdaInt,
					  RealVector& rhoUdiff)
{
   throw NotImplementedException(FromHere(),"MutationLibraryppDebug::transportCoeffNEQ()");
}
        
//////////////////////////////////////////////////////////////////////////////
   
void MutationLibraryppDebug::getSourceEE(CFdouble& temperature,
				    RealVector& tVec,
				    CFdouble& pressure,
				    CFdouble& rho,
				    const RealVector& ys,
				    bool flagJac,
				    CFdouble& omegaEE)
{
  throw NotImplementedException(FromHere(),"MutationLibraryppDebug::getSourceEE()");
}
      
//////////////////////////////////////////////////////////////////////////////
      
} // namespace Mutationpp

} // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

