#include "MutationppI/MutationLibrarypp.hh"
#include "MutationppI/Mutationpp.hh"
#include "Common/CFLog.hh"
#include "Environment/ObjectProvider.hh"
#include "Common/StringOps.hh"
#include <fstream>

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;

using Mutation::Thermodynamics::X_TO_Y;
using Mutation::Thermodynamics::Y_TO_X;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace Mutationpp {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<MutationLibrarypp,
			    PhysicalPropertyLibrary,
			    MutationppModule,
			    1>
mutationppLibraryProvider("Mutationpp");

//////////////////////////////////////////////////////////////////////////////

void MutationLibrarypp::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< std::string >("mixtureName","Name of the mixture."); 
  options.addConfigOption< std::string >
    ("StateModelName","Name of the state model (e.g. \"Equil\", \"ChemNonEq1T\", \"ChemNonEq1TTv\").");
  options.addConfigOption< bool >("ShiftH0","Shift the formation enthalpy to have H(T=0K)=0."); 
}
      
//////////////////////////////////////////////////////////////////////////////

MutationLibrarypp::MutationLibrarypp(const std::string& name)
  : Framework::PhysicalChemicalLibrary(name),
    m_gasMixture(CFNULL),
    m_smType(),
    m_H0(0.),
    m_vecH0(),
    m_y(),
    m_x(),
    m_yn(),
    m_xn(),
    m_molarmassp(),
    m_df(),
    m_rhoivBkp(),
    m_rhoiv()
{
  addConfigOptionsTo(this);
  
  _mixtureName = "";
  setParameter("mixtureName",&_mixtureName);
  
  _stateModelName = "Equil"; // equilibrium by default
  setParameter("StateModelName",&_stateModelName);
  
  m_shiftHO = true;
  setParameter("ShiftH0",&m_shiftHO);
}

//////////////////////////////////////////////////////////////////////////////

MutationLibrarypp::~MutationLibrarypp()
{
}

//////////////////////////////////////////////////////////////////////////////

void MutationLibrarypp::configure ( Config::ConfigArgs& args )
{
  Framework::PhysicalChemicalLibrary::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void MutationLibrarypp::setup()
{    
  CFLog(VERBOSE, "MutationLibrarypp::setup() => start\n"); 
  
  Framework::PhysicalChemicalLibrary::setup();
  
  Mutation::MixtureOptions mo(_mixtureName);
  mo.setStateModel(_stateModelName);
  m_gasMixture.reset(new Mutation::Mixture(mo));
  
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
  
  if (_stateModelName == "Equil") m_smType = LTE;
  if (_stateModelName == "ChemNonEq1T") m_smType = CNEQ;
  if (_stateModelName == "ChemNonEq1TTv") m_smType = TCNEQ;
  
  // AL: check this for LTE
  _nbTvib = m_gasMixture->nEnergyEqns()-1;
  
  // Setup the charge array
  for (CFint i = 0; i < _NS; ++i) {
    m_charge[i] = m_gasMixture->speciesCharge(i);
  }
  
  _Rgas = Mutation::RU;
  cf_assert(_Rgas > 0.);
  
  if (m_shiftHO) { 
    /// computation of the H formation at (close to) 0K
    auto_ptr<Mutation::Mixture> gasMixtureLTE; 
    Mutation::MixtureOptions moLTE(_mixtureName);
    moLTE.setStateModel("Equil");
    gasMixtureLTE.reset(new Mutation::Mixture(moLTE)); 
    
    const CFreal P0 = 10.; // pressure can be whatever in this case
    CFreal T0 = 10.; // the equilibrium solver doesn't work with T<10K
    gasMixtureLTE->setState(&P0, &T0, 1);
    // composition at (close to) 0K
    RealVector y0(_NS, const_cast<CFreal*>(gasMixtureLTE->Y()));
    T0 = 0.00000000001; // set T0 closer to 0K
    m_gasMixture->setState(&y0[0], &T0, 1); // suppose rho=1 (rho value is irrelevant at 0K)
    m_H0 = m_gasMixture->mixtureHMass(T0); 
    cf_always_assert(std::abs(m_H0) > 0.);
    CFLog(VERBOSE, "MutationLibrarypp::setup() => y0 = [ " << y0 << "], H0 = [" << m_H0 << "]\n");
    
    // species formation enthalpy at T=0K (this needs to be rechecked for TCNEQ)
    m_gasMixture->speciesHOverRT(&m_vecH0[0], CFNULL, CFNULL, CFNULL, CFNULL);   
    const CFreal RT0 = _Rgas*T0;
    for (CFuint i = 0; i < _NS; ++i) {
      cf_assert(m_molarmassp[i] > 0.);
      const CFreal k = RT0/m_molarmassp[i];
      m_vecH0[i] *= k;
    }
    CFLog(VERBOSE, "MutationLibrarypp::setup() => m_vecH0 = [ " << m_vecH0 << "]\n");
    RealVector yH0(m_vecH0*y0);
    CFLog(VERBOSE, "MutationLibrarypp::setup() => " << yH0.sum() << " == " << m_H0 << "\n");
  }
  
  CFLog(VERBOSE, "MutationLibrarypp::setup() => end\n"); 
}
      
//////////////////////////////////////////////////////////////////////////////

void MutationLibrarypp::unsetup()
{
  CFLog(VERBOSE, "MutationLibrarypp::unsetup() => start\n"); 
  
  if(isSetup()) {
    Framework::PhysicalChemicalLibrary::unsetup();
  }
  
  CFLog(VERBOSE, "MutationLibrarypp::unsetup() => start\n"); 
}
      
//////////////////////////////////////////////////////////////////////////////

CFdouble MutationLibrarypp::lambdaNEQ(CFdouble& temperature,
				      CFdouble& pressure)
{
  CFreal k = m_gasMixture->frozenThermalConductivity();
  // RESET_TO_ZERO(k);
  CFLog(DEBUG_MAX, "Mutation::lambdaNEQ() => k = " << k << "\n");
  return k;
}

//////////////////////////////////////////////////////////////////////////////

void MutationLibrarypp::lambdaVibNEQ(CFreal& temperature,
				     RealVector& tVec,
				     CFdouble& pressure,
				     CFreal& lambdaTrRo,
				     RealVector& lambdaInt)
{
  throw NotImplementedException(FromHere(),"MutationLibrarypp::lambdaVibNEQ()");
}

//////////////////////////////////////////////////////////////////////////////

CFdouble MutationLibrarypp::sigma(CFdouble& temp, //electrical conductivity
				  CFdouble& pressure,
				  CFreal* tVec)
{
  if (temp < 100.) {temp = 100.;}
  // this needs to be modified for NEQ case unless one assumes to have called 
  // setState() before  
  m_gasMixture->setState(&pressure, &temp, 1);
  return m_gasMixture->electricConductivity();
}
      
//////////////////////////////////////////////////////////////////////////////

void MutationLibrarypp::gammaAndSoundSpeed(CFdouble& temp,
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

void MutationLibrarypp::frozenGammaAndSoundSpeed(CFdouble& temp,
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
      
CFdouble MutationLibrarypp::soundSpeed(CFdouble& temp, CFdouble& pressure)
{
  m_gasMixture->setState(&pressure, &temp, 1);
  return m_gasMixture->equilibriumSoundSpeed();
}

//////////////////////////////////////////////////////////////////////////////
      
void MutationLibrarypp::setComposition(CFdouble& temp,
				       CFdouble& pressure,
				       RealVector* x)
{
  if (temp < 100.) {temp = 100.;}

  m_gasMixture->setState(&pressure, &temp, 1);
  const double* xm = m_gasMixture->X();
  
  if (x != CFNULL) {
    for(CFint i = 0; i < _NS; ++i) {
      (*x)[i] = xm[i];
    }
  }
    
  m_gasMixture->convert<X_TO_Y>(xm, &m_y[0]);
  
  CFLog(DEBUG_MAX, "Mutation::setComposition() => m_y = " << m_y << "\n");
}
      
//////////////////////////////////////////////////////////////////////////////

void MutationLibrarypp::setDensityEnthalpyEnergy(CFdouble& temp,
						 CFdouble& pressure,
						 RealVector& dhe)
{
  CFLog(DEBUG_MAX, "Mutation::setDensityEnthalpyEnergy() => P = " 
	<< pressure << ", T = " << temp << "\n");
  
  dhe[0] = m_gasMixture->density();
  dhe[1] = m_gasMixture->mixtureHMass() - m_H0;
  dhe[2] = dhe[1]-pressure/dhe[0];
  
  CFLog(DEBUG_MAX, "Mutation::setDensityEnthalpyEnergy() => " << dhe << ", " <<  m_y << "\n");
}
      
//////////////////////////////////////////////////////////////////////////////

void MutationLibrarypp::setDensityEnthalpyEnergy(CFdouble& temp,
						 RealVector& tVec,
						 CFdouble& pressure,
						 RealVector& dhe,
						 bool storeExtraData)
{
  throw NotImplementedException(FromHere(),"MutationLibrarypp::setDensityEnthalpyEnergy()");
}
      
//////////////////////////////////////////////////////////////////////////////

CFdouble MutationLibrarypp::density(CFdouble& temp,
				    CFdouble& pressure,
				    CFreal* tVec)
{
  if (m_smType == LTE) {m_gasMixture->setState(&pressure, &temp, 1);}
  return m_gasMixture->density();
}

//////////////////////////////////////////////////////////////////////////////

CFdouble MutationLibrarypp::pressure(CFdouble& rho,
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
  
  CFLog(DEBUG_MAX, "Mutation::pressure() => " << p << "\n");
  return p;
}

//////////////////////////////////////////////////////////////////////////////
      
CFdouble MutationLibrarypp::electronPressure(CFreal rhoE,
					     CFreal tempE)
{
  throw NotImplementedException(FromHere(),"MutationLibrarypp::soundSpeed()");
}

//////////////////////////////////////////////////////////////////////////////

CFdouble MutationLibrarypp::energy(CFdouble& temp,
				   CFdouble& pressure)
  
{
  m_gasMixture->setState(&pressure, &temp, 1);
  return m_gasMixture->mixtureEnergyMass()- m_H0;
}
      
//////////////////////////////////////////////////////////////////////////////

CFdouble MutationLibrarypp::enthalpy(CFdouble& temp,
				     CFdouble& pressure)
{
  m_gasMixture->setState(&pressure, &temp, 1);
  return m_gasMixture->mixtureHMass() - m_H0;
}
      
//////////////////////////////////////////////////////////////////////////////

 void MutationLibrarypp::setElemFractions(const RealVector& yn)
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
   throw NotImplementedException(FromHere(),"MutationLibrarypp::setElemFractions()");
 }

//////////////////////////////////////////////////////////////////////////////

 void MutationLibrarypp::setElementXFromSpeciesY(const RealVector& ys)
 {
   for (CFint ic = 0; ic < _NC; ++ic) {
     m_y[ic] = ys[ic];
   }
   //m_gasMixture->convert<Y_TO_XE>(&m_y[0], &m_xn[0]);
   throw NotImplementedException(FromHere(),"MutationLibrarypp::setElementXFromSpeciesY()");
 }

// //////////////////////////////////////////////////////////////////////////////

 void MutationLibrarypp::setElectronFraction(RealVector& ys)
 {
   // charge neutrality: xEl = sum(xIon)
   CFdouble yEl = 0.0;
   for (CFint is = 0; is < _NS; ++is) {
     if (m_charge[is] > 0) {
       yEl += ys[is] / m_molarmassp[is];
     }
   }

   yEl *= m_molarmassp[0]; // 1st species: electron
   ys[0] = yEl; // overwrite electron mass fraction
   
   CFLog(DEBUG_MAX, "Mutation::setElectronFraction() => " << ys[0] << "\n");
 }

//////////////////////////////////////////////////////////////////////////////

void MutationLibrarypp::setSpeciesFractions(const RealVector& ys)
{
  if (m_gasMixture->hasElectrons()) {
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

void MutationLibrarypp::getSpeciesMolarFractions
(const RealVector& ys, RealVector& xs)
{
  RealVector& yss = const_cast<RealVector&>(ys);
  m_gasMixture->convert<Y_TO_X>(&yss[0], &xs[0]);
  CFLog(DEBUG_MAX, "Mutation::getSpeciesMolarFractions() => " << xs << "\n");
}
      
//////////////////////////////////////////////////////////////////////////////

 void MutationLibrarypp::getSpeciesMassFractions
 (const RealVector& xs, RealVector& ys)
 {
   RealVector& xss = const_cast<RealVector&>(xs);
   m_gasMixture->convert<X_TO_Y>(&xss[0], &ys[0]);
   CFLog(DEBUG_MAX, "Mutation::getSpeciesMolarFractions() => " << ys << "\n");
 }

//////////////////////////////////////////////////////////////////////////////

 void MutationLibrarypp::getSpeciesMassFractions(RealVector& ys)
 {
   for (CFint is = 0; is < _NS; ++is) {
     ys[is] = m_y[is];
   }
   CFLog(DEBUG_MAX, "Mutation::getSpeciesMolarFractions() => " << ys << "\n");
 }

//////////////////////////////////////////////////////////////////////////////

void MutationLibrarypp::getTransportCoefs(CFdouble& temp,
					  CFdouble& pressure,
					  CFdouble& lambda,
					  CFdouble& lambdacor,
					  RealVector& lambdael,
					  RealMatrix& eldifcoef,
					  RealVector& eltdifcoef)
{
 throw NotImplementedException(FromHere(),"MutationLibrarypp::getTransportCoefs()");
}

//////////////////////////////////////////////////////////////////////////////

void MutationLibrarypp::getMassProductionTerm(CFdouble& temperature,
					      RealVector& tVec,
					      CFdouble& pressure,
					      CFdouble& rho,
					      const RealVector& ys,
					      bool flagJac,
					      RealVector& omega,
					      RealMatrix& jacobian)
{
  // CFLog(DEBUG_MAX, "MutationLibrarypp::getMassProductionTerm() => P = " 
  // << pressure << ", T = " << temperature << "\n");
  
  // we assume setState() already called before
  if (!_freezeChemistry) {
    m_gasMixture->netProductionRates(&omega[0]);
  } 
  else {
    omega = 0.;
  }
  
  CFLog(DEBUG_MAX, "Mutation::getMassProductionTerm() => omega = " << omega << "\n\n");
  // TODO: Need to figure out how to do the Jacobian later
}
      
//////////////////////////////////////////////////////////////////////////////

void MutationLibrarypp::getSource(CFdouble& temperature,
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
}
      
//////////////////////////////////////////////////////////////////////////////

void MutationLibrarypp::getRhoUdiff(CFdouble& temperature,
				   CFdouble& pressure,
				   RealVector& normConcGradients,
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

void MutationLibrarypp::getDij_fick(RealVector& dx,
				    CFdouble& pressure,
				    CFdouble& temperature,
				    CFreal& Diff_coeff,
				    RealMatrix& Dij,
				    RealVector& rhoUdiff)
{
  throw NotImplementedException(FromHere(),"MutationLibrarypp::getDij_fick()");
}
      
//////////////////////////////////////////////////////////////////////////////

void MutationLibrarypp::getGammaN(CFreal& m_GN)
{
  throw NotImplementedException(FromHere(),"MutationLibrarypp::getGammaN()");
}
//////////////////////////////////////////////////////////////////////////////

void MutationLibrarypp::getGammaO(CFreal& m_GO)
{
  throw NotImplementedException(FromHere(),"MutationLibrarypp::getGammaO()");
}
				  
//////////////////////////////////////////////////////////////////////////////

 void MutationLibrarypp::setSpeciesMolarFractions(const RealVector& xs)
 {
   throw NotImplementedException(FromHere(),"MutationLibrarypp::setSpeciesMolarFractions()");
 }
      
//////////////////////////////////////////////////////////////////////////////

void MutationLibrarypp::getSpeciesTotEnthalpies(CFdouble& temp,
						RealVector& tVec,
						CFdouble& pressure,
						RealVector& hsTot,
						RealVector* hsVib,
						RealVector* hsEl)
{
  // CFLog(DEBUG_MAX, "Mutation::getSpeciesTotEnthalpies()\n");
  
  // recheck this with 2-temperature
  CFreal* hv = (hsVib != CFNULL) ? &(*hsVib)[0] : CFNULL;
  CFreal* he = (hsEl  != CFNULL) ?  &(*hsEl)[0] : CFNULL;
  m_gasMixture->speciesHOverRT(&hsTot[0], CFNULL, CFNULL, hv, he); 
  
  const CFreal RT = _Rgas*temp;
  for (CFuint i = 0; i < _NS; ++i) {
    const CFreal k = RT/m_molarmassp[i];
    hsTot[i] *= k;
    if (hsVib != CFNULL) {hsVib[i] *= k;}
    if (hsEl != CFNULL) {hsEl[i]  *= k;}
  }
  
  hsTot -= m_vecH0; // shift on the formation enthalpy (considering H(T=0K)=0)
  
  CFLog(DEBUG_MAX, "Mutation::getSpeciesTotEnthalpies() => hsTot = " << hsTot << "\n");
}
      
//////////////////////////////////////////////////////////////////////////////

void MutationLibrarypp::getSourceTermVT(CFdouble& temperature,
					RealVector& tVec,
					CFdouble& pressure,
					CFdouble& rho,
					RealVector& omegav,
					CFdouble& omegaRad)
{
  throw NotImplementedException(FromHere(),"MutationLibrarypp::getSourceTermVT()");
}

//////////////////////////////////////////////////////////////////////////////

void MutationLibrarypp::getMolarMasses(RealVector& mm)
{
  assert(mm.size() == static_cast<CFuint>(_NS));

  // check the units
  for (CFint i = 0; i < _NS; ++i) {
    mm[i] = m_molarmassp[i];
  }
  
  CFLog(DEBUG_MAX, "Mutation::getMolarMasses() => " << mm << "\n");
}

//////////////////////////////////////////////////////////////////////////////
      
void MutationLibrarypp::transportCoeffNEQ(CFreal& temperature, 
					 CFdouble& pressure,
					 CFreal* tVec, 
					 RealVector& normConcGradients,
					 CFreal& eta,
					 CFreal& lambdaTrRo, 
					 RealVector& lambdaInt,
					 RealVector& rhoUdiff)
{
   throw NotImplementedException(FromHere(),"MutationLibrarypp::transportCoeffNEQ()");
}
        
//////////////////////////////////////////////////////////////////////////////
   
void MutationLibrarypp::getSourceEE(CFdouble& temperature,
				    RealVector& tVec,
				    CFdouble& pressure,
				    CFdouble& rho,
				    const RealVector& ys,
				    bool flagJac,
				    CFdouble& omegaEE)
{
  throw NotImplementedException(FromHere(),"MutationLibrarypp::getSourceEE()");
}
      
//////////////////////////////////////////////////////////////////////////////
      
} // namespace Mutationpp

} // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

