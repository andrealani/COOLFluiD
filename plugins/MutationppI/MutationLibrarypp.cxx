#include "MutationppI/MutationLibrarypp.hh"
#include "MutationppI/Mutationpp.hh"
#include "Common/CFLog.hh"
#include "Environment/ObjectProvider.hh"
#include "Common/BadValueException.hh"
#include "Framework/MeshData.hh"
#include "Framework/PhysicalModel.hh"
#include "Common/Stopwatch.hh"
#include "Environment/DirPaths.hh"
#include "Common/OSystem.hh"
#include "Common/PE.hh"
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
}

//////////////////////////////////////////////////////////////////////////////

MutationLibrarypp::MutationLibrarypp(const std::string& name)
  : Framework::PhysicalChemicalLibrary(name),
    m_gasMixture(CFNULL),
    m_y(),
    m_x(),
    m_yn(),
    m_xn(),
    m_molarmassp()
{
  addConfigOptionsTo(this);
  
  _mixtureName = "";
  setParameter("mixtureName",&_mixtureName);
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
  Framework::PhysicalChemicalLibrary::setup();
  
  // if this is a parallel simulation, only ONE process at a time
  // sets the library
  const std::string nsp = MeshDataStack::getActive()->getPrimaryNamespace();
  
  if (PE::GetPE().IsParallel()) {
    PE::GetPE().setBarrier(nsp);
    
    for (CFuint i = 0; i < PE::GetPE().GetProcessorCount(nsp); ++i) {
      if (i == PE::GetPE().GetRank (nsp)) {
        setLibrarySequentially2();
      }
      
      PE::GetPE().setBarrier(nsp);
    }
  }
  else {
    setLibrarySequentially2();
  }
}
      
//////////////////////////////////////////////////////////////////////////////

void MutationLibrarypp::setLibrarySequentially2()
{ 
  Mutation::MixtureOptions mo(_mixtureName);
  mo.setStateModel("Equil"); // equilibrium
  m_gasMixture.reset(new Mutation::Mixture(mo));
  
  _NS = m_gasMixture->nSpecies();
  
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
  
  // Setup the charge array
  for (CFint i = 0; i < _NS; ++i) {
    m_charge[i] = m_gasMixture->speciesCharge(i);
  }
}

//////////////////////////////////////////////////////////////////////////////

void MutationLibrarypp::unsetup()
{
  if(isSetup()) {
    Framework::PhysicalChemicalLibrary::unsetup();
  }
}
      
//////////////////////////////////////////////////////////////////////////////

CFdouble MutationLibrarypp::lambdaNEQ(CFdouble& temperature,
				      CFdouble& pressure)
{
  return m_gasMixture->frozenThermalConductivity();
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
  m_gasMixture->setState(&pressure, &temp, 1);
  return m_gasMixture->sigma();
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
  soundSpeed = m_gasMixture->frozenSoundSpeed();
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
}

//////////////////////////////////////////////////////////////////////////////

void MutationLibrarypp::setDensityEnthalpyEnergy(CFdouble& temp,
						 CFdouble& pressure,
						 RealVector& dhe)
{
  dhe[0] = m_gasMixture->density();
  dhe[1] = m_gasMixture->mixtureHMass() - m_gasMixture->mixtureHMass(298.15);
  dhe[2] = m_gasMixture->mixtureEnergyMass();
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
  m_gasMixture->setState(&pressure, &temp, 1);
  return m_gasMixture->density();
}

//////////////////////////////////////////////////////////////////////////////

CFdouble MutationLibrarypp::pressure(CFdouble& rho,
				     CFdouble& temp,
				     CFreal* tVec)
{
  return m_gasMixture->pressure(temp, rho, &m_y[0]);
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
  return m_gasMixture->mixtureEnergyMass();
}
      
//////////////////////////////////////////////////////////////////////////////

CFdouble MutationLibrarypp::enthalpy(CFdouble& temp,
				     CFdouble& pressure)
{
  m_gasMixture->setState(&pressure, &temp, 1);
  return m_gasMixture->mixtureHMass() - m_gasMixture->mixtureHMass(298.15);
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

// //////////////////////////////////////////////////////////////////////////////

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
 }

// //////////////////////////////////////////////////////////////////////

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
 }

// //////////////////////////////////////////////////////////////////////////////

 void MutationLibrarypp::getSpeciesMolarFractions
 (const RealVector& ys, RealVector& xs)
 {
    RealVector& yss = const_cast<RealVector&>(ys);
     m_gasMixture->convert<Y_TO_X>(&yss[0], &xs[0]);
 }
      
// //////////////////////////////////////////////////////////////////////////////

 void MutationLibrarypp::getSpeciesMassFractions
 (const RealVector& xs, RealVector& ys)
 {
    RealVector& xss = const_cast<RealVector&>(xs);
    m_gasMixture->convert<X_TO_Y>(&xss[0], &ys[0]);
 }

// //////////////////////////////////////////////////////////////////////////////

 void MutationLibrarypp::getSpeciesMassFractions(RealVector& ys)
 {
   for (CFint is = 0; is < _NS; ++is) {
     ys[is] = m_y[is];
   }
   std::cout << "getSpeciesMassFractions() WAS CALLED!!!!" << std::endl;
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
   RealVector& yss = const_cast<RealVector&>(ys);
   CFreal tp[2]; tp[0]= temperature; tp[1]= pressure;
   m_gasMixture->setState(&tp[0], &yss[0]);
   m_gasMixture->netProductionRates(&omega[0]);
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
   RealVector& yss = const_cast<RealVector&>(ys);
   CFreal tp[2]; tp[0]= temperature; tp[1]= pressure;
   m_gasMixture->setState(&tp[0], &yss[0]);
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
  
  double E = 0.0;
  m_gasMixture->stefanMaxwell(&m_df[0], &rhoUdiff[0], E);
  
  CFreal density = m_gasMixture->density();
  for (CFint is = 0; is < _NS; ++is) {
    rhoUdiff[is] *= m_y[is]*density;
  }
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

