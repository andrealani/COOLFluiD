#include "PlatoI/PlatoLibrary.hh"
#include "PlatoI/Plato.hh"
#include "Common/CFLog.hh"
#include "Environment/ObjectProvider.hh"
#include "Common/BadValueException.hh"
#include "Framework/MeshData.hh"
#include "Framework/PhysicalModel.hh"
#include "Common/Stopwatch.hh"
#include "Environment/DirPaths.hh"
#include "Common/OSystem.hh"
#include "Common/StringOps.hh"
#include <fstream>

#include <plato_fortran_Cpp.h>

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

   namespace Plato {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<PlatoLibrary,
			    PhysicalPropertyLibrary,
			    PlatoModule,
			    1>
platoLibraryProvider("Plato");

//////////////////////////////////////////////////////////////////////////////

void PlatoLibrary::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< std::string >("mixtureName","Name of the mixture.");
  options.addConfigOption< std::string >("reactionName","Name of the reaction.");
  options.addConfigOption< std::string >("transfName","Name of the transfer file.");
}

//////////////////////////////////////////////////////////////////////////////

PlatoLibrary::PlatoLibrary(const std::string& name)
  : Framework::PhysicalChemicalLibrary(name),
    m_y(),
    m_x(),
    m_yn(),
    m_xn(),
    m_molarmassp()
{
  addConfigOptionsTo(this);

  _mixtureName = "";
  setParameter("mixtureName",&_mixtureName);

  _reactionName = "empty";
  setParameter("reactionName",&_reactionName);
  
  _transfName = "empty";
  setParameter("transfName",&_transfName);
}

//////////////////////////////////////////////////////////////////////////////

PlatoLibrary::~PlatoLibrary()
{
}

//////////////////////////////////////////////////////////////////////////////

void PlatoLibrary::configure ( Config::ConfigArgs& args )
{
  Framework::PhysicalChemicalLibrary::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void PlatoLibrary::setup()
{    
  CFLog(VERBOSE, "PlatoLibrary::setup() => start\n"); 
  
  Framework::PhysicalChemicalLibrary::setup();
  
  // if this is a parallel simulation, only ONE process at a time
  // sets the library
  setLibrarySequentially();
  
  CFLog(VERBOSE, "PlatoLibrary::setup() => end\n"); 
}
      
//////////////////////////////////////////////////////////////////////////////

void PlatoLibrary::setLibrarySequentially()
{ 
  const char* solver = "COOLFluiD";
  
  const size_t lsolver   = strlen(solver);
  const size_t lmixture  = strlen(_mixtureName.c_str());
  const size_t lreaction = strlen(_reactionName.c_str());
  const size_t ltransfer = strlen(_transfName.c_str());
  
  if (m_libPath == "") {
    const char* basePath = getenv("PLATO_DIR");
    m_libPath = string(basePath) + "/database";
  }
  
  const size_t lpath = strlen(m_libPath.c_str());
  
  /*Initialize the library*/
  initializeC(solver, _mixtureName.c_str(), _reactionName.c_str(), _transfName.c_str(), 
	      m_libPath.c_str(), lsolver, lmixture, lreaction, ltransfer, lpath);
  
  // _NS = m_gasMixture->nSpecies();
  
  // m_y.resize(_NS);
  // m_x.resize(_NS);
  // m_yn.resize(m_gasMixture->nElements());
  // m_xn.resize(m_gasMixture->nElements());
  // m_charge.resize(_NS);
  // m_df.resize(_NS);
 
  // m_molarmassp.resize(_NS);
  // for (CFint i = 0; i < _NS; ++i) {
  //   m_molarmassp[i] = m_gasMixture->speciesMw(i);
  // }
  
  // // Setup the charge array
  // for (CFint i = 0; i < _NS; ++i) {
  //   m_charge[i] = m_gasMixture->speciesCharge(i);
  // }
  
  // _Rgas = Mutation::RU;
  // cf_assert(_Rgas > 0.);
}
      
//////////////////////////////////////////////////////////////////////////////

void PlatoLibrary::unsetup()
{
  if(isSetup()) {
    Framework::PhysicalChemicalLibrary::unsetup();
  }
}
      
//////////////////////////////////////////////////////////////////////////////

CFdouble PlatoLibrary::lambdaNEQ(CFdouble& temperature,
				 CFdouble& pressure)
{
  return 0.;
}

//////////////////////////////////////////////////////////////////////////////

void PlatoLibrary::lambdaVibNEQ(CFreal& temperature,
				     RealVector& tVec,
				     CFdouble& pressure,
				     CFreal& lambdaTrRo,
				     RealVector& lambdaInt)
{
  throw NotImplementedException(FromHere(),"PlatoLibrary::lambdaVibNEQ()");
}

//////////////////////////////////////////////////////////////////////////////

CFdouble PlatoLibrary::sigma(CFdouble& temp, //electrical conductivity
			     CFdouble& pressure,
			     CFreal* tVec)
{
  if (temp < 100.) {temp = 100.;}
  return 0.;
}
      
//////////////////////////////////////////////////////////////////////////////

void PlatoLibrary::gammaAndSoundSpeed(CFdouble& temp,
				      CFdouble& pressure,
				      CFdouble& rho,
				      CFdouble& gamma,
				      CFdouble& soundSpeed)
{
}
      
//////////////////////////////////////////////////////////////////////////////

void PlatoLibrary::frozenGammaAndSoundSpeed(CFdouble& temp,
						CFdouble& pressure,
						CFdouble& rho,
						CFdouble& gamma,
						CFdouble& soundSpeed,
						RealVector* tVec)
{
}
      
//////////////////////////////////////////////////////////////////////////////
      
CFdouble PlatoLibrary::soundSpeed(CFdouble& temp, CFdouble& pressure)
{
  return 0.;
}

//////////////////////////////////////////////////////////////////////////////
      
void PlatoLibrary::setComposition(CFdouble& temp,
				       CFdouble& pressure,
				       RealVector* x)
{
  if (temp < 100.) {temp = 100.;}
}
      
//////////////////////////////////////////////////////////////////////////////

void PlatoLibrary::setDensityEnthalpyEnergy(CFdouble& temp,
						 CFdouble& pressure,
						 RealVector& dhe)
{
}

//////////////////////////////////////////////////////////////////////////////

void PlatoLibrary::setDensityEnthalpyEnergy(CFdouble& temp,
						 RealVector& tVec,
					    CFdouble& pressure,
					    RealVector& dhe,
					    bool storeExtraData)
{
  throw NotImplementedException(FromHere(),"PlatoLibrary::setDensityEnthalpyEnergy()");
}
      
//////////////////////////////////////////////////////////////////////////////

CFdouble PlatoLibrary::density(CFdouble& temp,
			       CFdouble& pressure,
			       CFreal* tVec)
{
  return 0.;
}

//////////////////////////////////////////////////////////////////////////////

CFdouble PlatoLibrary::pressure(CFdouble& rho,
				CFdouble& temp,
				CFreal* tVec)
{
  return 0.;
}

//////////////////////////////////////////////////////////////////////////////

CFdouble PlatoLibrary::electronPressure(CFreal rhoE,
					CFreal tempE)
{
  throw NotImplementedException(FromHere(),"PlatoLibrary::soundSpeed()");
}

//////////////////////////////////////////////////////////////////////////////

CFdouble PlatoLibrary::energy(CFdouble& temp,
			      CFdouble& pressure)
  
{
  return 0.;
}
      
//////////////////////////////////////////////////////////////////////////////

CFdouble PlatoLibrary::enthalpy(CFdouble& temp,
				     CFdouble& pressure)
{
  return 0.;
}
      
//////////////////////////////////////////////////////////////////////////////

 void PlatoLibrary::setElemFractions(const RealVector& yn)
 {
   //m_gasMixture->convert<YE_TO_XE>(&m_yn[0], &m_xn[0]);
   throw NotImplementedException(FromHere(),"PlatoLibrary::setElemFractions()");
 }

// //////////////////////////////////////////////////////////////////////////////

 void PlatoLibrary::setElementXFromSpeciesY(const RealVector& ys)
 {
   //m_gasMixture->convert<Y_TO_XE>(&m_y[0], &m_xn[0]);
   throw NotImplementedException(FromHere(),"PlatoLibrary::setElementXFromSpeciesY()");
 }

// //////////////////////////////////////////////////////////////////////////////

 void PlatoLibrary::setElectronFraction(RealVector& ys)
 {
   // // charge neutrality: xEl = sum(xIon)
   // CFdouble yEl = 0.0;
   // for (CFint is = 0; is < _NS; ++is) {
   //   if (m_charge[is] > 0) {
   //     yEl += ys[is] / m_molarmassp[is];
   //   }
   // }

   // yEl *= m_molarmassp[0]; // 1st species: electron
   // ys[0] = yEl; // overwrite electron mass fraction
 }

// //////////////////////////////////////////////////////////////////////

 void PlatoLibrary::setSpeciesFractions(const RealVector& ys)
 {
   // if (m_gasMixture->hasElectrons()) {
   //   setElectronFraction(const_cast<RealVector&>(ys));
   // }
  
   // for (CFint is = 0; is < _NS; ++is) {
   //   m_y[is] = ys[is];

   //   if (m_y[is] < 0.0) m_y[is] = 0.0;
   //   cf_assert(m_y[is] < 1.1);
   // }
 }

//////////////////////////////////////////////////////////////////////////////

 void PlatoLibrary::getSpeciesMolarFractions
 (const RealVector& ys, RealVector& xs)
 {
    // RealVector& yss = const_cast<RealVector&>(ys);
    //  m_gasMixture->convert<Y_TO_X>(&yss[0], &xs[0]);
 }
      
//////////////////////////////////////////////////////////////////////////////

 void PlatoLibrary::getSpeciesMassFractions
 (const RealVector& xs, RealVector& ys)
 {
    // RealVector& xss = const_cast<RealVector&>(xs);
    // m_gasMixture->convert<X_TO_Y>(&xss[0], &ys[0]);
 }

// //////////////////////////////////////////////////////////////////////////////

 void PlatoLibrary::getSpeciesMassFractions(RealVector& ys)
 {
   // for (CFint is = 0; is < _NS; ++is) {
   //   ys[is] = m_y[is];
   // }
   // std::cout << "getSpeciesMassFractions() WAS CALLED!!!!" << std::endl;
 }

//////////////////////////////////////////////////////////////////////////////

void PlatoLibrary::getTransportCoefs(CFdouble& temp,
					  CFdouble& pressure,
					  CFdouble& lambda,
					  CFdouble& lambdacor,
					  RealVector& lambdael,
					  RealMatrix& eldifcoef,
					  RealVector& eltdifcoef)
{
 throw NotImplementedException(FromHere(),"PlatoLibrary::getTransportCoefs()");
}

//////////////////////////////////////////////////////////////////////////////

void PlatoLibrary::getMassProductionTerm(CFdouble& temperature,
					      RealVector& tVec,
					      CFdouble& pressure,
					      CFdouble& rho,
					      const RealVector& ys,
					      bool flagJac,
					      RealVector& omega,
					      RealMatrix& jacobian)
{
   // RealVector& yss = const_cast<RealVector&>(ys);
   // CFreal tp[2]; tp[0]= temperature; tp[1]= pressure;
   // m_gasMixture->setState(&tp[0], &yss[0]);
   // m_gasMixture->netProductionRates(&omega[0]);
  // TODO: Need to figure out how to do the Jacobian later
}

//////////////////////////////////////////////////////////////////////////////

void PlatoLibrary::getSource(CFdouble& temperature,
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
   // RealVector& yss = const_cast<RealVector&>(ys);
   // CFreal tp[2]; tp[0]= temperature; tp[1]= pressure;
   // m_gasMixture->setState(&tp[0], &yss[0]);
   // m_gasMixture->netProductionRates(&omega[0]);
}
      
//////////////////////////////////////////////////////////////////////////////

void PlatoLibrary::getRhoUdiff(CFdouble& temperature,
				   CFdouble& pressure,
				   RealVector& normConcGradients,
				   CFreal* tVec,
				   RealVector& rhoUdiff,
				   bool fast)
{  
  // // Set driving forces as gradients of molar fractions
  // CFreal MMass = m_gasMixture->mixtureMw();
  // CFreal normMMassGradient = 0.0;
  // for (CFint is = 0; is < _NS; ++is) {
  //   normMMassGradient += normConcGradients[is] / m_molarmassp[is];
  // }
  // normMMassGradient *= -MMass*MMass;
  
  // for (CFint is = 0; is < _NS; ++is) {
  //   m_df[is] = (MMass*normConcGradients[is] + m_y[is]*normMMassGradient) / m_molarmassp[is];
  // }
  
  // double E = 0.0;
  // m_gasMixture->stefanMaxwell(&m_df[0], &rhoUdiff[0], E);
  
  // CFreal density = m_gasMixture->density();
  // for (CFint is = 0; is < _NS; ++is) {
  //   rhoUdiff[is] *= m_y[is]*density;
  // }
}

//////////////////////////////////////////////////////////////////////////////

void PlatoLibrary::getDij_fick(RealVector& dx,
				    CFdouble& pressure,
				    CFdouble& temperature,
				    CFreal& Diff_coeff,
				    RealMatrix& Dij,
				    RealVector& rhoUdiff)
{
  throw NotImplementedException(FromHere(),"PlatoLibrary::getDij_fick()");
}
      
//////////////////////////////////////////////////////////////////////////////

void PlatoLibrary::getGammaN(CFreal& m_GN)
{
  throw NotImplementedException(FromHere(),"PlatoLibrary::getGammaN()");
}
//////////////////////////////////////////////////////////////////////////////

void PlatoLibrary::getGammaO(CFreal& m_GO)
{
  throw NotImplementedException(FromHere(),"PlatoLibrary::getGammaO()");
}
				  
//////////////////////////////////////////////////////////////////////////////

 void PlatoLibrary::setSpeciesMolarFractions(const RealVector& xs)
 {
   throw NotImplementedException(FromHere(),"PlatoLibrary::setSpeciesMolarFractions()");
 }
      
//////////////////////////////////////////////////////////////////////////////

void PlatoLibrary::getSpeciesTotEnthalpies(CFdouble& temp,
						RealVector& tVec,
						CFdouble& pressure,
						RealVector& hsTot,
						RealVector* hsVib,
						RealVector* hsEl)
{
  // CFreal* hv = (hsVib != CFNULL) ? &(*hsVib)[0] : CFNULL;
  // CFreal* he = (hsEl  != CFNULL) ?  &(*hsEl)[0] : CFNULL;
  // m_gasMixture->speciesHOverRT(&hsTot[0], CFNULL, CFNULL, hv, he);
  // const CFreal RT = _Rgas*temp;
  // for (CFuint i = 0; i < _NS; ++i) {
  //   const CFreal k = RT/m_molarmassp[i];
  //   hsTot[i] *= k;
  //   if (hsVib != CFNULL) {hsVib[i] *= k;}
  //   if (hsEl != CFNULL) {hsEl[i]  *= k;}
  // }
}
      
//////////////////////////////////////////////////////////////////////////////

void PlatoLibrary::getSourceTermVT(CFdouble& temperature,
					RealVector& tVec,
					CFdouble& pressure,
					CFdouble& rho,
					RealVector& omegav,
					CFdouble& omegaRad)
{
  throw NotImplementedException(FromHere(),"PlatoLibrary::getSourceTermVT()");
}

//////////////////////////////////////////////////////////////////////////////

void PlatoLibrary::getMolarMasses(RealVector& mm)
{
  // assert(mm.size() == static_cast<CFuint>(_NS));

  // // check the units
  // for (CFint i = 0; i < _NS; ++i) {
  //   mm[i] = m_molarmassp[i];
  // }
}

//////////////////////////////////////////////////////////////////////////////
      
void PlatoLibrary::transportCoeffNEQ(CFreal& temperature, 
					 CFdouble& pressure,
					 CFreal* tVec, 
					 RealVector& normConcGradients,
					 CFreal& eta,
					 CFreal& lambdaTrRo, 
					 RealVector& lambdaInt,
					 RealVector& rhoUdiff)
{
   throw NotImplementedException(FromHere(),"PlatoLibrary::transportCoeffNEQ()");
}
        
//////////////////////////////////////////////////////////////////////////////
   
void PlatoLibrary::getSourceEE(CFdouble& temperature,
			       RealVector& tVec,
			       CFdouble& pressure,
			       CFdouble& rho,
			       const RealVector& ys,
			       bool flagJac,
			       CFdouble& omegaEE)
{
  throw NotImplementedException(FromHere(),"PlatoLibrary::getSourceEE()");
}
      
//////////////////////////////////////////////////////////////////////////////
      
} // namespace Plato

} // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

