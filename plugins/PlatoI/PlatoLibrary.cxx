#include "PlatoI/PlatoLibrary.hh"
#include "PlatoI/Plato.hh"
#include <plato_constants_Cpp.h>
#include <plato_fortran_Cpp.h>
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
    m_Xc(),
    m_xn(),
    m_mm()
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

  /*Number of species*/
  _NS = get_nb_species();
    
  m_y.resize(_NS);
  m_x.resize(_NS);
  m_Xc.resize(get_nb_comp());
  m_yn.resize(get_nb_elem());
  m_xn.resize(get_nb_elem());
  m_charge.resize(_NS);
  m_df.resize(_NS);
  m_mm.resize(_NS);

  /*Molar masses*/
  get_mi(&m_mm[0]);
  
  /*Charges*/ 
  get_qi(&m_charge[0]);
  
  /*Perfect gas constant*/
   _Rgas = URU;
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
  cout << "PLATO interface STOPL:: lambdaNEQ\n"; 
  throw NotImplementedException(FromHere(),"PlatoLibrary::lambdaNEQ()");
  return 0.;
}

//////////////////////////////////////////////////////////////////////////////

void PlatoLibrary::lambdaVibNEQ(CFreal& temperature,
				RealVector& tVec,
				CFdouble& pressure,
				CFreal& lambdaTrRo,
				RealVector& lambdaInt)
{
  cout << "PLATO interface STOP:: lambdaVibNEQ\n";
  throw NotImplementedException(FromHere(),"PlatoLibrary::lambdaVibNEQ()");
}

//////////////////////////////////////////////////////////////////////////////

CFdouble PlatoLibrary::sigma(CFdouble& temp, //electrical conductivity
			     CFdouble& pressure,
			     CFreal* tVec)
{
  cout << "PLATO interface STOP:: sigma\n";
  throw NotImplementedException(FromHere(),"PlatoLibrary::sigma()");
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
  cout << "PLATO interface STOP:: gammaAndSoundSpeed\n";
  throw NotImplementedException(FromHere(),"PlatoLibrary::gammaAndSoundSpeed()");
}
      
//////////////////////////////////////////////////////////////////////////////

void PlatoLibrary::frozenGammaAndSoundSpeed(CFdouble& temp,
				            CFdouble& pressure,
					    CFdouble& rho,
				            CFdouble& gamma,
				            CFdouble& soundSpeed,
					    RealVector* tVec)
{
  cout << "PLATO interface STOP:: frozenGammaAndSoundSpeed\n";
  throw NotImplementedException(FromHere(),"PlatoLibrary::frozenGammaAndSoundSpeed()");
}
      
//////////////////////////////////////////////////////////////////////////////
      
CFdouble PlatoLibrary::soundSpeed(CFdouble& temp, CFdouble& pressure)
{
  cout << "PLATO interface STOP:: soundSpeed\n";
  throw NotImplementedException(FromHere(),"PlatoLibrary::soundSpeed()");
  return 0.;
}

//////////////////////////////////////////////////////////////////////////////
      
void PlatoLibrary::setComposition(CFdouble& temp,
				  CFdouble& pressure,
				  RealVector* x)
{
  const size_t flag = 0;

  /*Temperature fix*/
  if (temp < 300.) {temp = 300.;}
  
  /*Compute mole fractions given pressure and temperature*/
  get_eq_composition_mole(&pressure, &temp, &m_Xc[0], &m_x[0], &flag);

  if (x != CFNULL) {
    for(CFint i = 0; i < _NS; ++i) {
      (*x)[i] = static_cast<CFreal>(m_x[i]);
    }
  }

  /*Get mass fractions from mole fractions (call to be implemented!)*/
  cout << "PLATO interface STOP:: setComposition\n";
  throw NotImplementedException(FromHere(),"PlatoLibrary::setComposition()");
}
      
//////////////////////////////////////////////////////////////////////////////

void PlatoLibrary::setDensityEnthalpyEnergy(CFdouble& temp,
				            CFdouble& pressure,
					    RealVector& dhe)
{
  cout << "PLATO interface STOP:: setDensityEnthalpyEnergy\n";
  throw NotImplementedException(FromHere(),"PlatoLibrary::setDensityEnthalpyEnergy()");
}

//////////////////////////////////////////////////////////////////////////////

void PlatoLibrary::setDensityEnthalpyEnergy(CFdouble& temp,
					    RealVector& tVec,
					    CFdouble& pressure,
					    RealVector& dhe,
					    bool storeExtraData)
{
  cout << "PLATO interface STOP:: setDensityEnthalpyEnergy\n";
  throw NotImplementedException(FromHere(),"PlatoLibrary::setDensityEnthalpyEnergy()");
}
      
//////////////////////////////////////////////////////////////////////////////

CFdouble PlatoLibrary::density(CFdouble& temp,
			       CFdouble& pressure,
			       CFreal* tVec)
{
  cout << "PLATO interface STOP:: density\n";
  throw NotImplementedException(FromHere(),"PlatoLibrary::density()");
  return 0.;
}

//////////////////////////////////////////////////////////////////////////////

CFdouble PlatoLibrary::pressure(CFdouble& rho,
				CFdouble& temp,
				CFreal* tVec)
{
  cout << "PLATO interface STOP:: pressure\n";
  throw NotImplementedException(FromHere(),"PlatoLibrary::pressure()");
  return 0.;
}

//////////////////////////////////////////////////////////////////////////////

CFdouble PlatoLibrary::electronPressure(CFreal rhoE,
					CFreal tempE)
{
  cout << "PLATO interface STOP:: soundSpeed\n";
  throw NotImplementedException(FromHere(),"PlatoLibrary::soundSpeed()");
}

//////////////////////////////////////////////////////////////////////////////

CFdouble PlatoLibrary::energy(CFdouble& temp,
			      CFdouble& pressure)
  
{
  cout << "PLATO interface STOP:: energy\n";
  throw NotImplementedException(FromHere(),"PlatoLibrary::energy()");
  return 0.;
}
      
//////////////////////////////////////////////////////////////////////////////

CFdouble PlatoLibrary::enthalpy(CFdouble& temp,
				CFdouble& pressure)
{
  cout << "PLATO interface STOP:: enthalpy\n";
  throw NotImplementedException(FromHere(),"PlatoLibrary::enthalpy()");
  return 0.;
}
      
//////////////////////////////////////////////////////////////////////////////

 void PlatoLibrary::setElemFractions(const RealVector& yn)
 {
   cout << "PLATO interface STOP:: setElemFractions\n";
   //m_gasMixture->convert<YE_TO_XE>(&m_yn[0], &m_xn[0]);
   throw NotImplementedException(FromHere(),"PlatoLibrary::setElemFractions()");
 }

// //////////////////////////////////////////////////////////////////////////////

 void PlatoLibrary::setElementXFromSpeciesY(const RealVector& ys)
 {
   cout << "PLATO interface STOP:: setElementXFromSpeciesY\n";
   //m_gasMixture->convert<Y_TO_XE>(&m_y[0], &m_xn[0]);
   throw NotImplementedException(FromHere(),"PlatoLibrary::setElementXFromSpeciesY()");
 }

// //////////////////////////////////////////////////////////////////////////////

 void PlatoLibrary::setElectronFraction(RealVector& ys)
 {
   cout << "PLATO interface STOP:: setElementXFromSpeciesY\n";
   throw NotImplementedException(FromHere(),"PlatoLibrary::setElementXFromSpeciesY()");
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
   cout << "PLATO interface STOP:: setSpeciesFractions\n";
   throw NotImplementedException(FromHere(),"PlatoLibrary::setSpeciesFractions()");
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
    cout << "PLATO interface STOP:: getSpeciesMolarFractions\n";
    throw NotImplementedException(FromHere(),"PlatoLibrary::getSpeciesMolarFractions()");
    // RealVector& yss = const_cast<RealVector&>(ys);
    //  m_gasMixture->convert<Y_TO_X>(&yss[0], &xs[0]);
 }
      
//////////////////////////////////////////////////////////////////////////////

 void PlatoLibrary::getSpeciesMassFractions
 (const RealVector& xs, RealVector& ys)
 {
    cout << "PLATO interface STOP:: getSpeciesMassFractions\n";
    throw NotImplementedException(FromHere(),"PlatoLibrary::getSpeciesMassFractions()");
    // RealVector& xss = const_cast<RealVector&>(xs);
    // m_gasMixture->convert<X_TO_Y>(&xss[0], &ys[0]);
 }

// //////////////////////////////////////////////////////////////////////////////

 void PlatoLibrary::getSpeciesMassFractions(RealVector& ys)
 {
  cout << "PLATO interface STOP:: getSpeciesMassFractions\n";
  throw NotImplementedException(FromHere(),"PlatoLibrary::getSpeciesMassFractions()");
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
 cout << "PLATO interface STOP:: getTransportCoefs\n";
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
  cout << "PLATO interface STOP:: getMassProductionTerm\n";
  throw NotImplementedException(FromHere(),"PlatoLibrary::getMassProductionTerm()");
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
  cout << "PLATO interface STOP:: getSource\n";
  throw NotImplementedException(FromHere(),"PlatoLibrary::getSource()");
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
  cout << "PLATO interface STOP:: getRhoUDiff\n";  
  throw NotImplementedException(FromHere(),"PlatoLibrary::getRhoUdiff()");
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
  cout << "PLATO interface STOP:: getDij_fick\n";
  throw NotImplementedException(FromHere(),"PlatoLibrary::getDij_fick()");
}
      
//////////////////////////////////////////////////////////////////////////////

void PlatoLibrary::getGammaN(CFreal& m_GN)
{
  cout << "PLATO interface STOP:: getGammaN\n";
  throw NotImplementedException(FromHere(),"PlatoLibrary::getGammaN()");
}
//////////////////////////////////////////////////////////////////////////////

void PlatoLibrary::getGammaO(CFreal& m_GO)
{
  cout << "PLATO interface STOP:: getGammaO\n";
  throw NotImplementedException(FromHere(),"PlatoLibrary::getGammaO()");
}
				  
//////////////////////////////////////////////////////////////////////////////

 void PlatoLibrary::setSpeciesMolarFractions(const RealVector& xs)
 {
   cout << "PLATO interface STOP:: setSpeciesMolarFractions\n";
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
  cout << "PLATO interface STOP:: getSpeciesTotEnthalpies\n";
  throw NotImplementedException(FromHere(),"PlatoLibrary::getSpeciesTotEnthalpies()");
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
  cout << "PLATO interface STOP:: getSourceTermVT\n";
  throw NotImplementedException(FromHere(),"PlatoLibrary::getSourceTermVT()");
}

//////////////////////////////////////////////////////////////////////////////

void PlatoLibrary::getMolarMasses(RealVector& mm)
{
  assert(mm.size() == static_cast<CFuint>(_NS));

  for (CFint i = 0; i < _NS; ++i) {
     mm[i] = m_mm[i];
  }
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
   cout << "PLATO interface STOP:: transportCoeffNEQ\n";
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
  cout << "PLATO interface STOP:: getSourceEE\n";
  throw NotImplementedException(FromHere(),"PlatoLibrary::getSourceEE()");
}
      
//////////////////////////////////////////////////////////////////////////////
      
} // namespace Plato

} // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

