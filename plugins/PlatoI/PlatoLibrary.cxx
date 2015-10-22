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
    _Xc(),
    _Yc(),
    _Xi(),
    _Yi(),
    _Xn(),
    _Yn(),
    _mmi(),
    _Ri(),
    _hi(),
    _qi(),
    _rhoi(),
    _tvec(),
    _dfi()
{
  addConfigOptionsTo(this);

  _path = "empty";
  setParameter("path",&_path);

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
/*
 * This function takes care of the initialization of the PLATO library
 */
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
/*
 * This function initializes the PLATO library
 */
void PlatoLibrary::setLibrarySequentially()
{ 
  const char* solver = "COOLFluiD";
  
  const size_t lsolver   = strlen(solver);
  const size_t lmixture  = strlen(_mixtureName.c_str());
  const size_t lreaction = strlen(_reactionName.c_str());
  const size_t ltransfer = strlen(_transfName.c_str());
  const size_t lpath     = strlen(_path.c_str());
  
  /*Initialize the PLATO library*/
  initializeC(solver, _mixtureName.c_str(), _reactionName.c_str(), _transfName.c_str(), 
	      _path.c_str(), lsolver, lmixture, lreaction, ltransfer, lpath);

  /*Call PLATO library functions to get parameters about the physical model (e.g. number of species, temperatures)*/
  /*Number of elements (or nuclei)*/
  _NC = get_nb_elem();

  /*Number of species (total, atomic and molecular)*/
  _NS   = get_nb_species();
  _nAt  = get_nb_at_species();
  _nMol = get_nb_mol_species();

  /*Number of temperatures (total, vibrationa, free-electron-electronic)*/
  _nTemp  = get_nb_temp();
  _nbTvib = get_nb_vib_temp();
  _nbTe   = get_nb_el_temp();

  /*Allocate working arrays*/
  _Xi.resize(_NS);
  _Yi.resize(_NS);
  _Xc.resize(get_nb_comp());
  _Yc.resize(get_nb_comp());
  _Xn.resize(_NC);
  _Yn.resize(_NC);
  _hi.resize(_NS);
  _mmi.resize(_NS);
  _Ri.resize(_NS);
  _qi.resize(_NS);
  _rhoi.resize(_NS);
  _dfi.resize(_NS);
  _tvec.resize(_nTemp);

  /*Call PLATO library functions to get molar mass, gas constant and charge arrays*/
  /*Molar masses [kg/mol]*/
  get_mi(&_mmi[0]);
 
  /*Gas constants [J/(kg*K)]*/
  get_Ri(&_Ri[0]);
 
  /*Charges*/ 
  get_qi(&_qi[0]);
  
  /*Universla gas constant [J/(mol*K)]*/
   _Rgas = URU;
}
 
//////////////////////////////////////////////////////////////////////////////     
/*
 * This function shuts down the PLATO library and frees the memory allocated for local arrays
 */
void PlatoLibrary::unsetup()
{
  if(isSetup()) {

    /*Shut down PLATO library*/
    finalize();

    /*Deallocate working arrays*/
    _Xi.resize(0);
    _Yi.resize(0);
    _Xc.resize(0);
    _Yc.resize(0);
    _Xn.resize(0);
    _Yn.resize(0);
    _hi.resize(0);
    _mmi.resize(0);
    _Ri.resize(0);
    _qi.resize(0);
    _rhoi.resize(0);
    _dfi.resize(0);
    _tvec.resize(0);

    Framework::PhysicalChemicalLibrary::unsetup();
  }
}

//////////////////////////////////////////////////////////////////////////////
/*
 * This function selects the free-electron temperature
 */     
CFdouble PlatoLibrary::get_el_temp(CFdouble &temp, CFreal* tVec) {

 CFdouble Te;

 Te = temp;
 if (tVec !=NULL) Te = tVec[_nTemp - 2];  

 return Te;
}
    
////////////////////////////////////////////////////////////////////////////// 
/*
 * This functions returns the specific gas constants [J/(kg*K)]
 */
void PlatoLibrary::setRiGas(RealVector& Ri)
{
  for(CFint is = 0; is < _NS; ++is) {
    Ri[is] = _Ri[is];
  }
}

//////////////////////////////////////////////////////////////////////////////
/*
 * This function returns the molar masses [kg/mol]
 */
void PlatoLibrary::getMolarMasses(RealVector& mm) 
{
  for(CFint is = 0; is < _NS; ++is) {
    mm[is] = _mmi[is];
  }
}

//////////////////////////////////////////////////////////////////////////////
/*
 * This function returns the molecular species IDs
 */
void PlatoLibrary::setMoleculesIDs(std::vector<CFuint>& v)
{
 if (_nMol > 0) {
   v.reserve(_nMol);

   /*Get molecular IDS (FORTRAN array)*/
   get_mol_ids(&v[0]);

   /*Subtract 1 to be consistent with C/C++ arrays*/
   for (CFint i = 0; i < _nMol; ++i) {
     v[i] -= 1;
   }

 }
}

//////////////////////////////////////////////////////////////////////////////
/*
 * This function returns the dynamic viscosity given the pressure and the
 * temperatures (the mole fractions are stored in the vector "_Xi" which has 
 * to be filled before invoking this function)
 */
CFdouble PlatoLibrary::eta(CFdouble& temp, CFdouble& pressure, CFreal* tVec)
{
  cout << "PLATO interface STOP::eta\n";
  throw NotImplementedException(FromHere(),"PlatoLibrary::eta()");
  return 0.;
}

//////////////////////////////////////////////////////////////////////////////
/*
 * This function returns the equilibrium totoal thermal conductivity given the pressure and the
 * temperature (the mole fractions are stored in the vector "_Xi" which has 
 * to be filled before invoking this function)
 */
CFdouble PlatoLibrary::lambdaEQ(CFdouble& temp, CFdouble& pressure) 
{
  cout << "PLATO interface STOP::lambda\n";
  throw NotImplementedException(FromHere(),"PlatoLibrary::lambda()");
  return 0.;
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
/*
 * This function returns the electrical conductivity given the pressure and the
 * temperatures (the mole fractions are stored in the vector "_Xi" which has 
 * to be filled before invoking this function)
 */
CFdouble PlatoLibrary::sigma(CFdouble& temp, CFdouble& pressure, CFreal* tVec)
{
  /*Heavy-particle and free-electron temperatures*/
  CFdouble Th = temp;
  CFdouble Te = get_el_temp(temp, tVec);
  
  /*Electron mole fraction*/
  CFdouble Xe = _Xi[0];
 
  /*Compute number density*/
  CFdouble nd = get_nb_density(&pressure, &Th, &Te, &Xe);

  /*Compute electrical conductivity*/
  CFdouble sigma = compute_sigma_e(&nd, &Th, &Te, &_Xi[0]);

  return sigma;
}

////////////////////////////////////////////////////////////////////////////// 
/*
 * This function returns the equilibrium speed of sound and specific heat ratio
 * (the mole fractions are stored in the vector "_Xi" which has to be filled 
 * before invoking this function)
 */
void PlatoLibrary::gammaAndSoundSpeed(CFdouble& temp, CFdouble& pressure, CFdouble& rho, CFdouble& gamma, CFdouble& soundSpeed)
{
  get_eq_gamma_sound_speed(&pressure, &temp, &_Xi[0], &gamma, &soundSpeed);  
}

//////////////////////////////////////////////////////////////////////////////      
/*
 * This function returns the frozen speed of sound and the specific heat ratio 
 * (the mole fractions are stored in the vector "_Xi" which has to be filled before invoking this functio)
 */
void PlatoLibrary::frozenGammaAndSoundSpeed(CFdouble& temp, CFdouble& pressure, CFdouble& rho, CFdouble& gamma, CFdouble& soundSpeed, RealVector* tVec)
{
  /*Temperature vector*/
  _tvec[0] = temp;
  for (CFint i = 1; i < _nTemp; ++i) {
    _tvec[i] = (*tVec)[i - 1];
  }

  /*Partial densities*/
  for (CFint i = 0; i < _NS; ++i) {
    _rhoi[i] = rho*_Yi[i];
  }

  /*Get frozen specific heat ratio and speed of sound*/  
  get_frozen_gamma_sound_speed(&pressure, &rho, &_rhoi[0], &_tvec[0], &gamma, &soundSpeed);
}
      
//////////////////////////////////////////////////////////////////////////////
/*
 * This function returns the equilibrium speed of sound given the pressure and the temperature
 * (the mole fractions are stored in the vector "_Xi" which has to be filled before invoking this function)
 */      
CFdouble PlatoLibrary::soundSpeed(CFdouble& temp, CFdouble& pressure)
{
  CFdouble gamma, soundSpeed;

  /*Compute equilibrium sound speed*/
  get_eq_gamma_sound_speed(&pressure, &temp, &_Xi[0], &gamma, &soundSpeed);

  return soundSpeed;
}

//////////////////////////////////////////////////////////////////////////////
/*
 * This function returns the equilibrium chemical composition given the pressure and 
 * temperature (the mole and mass fractions computed here are stored in the arrays 
 * "_Xi" and "_Yi" which have to be filled before calling this function) 
 */
void PlatoLibrary::setComposition(CFdouble& temp, CFdouble& pressure, RealVector* x)
{
  const size_t flag = 0;

  /*Temperature fix*/
  CFdouble T = temp;
  if (temp < 300.) {T = 300.;}
  
  /*Compute mole fractions given pressure and temperature*/
  get_eq_composition_mole(&pressure, &T, &_Xc[0], &_Xi[0], &flag);

  if (x != CFNULL) {
    for(CFint i = 0; i < _NS; ++i) {
      (*x)[i] = static_cast<CFreal>(_Xi[i]);
    }
  }
 
  /*Get mass fractions from mole fractions (which will be used later)*/
  mole_to_mass_fractions(&_Xi[0], &_Yi[0]);

  /*Get mass fractions of chemical components*/
  get_comp_fractions(&_Yi[0], &_Yc[0]);
}
 
//////////////////////////////////////////////////////////////////////////////     
/*
 * This function returns the gas density, specific enthalpy and energy given 
 * temperature and pressure (the mole fractions are stored in vector "_Xi" which 
 * has to be filled before calling this function)
 */
void PlatoLibrary::setDensityEnthalpyEnergy(CFdouble& temp, CFdouble& pressure, RealVector& dhe)
{
  /*Set temperature vector*/
  for (CFint i = 0; i < _nTemp; ++i) {
    _tvec[i] = temp;
  }

  /*Compute species enthalpies*/
  species_enthalpy(&_tvec[0], &_hi[0]);
  
  /*Compute gas specific enthalpy and density*/
  CFdouble rho = 0.;  CFdouble h = 0.; 
  for (CFint i = 0; i < _NS; ++i) {
    h   += _Yi[i]*_hi[i];
    rho += _Xi[i]*_mmi[i];
  }
  rho *= pressure/(_Rgas*temp);

  /*Density, specific enthalpy and specific energy*/
  dhe[0] = rho;
  dhe[1] = h;
  dhe[2] = h - pressure/rho;
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
/*
 * This function returns the density given the pressure and temperatures (the mole fractions 
 * are stored in vector "_Xi" which has to be filled before calling this function)
 */
CFdouble PlatoLibrary::density(CFdouble& temp, CFdouble& pressure, CFreal* tVec)
{
  /*Heavy-particle and free-electron temperatures*/
  CFdouble Th = temp;
  CFdouble Te = get_el_temp(temp, tVec);
 
  /*Electron mole fraction*/
  CFdouble Xe = _Xi[0];

  /*Compute number density*/
  CFdouble nd = get_nb_density(&pressure, &Th, &Te, &Xe);

  /*Compute density*/
  CFdouble rho = get_density(&nd, &_Xi[0]);

  return rho;
}

//////////////////////////////////////////////////////////////////////////////
/*
 * This function returns the gas pressure given the density and the temperatures
 * (note that the component mass fractions "_Yc" vector has to be filled before calling this function)
 */
CFdouble PlatoLibrary::pressure(CFdouble& rho, CFdouble& temp, CFreal* tVec)
{
  /*Heavy-particle and free-electron temperatures*/
  CFdouble Th = temp;
  CFdouble Te = get_el_temp(temp, tVec);
  
  /*Compute pressure (note that we multiply by the density since the PLATO function
   *assumes as input arguments the partial densities of the chemical components)*/
  CFdouble pressure = rho*get_pressure(&Th, &Te, &_Yc[0]);
  
  return pressure;
}

//////////////////////////////////////////////////////////////////////////////
CFdouble PlatoLibrary::electronPressure(CFreal rhoE,
					CFreal tempE)
{
  cout << "PLATO interface STOP:: soundSpeed\n";
  throw NotImplementedException(FromHere(),"PlatoLibrary::soundSpeed()");
  return 0.;
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
   //m_gasMixture->convert<YE_TO_XE>(&_Yn[0], &_Xn[0]);
   throw NotImplementedException(FromHere(),"PlatoLibrary::setElemFractions()");
 }

// //////////////////////////////////////////////////////////////////////////////
 void PlatoLibrary::setElementXFromSpeciesY(const RealVector& ys)
 {
   cout << "PLATO interface STOP:: setElementXFromSpeciesY\n";
   //m_gasMixture->convert<Y_TO_XE>(&_Yi[0], &m_xn[0]);
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
   //   if (_qi[is] > 0) {
   //     yEl += ys[is] / m_molarmassp[is];
   //   }
   // }

   // yEl *= m_molarmassp[0]; // 1st species: electron
   // ys[0] = yEl; // overwrite electron mass fraction
 }

//////////////////////////////////////////////////////////////////////////////
/*
 * This function fills the mass and mole fractions working vectors "_Yi, "_Yc", "_Xi" and "_Xc", which are 
 * used in other function calls, given the (species) mass fractions "ys" 
 */
void PlatoLibrary::setSpeciesFractions(const RealVector& ys)
{
  /*Set mass fractions (species)*/
  for (CFint is = 0; is < _NS; ++is) {
    _Yi[is] = ys[is];
    /*Fix application*/
    if (_Yi[is] < 0.0) _Yi[is] = 1.e-12;
       cf_assert(_Yi[is] < 1.1);
  }

  /*Set mass fractions of chemical components*/
  get_comp_fractions(&_Yi[0], &_Yc[0]);

  /*Compute mole fractions based on mass fractions*/
  mass_to_mole_fractions(&_Yi[0], &_Xi[0]);

  /*Get mole fractions of chemical components*/
  get_comp_fractions(&_Xi[0], &_Xc[0]);
}

//////////////////////////////////////////////////////////////////////////////
void PlatoLibrary::getSpeciesMolarFractions(const RealVector& ys, RealVector& xs)
{
  RealVector& yss = const_cast<RealVector&>(ys);

  mass_to_mole_fractions(&yss[0], &xs[0]);

  cout << "PLATO interface STOP:: getSpeciesMolarFractions\n";
  throw NotImplementedException(FromHere(),"PlatoLibrary::getSpeciesMolarFractions()");
}
      
//////////////////////////////////////////////////////////////////////////////
void PlatoLibrary::getSpeciesMassFractions(const RealVector& xs, RealVector& ys)
{
  RealVector& xss = const_cast<RealVector&>(xs);

  mole_to_mass_fractions(&xss[0], &ys[0]);

  cout << "PLATO interface STOP:: getSpeciesMassFractions\n";
  throw NotImplementedException(FromHere(),"PlatoLibrary::getSpeciesMassFractions()");
}

////////////////////////////////////////////////////////////////////////////////
void PlatoLibrary::getSpeciesMassFractions(RealVector& ys)
{
  for (CFint is = 0; is < _NS; ++is) {
    ys[is] = _Yi[is];
   }
  cout << "PLATO interface STOP:: getSpeciesMassFractions\n";
  throw NotImplementedException(FromHere(),"PlatoLibrary::getSpeciesMassFractions()");
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
  /*Partial densities*/

  /*Temperature vector*/

  /*Compute source term*/
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
  //   m_df[is] = (MMass*normConcGradients[is] + _Yi[is]*normMMassGradient) / m_molarmassp[is];
  // }
  
  // double E = 0.0;
  // m_gasMixture->stefanMaxwell(&m_df[0], &rhoUdiff[0], E);
  
  // CFreal density = m_gasMixture->density();
  // for (CFint is = 0; is < _NS; ++is) {
  //   rhoUdiff[is] *= _Yi[is]*density;
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
void PlatoLibrary::getTransportCoefs(CFdouble& temp,
			             CFdouble& pressure,
			             CFdouble& lambda,
			             CFdouble& lambdacor,
			             RealVector& lambdael,
			             RealMatrix& eldifcoef,
			             RealVector& eltdifcoef)
{
  cout << "PLATO interface STOP::getTransportCoefs \n";
  throw NotImplementedException(FromHere(),"PlatoLibrary::getTransportCoefs()");
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

