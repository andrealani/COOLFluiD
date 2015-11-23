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
    _Xitol(),
    _Yi(),
    _Xn(),
    _Yn(),
    _mmi(),
    _Ri(),
    _hi(),
    _molIDs(),
    _hiVib(),
    _hiEl(),
    _qi(),
    _rhoi(),
    _rhoe(),
    _tvec(),
    _lambdavec(),
    _prodterm(),
    _jprodterm(),
    _Di(),
    _Dij(),
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
/*!
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
/*!
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

  /*Number of temperatures (total, vibrationa;, free-electron-electronic)*/
  _nTemp  = get_nb_temp();
  _nbTvib = get_nb_vib_temp();
  _nbTe   = get_nb_el_temp();

  /*Number of dimensions*/
  _nDim = PhysicalModelStack::getActive()->getDim();

  /*Number of governing equations*/
  _nEqs = _NS + _nDim + _nTemp;

  /*Allocate working arrays*/
  _Xi.resize(_NS);
  _Xitol.resize(_NS);
  _Yi.resize(_NS);
  _Xc.resize(get_nb_comp());
  _Yc.resize(get_nb_comp());
  _Xn.resize(_NC);
  _Yn.resize(_NC);
  _molIDs.resize(_NS);
  _hi.resize(_NS);
  _hiVib.resize(_NS);
  _hiEl.resize(_NS);
  _mmi.resize(_NS);
  _Ri.resize(_NS);
  _qi.resize(_NS);
  _rhoi.resize(_NS);
  _rhoe.resize(_nTemp);
  _Di.resize(_NS);
  _Dij.resize((get_nb_comp())*((get_nb_comp()) + 1)/2);
  _dfi.resize(_NS);
  _tvec.resize(_nTemp);
  _lambdavec.resize(_nTemp);
  _prodterm.resize(_nEqs);
  _jprodterm.resize(_nEqs,_nEqs);

  /*Call PLATO library functions/subroutines to get molar mass, gas constant, charge and molecule IDs arrays*/
  /*Molar masses [kg/mol]*/
  get_mi(&_mmi[0]);
 
  /*Gas constants [J/(kg*K)]*/
  get_Ri(&_Ri[0]);
 
  /*Charges*/ 
  get_qi(&_qi[0]);
  
  /*Molecular IDs*/
  get_mol_ids(&_molIDs[0]);

  /*Subtract 1 to be consistent with C/C++ arrays*/
  for (CFint i = 0; i < _nMol; ++i) {
     _molIDs[i] -= 1;
  } 

  /*Universal gas constant [J/(mol*K)]*/
   _Rgas = URU;

  /*Set flag to indicate is the mixture is neutral or ionized*/
  _hasElectrons = (get_nb_e() == 1) ? true : false;
}
 
//////////////////////////////////////////////////////////////////////////////     
/*!
 * This function shuts down the PLATO library and frees the memory allocated for local arrays
 */
void PlatoLibrary::unsetup()
{
  if(isSetup()) {

    /*Shut down PLATO library*/
    finalize();

    /*Deallocate working arrays*/
    _Xi.resize(0);
    _Xitol.resize(0);
    _Yi.resize(0);
    _Xc.resize(0);
    _Yc.resize(0);
    _Xn.resize(0);
    _Yn.resize(0);
    _hi.resize(0);
    _molIDs.resize(0);
    _hiVib.resize(0);
    _hiEl.resize(0);
    _mmi.resize(0);
    _Ri.resize(0);
    _qi.resize(0);
    _rhoi.resize(0);
    _rhoe.resize(0);
    _Di.resize(0);
    _Dij.resize(0);
    _dfi.resize(0);
    _tvec.resize(0);
    _lambdavec.resize(0); 
    _prodterm.resize(0);
    _jprodterm.resize(0,0);

    Framework::PhysicalChemicalLibrary::unsetup();
  }
}

//////////////////////////////////////////////////////////////////////////////
/*!
 * This function selects the free-electron temperature
 */     
CFdouble PlatoLibrary::get_el_temp(CFdouble &temp, CFreal* tVec) {

 CFdouble Te;

 Te = temp;
 if (tVec !=NULL) Te = tVec[_nTemp - 2];  

 return Te;
}

//////////////////////////////////////////////////////////////////////////////
/*!
 * This function applies a small tolerance on the mole fractions by ensuring that the 
 * sum of mole fractions is still equal to one.
 */
void PlatoLibrary::comp_tol(RealVector& X, RealVector& Xtol) {
 
  CFdouble tol = 1.e-12;
  CFdouble sum_X = 0.;

  for (CFint i = 0; i < _NS; ++i ) {
    Xtol[i]  = X[i] + tol;
    sum_X   += Xtol[i];
  }
   
  sum_X = 1./sum_X;
  for (CFint i = 0; i < _NS; ++i) {
    Xtol[i] *= sum_X;
  }

}
    
////////////////////////////////////////////////////////////////////////////// 
/*!
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
/*!
 * This function returns the molecular species IDs
 */
void PlatoLibrary::setMoleculesIDs(std::vector<CFuint>& v)
{
 if (_nMol > 0) {
   v.reserve(_nMol);

   for (CFint i = 0; i < _nMol; ++i) {
     v.push_back(_molIDs[i]);
   }

 }
}

//////////////////////////////////////////////////////////////////////////////
/*!
 * This function returns the dynamic viscosity given the pressure and the
 * temperatures (the mole fractions are stored in the vector "_Xi" which has 
 * to be filled before invoking this function)
 */
CFdouble PlatoLibrary::eta(CFdouble& temp, CFdouble& pressure, CFreal* tVec)
{
  /*Heavy-particle and free-electron temperatures*/
  CFdouble Th = temp;
  CFdouble Te = get_el_temp(temp, tVec);
  
  /*Electron mole fraction*/
  CFdouble Xe = _Xi[0];
 
  /*Compute number density*/
  CFdouble nd = get_nb_density(&pressure, &Th, &Te, &Xe);

  /*Compute dynamic viscosity*/
  CFdouble eta = get_dyn_vis(&nd, &Th, &Te, &_Xi[0]);
  
  return eta;
}

//////////////////////////////////////////////////////////////////////////////
/*!
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
/*!
 * This function returns the thermal conductivity (the mole fractions
 * are stored in the vector "_Xi" which has to be filled before invoking this function)
 */
CFdouble PlatoLibrary::lambdaNEQ(CFdouble& temp, CFdouble& pressure)
{
  /*Temperature vector*/
  for (CFint i = 0; i < _nTemp; ++i) {
    _tvec[i] = temp;
  }

  /*Number density*/
  CFdouble nd = pressure/(UKB*temp);

  /*Compute thermal conductivity*/
  get_lambda_vec(&nd, &_Xi[0], &_tvec[0], &_lambdavec[0]);

  return (_lambdavec[0]);
}

//////////////////////////////////////////////////////////////////////////////
/*!
 * This function returns the thermal conductivity components (the mole fractions
 * are stored in the vector "_Xi" which has to be filled before invoking this function)
 */
void PlatoLibrary::lambdaVibNEQ(CFreal& temp, RealVector& tVec,  CFdouble& pressure, CFreal& lambdaTrRo, RealVector& lambdaInt)
{
  /*Heavy-particle and free-electron temperatures*/
  CFdouble Th = temp;
  CFdouble Te = tVec[_nTemp - 2];

  /*Set temperature vector*/
  _tvec[0] = Th;
  for (CFint i = 1; i < _nTemp; ++i) {
    _tvec[i] = tVec[i - 1]; 
  }
 
  /*Electron mole fraction*/
  CFdouble Xe = _Xi[0];

  /*Compute number density*/
  CFdouble nd = get_nb_density(&pressure, &Th, &Te, &Xe);

  /*Compute thermal conductivity components*/
  get_lambda_vec(&nd, &_Xi[0], &_tvec[0], &_lambdavec[0]);

  lambdaTrRo = _lambdavec[0];
  for (CFint i = 1; i < _nTemp; ++i) {
    lambdaInt[i - 1] = _lambdavec[i];
  }
}

//////////////////////////////////////////////////////////////////////////////
/*!
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
  CFdouble sigma = get_sigma_e(&nd, &Th, &Te, &_Xi[0]);

  return sigma;
}

////////////////////////////////////////////////////////////////////////////// 
/*!
 * This function returns the equilibrium speed of sound and specific heat ratio
 * (the mole fractions are stored in the vector "_Xi" which has to be filled 
 * before invoking this function)
 */
void PlatoLibrary::gammaAndSoundSpeed(CFdouble& temp, CFdouble& pressure, CFdouble& rho, CFdouble& gamma, CFdouble& soundSpeed)
{
  get_eq_gamma_sound_speed(&pressure, &temp, &_Xi[0], &gamma, &soundSpeed);  
}

//////////////////////////////////////////////////////////////////////////////      
/*!
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
/*!
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
/*!
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
/*!
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
/*!
 * This function returns the gas density, specific enthalpy and energy, and the specific 
 * non-equilibrium energies given temperature and pressure (the mole fractions are stored in vector "_Xi" which 
 * has to be filled before calling this function)
 */
void PlatoLibrary::setDensityEnthalpyEnergy(CFdouble& temp,
					    RealVector& tVec,
					    CFdouble& pressure,
					    RealVector& dhe,
					    bool storeExtraData)
{
  /*Heavy-particle and free-electron temperatures*/
  CFdouble Th = temp;
  CFdouble Te = get_el_temp(temp, &tVec[0]);
 
  /*Electron mole fraction*/
  CFdouble Xe = _Xi[0];

  /*Compute number  and mass densities*/
  CFdouble nd  = get_nb_density(&pressure, &Th, &Te, &Xe);
  CFdouble rho = get_density(&nd, &_Xi[0]);

  /*Set temperature vector*/
  _tvec[0] = Th;
  for (CFint i = 1; i < _nTemp; ++i) {
    _tvec[i] = tVec[i - 1]; 
  }

  /*Compute energy densities. By passing the mass fractions we obtain quantities per unit mass*/
  get_energy_densities(&_Yi[0], &_tvec[0], &_rhoe[0]);
  CFdouble et = _rhoe[0];

  /*Density, mixture enthalpy and specific energy, and non-equilibrium energies*/
  dhe[0] = rho;
  dhe[1] = et + pressure/rho;
  dhe[2] = et;
  for (CFint i = 1; i < _nTemp; ++i) {
    dhe[2 + i] = _rhoe[i];
  }

  if (storeExtraData) {
     cout << "Extra data not available:: PlatoLibrary::setDensityEnthalpyEnergy()\n";
     throw NotImplementedException(FromHere(),"PlatoLibrary::setDensityEnthalpyEnergy()");
  }
}

//////////////////////////////////////////////////////////////////////////////      
/*!
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
/*!
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
  CFdouble pressure = rho*(get_pressure(&Th, &Te, &_Yc[0]));

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
/*!
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
/*!
 * This function computes the mole fractions given the mass fractions
 */
void PlatoLibrary::getSpeciesMolarFractions(const RealVector& ys, RealVector& xs)
{
  RealVector& yss = const_cast<RealVector&>(ys);

  mass_to_mole_fractions(&yss[0], &xs[0]);
}
      
//////////////////////////////////////////////////////////////////////////////
/*!
 * This function computes the mass fractions given the mole fractions
 */
void PlatoLibrary::getSpeciesMassFractions(const RealVector& xs, RealVector& ys)
{
  RealVector& xss = const_cast<RealVector&>(xs);

  mole_to_mass_fractions(&xss[0], &ys[0]);
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
/*! 
 * This function returns the mass production terms given the density, the mass fractions and 
 * the temperature  
 */
void PlatoLibrary::getMassProductionTerm(CFdouble& temp, RealVector& tVec, CFdouble& pressure,
					 CFdouble& rho, const RealVector& ys, bool flagJac,  
                                         RealVector& omega, RealMatrix& jacobian)
{
  /*Partial densities*/
  for (CFint i = 0; i < _NS; ++i) {
    _rhoi[i] = rho*ys[i];
  }

  /*Temperature vector*/
  _tvec[0] = temp;
  for (CFint i = 1; i < _nTemp; ++i) {
    _tvec[i] = tVec[i - 1];
  }

  /*Compute source term and the related Jacobian with respect to natural variables*/
  if (flagJac) {

    source_jac(&_nDim, &_rhoi[0], &_tvec[0], &_prodterm[0], &_jprodterm[0]);

    /*Transpose matrix (Fortran stores arrays by columns, C/C++ by rows)*/
    for (CFint i = 0; i < _nEqs; ++i) {
      for (CFint j = 0; j < _nEqs; ++j) {
        jacobian(i,j) = _jprodterm(j,i);
      }
    }

  /*Compute source term only*/  
  } else {

    source(&_nDim, &_rhoi[0], &_tvec[0], &_prodterm[0]); 

  }

  /*Fill the arrays as used by COOLFluiD*/
  /*Mass production terms*/
  for (CFint i = 0; i < _NS ; ++i) {
    omega[i] = _prodterm[i];
  }

}

//////////////////////////////////////////////////////////////////////////////
/*! 
 * This function returns the mass production and energy transfer terms given the density, the mass fractions and 
 * the temperature  
 */
void PlatoLibrary::getSource(CFdouble& temp, RealVector& tVec, CFdouble& pressure, CFdouble& rho,
			     const RealVector& ys, bool flagJac, RealVector& omega,
			     RealVector& omegav, CFdouble& omegaRad, RealMatrix& jacobian)   
  
{
  /*Partial densities*/
  for (CFint i = 0; i < _NS; ++i) {
    _rhoi[i] = rho*ys[i];
  }

  /*Temperature vector*/
  _tvec[0] = temp;
  for (CFint i = 1; i < _nTemp; ++i) {
    _tvec[i] = tVec[i - 1];
  }

  /*Compute source term and the related Jacobian with respect to natural variables*/ 
  if (flagJac) {

    source_jac(&_nDim, &_rhoi[0], &_tvec[0], &_prodterm[0], &_jprodterm[0]);

    /*Transpose matrix (Fortran stores arrays by columns, C/C++ by rows)*/
    for (CFint i = 0; i < _nEqs; ++i) {
      for (CFint j = 0; j < _nEqs; ++j) {
        jacobian(i,j) = _jprodterm(j,i);
      }
    }

  /*Compute source term only*/ 
  } else {

    source(&_nDim, &_rhoi[0], &_tvec[0], &_prodterm[0]);

  }

  /*Fill the arrays as used by COOLFluiD*/
  /*Mass production terms*/
  for (CFint i = 0; i < _NS ; ++i) {
    omega[i] = _prodterm[i];
  }

  /*Energy transfer terms*/
  omegaRad = _prodterm[_NS + _nDim];
  for (CFint i = 0; i < (_nTemp - 1); i++) {
    omegav[i] = _prodterm[_NS + _nDim + 1 + i];
  }

}
      
//////////////////////////////////////////////////////////////////////////////
/*!
 * This function computes the mass diffusion fluxes given the temperature, the pressure and the mass fraction 
 * gradients (note that the species mole and mass fraction , "_Xi" and "_Yi" vectors, have to be filled before calling this function)
 */
void PlatoLibrary::getRhoUdiff(CFdouble& temp, CFdouble& pressure, RealVector& normConcGradients, CFreal* tVec, RealVector& rhoUdiff, bool fast)
{
  /*Heavy particle and free-electron temperatures*/
  CFdouble Th = temp;
  CFdouble Te = get_el_temp(temp, tVec);

  /*Electron mole fraction*/
  CFdouble Xe = _Xi[0];

  /*Compute number density*/
  CFdouble nd = get_nb_density(&pressure, &Th, &Te, &Xe);

  /*Compute binary diffusion coefficients*/
  get_bin_diff_coeff(&nd, &Th, &Te, &_Xi[0], &_Dij[0]);

  /*Diffusion driving forces (gradients of mole fractions*/
  if (!fast) {
    CFdouble mm = 0.0;
    CFdouble normMMassGradient = 0.0;
    for (CFint is = 0; is < _NS; ++is) {
        mm += _Xi[is]*_mmi[is];
        normMMassGradient += normConcGradients[is]/_mmi[is];
    }
    normMMassGradient *= -(mm*mm);
  
    for (CFint is = 0; is < _NS; ++is) {
      _dfi[is] = (mm*normConcGradients[is] + _Yi[is]*normMMassGradient)/_mmi[is];
    }
  } else {
    for (CFint is = 0; is < _NS; ++is) {
      _dfi[is] = normConcGradients[is];
    }
  } 
  
  /*Apply tolerance on mole fractions (this step is mandatory for stability)*/
  comp_tol(_Xi, _Xitol);

  /*Compute mass diffusion fluxes (solve the Stefan-Maxwell equations)*/
  get_species_diff_flux(&Th, &Te, &nd, &_Xitol[0], &_Dij[0], &_dfi[0], &rhoUdiff[0]);
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
/*!
 * This function returns the species total, vibrational and electronic enthalpies given the temperatures.
 * Note that the entahlpy of free-electrons is also stored in array hsEl
 */
void PlatoLibrary::getSpeciesTotEnthalpies(CFdouble& temp,
				           RealVector& tVec,
					   CFdouble& pressure,
					   RealVector& hsTot,
					   RealVector* hsVib,
					   RealVector* hsEl)
{
  /*Set temperature vector*/
  _tvec[0] = temp;
  for (CFint i = 1; i < _nTemp; ++i) {
    _tvec[i] = tVec[i - 1];
  }

  /* Thermo-chemical non-equilibrium case*/
  if (_nTemp > 1) {

     /*Total enthalpies*/
     species_tot_vib_el_enthalpy(&_tvec[0], &hsTot[0], &_hiVib[0], &_hiEl[0]);

     /*Vibrational energies*/
     for (CFint i = 0; i < _nMol; ++i) {
       (*hsVib)[i] = _hiVib[_molIDs[i]];
     }

     /*Free-electron/electronic energies*/
     for (CFint i = 0; i < _NS; ++i) {
       (*hsEl)[i] = _hiEl[i];
     }
     
  /*Chemical non-equilibrium case*/
  } else {

    /*Species total enthalpies*/
    species_enthalpy(&_tvec[0], &hsTot[0]);

  }
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

