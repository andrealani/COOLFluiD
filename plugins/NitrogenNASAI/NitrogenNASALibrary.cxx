#include "NitrogenNASAI/NitrogenNASALibrary.hh"
#include "NitrogenNASAI/NitrogenNASA.hh"
#include "Common/CFLog.hh"
#include "Environment/ObjectProvider.hh"
#include "Common/BadValueException.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/MeshData.hh"
#include "Common/Stopwatch.hh"
#include "Environment/DirPaths.hh"
#include "Common/OSystem.hh"
#include "Common/PE.hh"
#include "Common/StringOps.hh"
#include <cmath>

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace NitrogenNASA {

//////////////////////////////////////////////////////////////////////////////

extern "C" {

  /**
   * Importing of the FORTRAN subroutine N2INITIALIZE
   */
  void FORTRAN_NAME(n2initialize)(FINT, FINT);

  /**
   * Importing of the FORTRAN subroutine N2LOOKUP
   */
  void FORTRAN_NAME(n2lookup)();
     
  /**
   * Importing of the FORTRAN subroutine N2CV
   */
  void FORTRAN_NAME(n2cv)(FDOUBLE, FDOUBLE);    

  /**
   * Importing of the FORTRAN subroutine N2CP
   */
  void FORTRAN_NAME(n2cp)(FDOUBLE, FDOUBLE);   
  
  /**
   * Importing of the FORTRAN subroutine N2ENERGY
   */
  void FORTRAN_NAME(n2energy)(FDOUBLE, FDOUBLE, FDOUBLE, FDOUBLE, FDOUBLE, FDOUBLE, FDOUBLE);   
  
  /**
   * Importing of the FORTRAN subroutine N2ENTHALPY
   */
  void FORTRAN_NAME(n2enthalpy)(FDOUBLE, FDOUBLE, FDOUBLE, FDOUBLE, FDOUBLE, FDOUBLE, FDOUBLE); 

  /**
   * Importing of the FORTRAN subroutine FROZENGAMMAN2
   */
  void FORTRAN_NAME(frozengamman2)(FDOUBLE, FDOUBLE, FDOUBLE);   

  /**
   * Importing of the FORTRAN subroutine SOUNDSPEED
   */
  void FORTRAN_NAME(soundspeed)(FDOUBLE, FDOUBLE, FDOUBLE, FDOUBLE); 
 
  /**
   * Importing of the FORTRAN subroutine OMEGAN2
   */
  void FORTRAN_NAME(omegan2)(FDOUBLE, FDOUBLE, FDOUBLE, FDOUBLE, FDOUBLE, FDOUBLE);   
                              
  /**
  * Importing of the FORTRAN subroutine EQCOMPN2	
  */
  void FORTRAN_NAME(eqcompn2)(FDOUBLE, FDOUBLE, FDOUBLE, FDOUBLE);  

  /**
  * Importing of the FORTRAN subroutine MOLARMASSES	
  */
  void FORTRAN_NAME(molarmasses)(FDOUBLE);  

  /**
  * Importing of the FORTRAN subroutine UNGASCONST	
  */
  void FORTRAN_NAME(ungasconst)(FDOUBLE);  

  /**
  * Importing of the FORTRAN subroutine SPECMASSTOMOLFRAC	
  */
  void FORTRAN_NAME(specmasstomolfrac)(FDOUBLE, FDOUBLE);
  
  /**
  * Importing of the FORTRAN subroutine DENSITY	
  */
  void FORTRAN_NAME(density)(FDOUBLE, FDOUBLE, FDOUBLE);

  /**
  * Importing of the FORTRAN subroutine NUMBERD	
  */
  void FORTRAN_NAME(numberd)(FDOUBLE, FDOUBLE, FDOUBLE);

  /**
  * Importing of the FORTRAN subroutine PRESSURE	
  */
  void FORTRAN_NAME(pressure)(FDOUBLE, FDOUBLE, FDOUBLE, FDOUBLE);

  /**
  * Importing of the FORTRAN subroutine MOLARMASS	
  */
  void FORTRAN_NAME(molarmass)(FDOUBLE, FDOUBLE, FDOUBLE);

}

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<NitrogenNASALibrary,
          PhysicalPropertyLibrary,
          NitrogenNASAModule,
          1>
nitrogenNASALibraryProvider("NitrogenNASA");

//////////////////////////////////////////////////////////////////////////////

void NitrogenNASALibrary::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< std::string >("path","Library path.");

  options.addConfigOption< CFdouble >("Xlim","Species molar fractions limit.");

  options.addConfigOption< CFdouble >("TminFix","Minimum allowed temperature for chemistry.");

  options.addConfigOption< std::string >("PhysModel","Physical model for nonequilibrium conditions.");
 
  options.addConfigOption< std::string >("FlagVibSTS","Physical model for nonequilibrium conditions.");
}

//////////////////////////////////////////////////////////////////////////////

NitrogenNASALibrary::NitrogenNASALibrary(const std::string& name)
  : PhysicalChemicalLibrary(name)
{
  addConfigOptionsTo(this);

  _libPath = "";
  setParameter("path",&_libPath);

  _Xlim = 1.0e-20;
  setParameter("Xlim",&_Xlim);

  _TminFix =100.;
  setParameter("TminFix",&_TminFix);

  _libPhysModel = "";
  setParameter("PhysModel",&_libPhysModel);

  _libFlagVibSTS = "";
  setParameter("FlagVibSTS",&_libFlagVibSTS);
  
  _EPS = 1e-6; // this value gives problems with air11 (composition, speed of sound) and LTE
  _TOL = 1e-9; // recommended to replace with the use of Xlim and COMPOTOL2
}

//////////////////////////////////////////////////////////////////////////////

NitrogenNASALibrary::~NitrogenNASALibrary()
{

}

//////////////////////////////////////////////////////////////////////////////

void NitrogenNASALibrary::configure ( Config::ConfigArgs& args )
{
  PhysicalChemicalLibrary::configure(args);

  if (_lambdaAlgoStr == "CG") _lambdaAlgo = LAMBDACG;
  if (_lambdaAlgoStr == "Direct") _lambdaAlgo = LAMBDAD;
  if (_etaAlgoStr == "CG") _etaAlgo = ETACG;
  if (_etaAlgoStr == "Direct") _etaAlgo = ETAD;
}

//////////////////////////////////////////////////////////////////////////////

void NitrogenNASALibrary::setup()
{
  PhysicalChemicalLibrary::setup();
  
  // if this is a parallel simulation, only ONE process at a time
  // sets the library
  copyDataFiles();
  setLibrarySequentially();
  deleteDataFiles();
}
      
//////////////////////////////////////////////////////////////////////////////

void NitrogenNASALibrary::copyDataFiles()
{
  // if nobody has configured the library path, it is done here
  // with default = basedir + "plugins/NitrogenNASAI/data/nitrogenNASA/"
  if (_libPath == "") {
    CFLog(NOTICE, "NitrogenNASALibrary::libpath set to default" << "\n");
    std::string baseDir = Environment::DirPaths::getInstance().getBaseDir().string();
    _libPath = baseDir + "/plugins/NitrogenNASAI/";
  }

  std::string command1 = "cp -R " + _libPath + "data/ .";
  Common::OSystem::getInstance().executeCommand(command1);


}

//////////////////////////////////////////////////////////////////////////////

void NitrogenNASALibrary::deleteDataFiles()
 {
  std::string command2 = "rm -fr data/ nitrogenNASA.in";
  Common::OSystem::getInstance().executeCommand(command2);
}

//////////////////////////////////////////////////////////////////////////////

void NitrogenNASALibrary::setLibrarySequentially()
{
  // Library inizialization
  // First one has to check which physical model is using
  if (_libPhysModel == "vibSTS") {
     if ((_libFlagVibSTS != "NASA") && (_libFlagVibSTS != "Capitelli") ) {
        cout << "NitrogenNASALibrary::Error in the option FlagVibSTS " << endl;
        cout << "NitrogenNASALibrary::please change the CFcase file " << endl;
        abort();
     } else if (_libFlagVibSTS == "NASA" ) {  
       _NS = 62; 
     } else if (_libFlagVibSTS == "Capitelli") { 
       _NS = 69;
     } 
     _nbTvib = 0;
  } else {
     _NS = 2;
     _nbTvib = 1; 
  }
  cout << "NitrogenNASALibrary::Physical model: " << _libPhysModel << endl;
  if (_libFlagVibSTS != "") { 
     cout << "NitrogenNASALibrary::Model for the rates: " << _libFlagVibSTS << endl;
  } 
  cout << "NitrogenNASALibrary::Number of species: " << _NS << endl;
  cout << "NitrogenNASALibrary::Number of vibrational temperatures: " << _nbTvib << endl;

  // Initializing the Nitrogen library
  FORTRAN_NAME(n2initialize)(&_NS, &_nbTvib);

  // Number of nuclei (N-N2 system)
  _NC = 1;

  _X    =  new CFdouble[_NS];
  _XTOL =  new CFdouble[_NS];
  _Xini =  new CFdouble[_NS];
  _Xn   =  new CFdouble[_NC];
  _Yn   =  new CFdouble[_NC];
  _Y    =  new CFdouble[_NS];
  _Yini =  new CFdouble[_NS];

  // universal gas constant
  FORTRAN_NAME(ungasconst)(&_Rgas);

  // initialization
  for(CFint i = 0; i < _NS; ++i) {
    _X[i] = 1.0;
    _Y[i] = 1.0;
  }

  _OMEGA = new CFdouble[_NS];
  _MOLARMASSP = new CFdouble[_NS];

   // species molar masses
  FORTRAN_NAME(molarmasses)(_MOLARMASSP);

  const CFuint sizeJ = _NS + _nbTvib + _nbTe + 1;
  const CFuint sizeTV = _nbTvib;

  _OMEGAJACOB = new CFdouble[_NS*sizeJ];
  _OMEGAVTJACOB = new CFdouble[(_nbTvib + _nbTe)*(sizeJ-1)];
  
  _HTOTALUM =  new CFdouble[_NS];
  _HTRANSUM =  new CFdouble[_NS];
  _HELECTUM =  new CFdouble[_NS];
  _HROTUM   =  new CFdouble[_NS];
  _HVIBRUM  =  new CFdouble[_NS];
  _HFORMUM  =  new CFdouble[_NS];
  _OMEGAVIB  =  new CFdouble[sizeTV];

  _CV =  new CFdouble[_NS];
  _CP =  new CFdouble[_NS];

  _extraData.enthalpyTt.resize(_NS);
  _extraData.energyTr.resize(_NS);
  _extraData.energyVib.resize(_NS);
  _extraData.cpVib.resize(_nbTvib);
  _extraData.enthalpyForm.resize(_NS);
  _extraData.dRhoEdRhoi.resize(_NS);
  _extraData.dRhoEvdRhoi.resize(_NS);
  _atomicityCoeff.resize(_NS);

  PhysicalChemicalLibrary::setup();

  // Calling the N2LOOKUp subroutine in order to 
  // perform bi-cubic spline fitting
  // (just for NASA Ames database)
  if ((_libPhysModel == "vibSTS") && (_libFlagVibSTS == "NASA")) {
     FORTRAN_NAME(n2lookup)();
  }
     
}

//////////////////////////////////////////////////////////////////////////////

void NitrogenNASALibrary::unsetup()
{
  if(isSetup()) {

    deletePtrArray(_X);
    deletePtrArray(_XTOL);
    deletePtrArray(_Xini);
    deletePtrArray(_Xn);
    deletePtrArray(_Yn);
    deletePtrArray(_Y);
    deletePtrArray(_Yini);

    deletePtrArray(_OMEGA);
    deletePtrArray(_MOLARMASSP);
    deletePtrArray(_OMEGAJACOB);

    deletePtrArray(_HTOTALUM);
    deletePtrArray(_HTRANSUM);
    deletePtrArray(_HELECTUM);
    deletePtrArray(_HROTUM);
    deletePtrArray(_HVIBRUM);
    deletePtrArray(_HFORMUM);
    deletePtrArray(_OMEGAVIB);

    deletePtrArray(_CP);
    deletePtrArray(_CV);

    PhysicalChemicalLibrary::unsetup();
  }
}

//////////////////////////////////////////////////////////////////////////////

//thermal conductivity by direct method
CFdouble NitrogenNASALibrary::lambdaD(CFdouble& temperature,
				   CFdouble& pressure)
{
  return 0;
}

//////////////////////////////////////////////////////////////////////////////

//dynamic viscosity by direct method
CFdouble NitrogenNASALibrary::etaD(CFdouble& temperature,
				CFdouble& pressure,
				CFreal* tVec)
{
  return 0;
}

//////////////////////////////////////////////////////////////////////////////

//thermal conductivity by conjugate gradient method
CFdouble NitrogenNASALibrary::lambdaCG(CFdouble& temperature,
                                   CFdouble& pressure)
{
  return 0;  
}

//////////////////////////////////////////////////////////////////////////////

CFdouble NitrogenNASALibrary::lambdaNEQ(CFdouble& temperature,
				     CFdouble& pressure)
{
  return 0;
}

//////////////////////////////////////////////////////////////////////////////

void NitrogenNASALibrary::lambdaVibNEQ(CFreal& temperature,
            RealVector& tVec,
            CFdouble& pressure,
            CFreal& lambdaTrRo,
            RealVector& lambdaInt)
{
  
}

//////////////////////////////////////////////////////////////////////////////

CFdouble NitrogenNASALibrary::etaCG(CFdouble& temperature,
				 CFdouble& pressure,
				 CFreal* tVec)
{
 return 0;  
}

//////////////////////////////////////////////////////////////////////////////

CFdouble NitrogenNASALibrary::sigma(CFdouble& temp,   //electrical conductivity
         CFdouble& pressure,
         CFreal* tVec)
{
 return 0;
}

//////////////////////////////////////////////////////////////////////////////

void NitrogenNASALibrary::gammaAndSoundSpeed(CFdouble& temp,
           CFdouble& pressure,
           CFdouble& rho,
           CFdouble& gamma,
           CFdouble& soundSpeed)
{
    
}

//////////////////////////////////////////////////////////////////////////////

void NitrogenNASALibrary::frozenGammaAndSoundSpeed(CFdouble& temp,
						CFdouble& pressure,
						CFdouble& rho,
						CFdouble& gamma,
						CFdouble& soundSpeed,
						RealVector* tVec)
{
  
  // frozen specifice heat ratio (gamma)
  FORTRAN_NAME(frozengamman2)(_Y, &temp, &gamma);
 
  // frozen speed of sound 
  soundSpeed = std::sqrt(gamma*pressure/rho);

}

//////////////////////////////////////////////////////////////////////////////

CFdouble NitrogenNASALibrary::soundSpeed(CFdouble& temp,
                                      CFdouble& pressure)
{
 
  CFdouble gamma = 0.0;
  CFdouble rho = 0.0; 

  // frozen specifice heat ratio (gamma)
  FORTRAN_NAME(frozengamman2)(_Y, &temp, &gamma);

  // density
  CFdouble mm = 0.0;  
  for (CFint is = 0; is < _NS; ++is) {
       mm += _Y[is]/_MOLARMASSP[is];
  }
  rho = pressure/(_Rgas*temp*mm);  

  return std::sqrt(gamma*pressure/rho); 
}

//////////////////////////////////////////////////////////////////////////////

void NitrogenNASALibrary::setComposition(CFdouble& temp,
             CFdouble& pressure,
             RealVector* x)
{
  // initialization of the molar fractions
  for(CFint i = 0; i < _NS; ++i) {
    _Xini[i] = 1.0;
  }

  // equilibrium composition in terms of molar fractions 
  FORTRAN_NAME(eqcompn2)(&temp, &pressure, &_Xlim, _X);

  if (x != CFNULL) {
    for(CFint i = 0; i < _NS; ++i) {
      (*x)[i] = static_cast<CFreal>(_X[i]);
    }
  }

}

//////////////////////////////////////////////////////////////////////////////

void NitrogenNASALibrary::setDensityEnthalpyEnergy(CFdouble& temp,
            CFdouble& pressure,
            RealVector& dhe)
{
  dhe = 0.0;
  FORTRAN_NAME(n2enthalpy) (&temp, &temp, _HTRANSUM, _HROTUM, _HVIBRUM, _HFORMUM, _HTOTALUM);

  CFreal rho = 0.0;
  for (CFuint i = 0; i < (CFuint)_NS; ++i) { 
      dhe[1] += _Y[i]*_HTOTALUM[i];
      rho += _X[i]*_MOLARMASSP[i]*pressure/(_Rgas*temp);
  }

  dhe[0] = rho;
  dhe[2] = dhe[1] - pressure/rho; 
  
}

//////////////////////////////////////////////////////////////////////////////

void NitrogenNASALibrary::setDensityEnthalpyEnergy(CFdouble& temp,
						RealVector& tVec,
						CFdouble& pressure,
						RealVector& dhe,
						bool storeExtraData)
{
  dhe = 0.0;
  CFdouble tv = 0.0;
  
  // vibrational temperature
  if (_nbTvib == 1) {
      tv = tVec[_nbTvib - 1];
  } else {
      tv = temp;
  } 

  FORTRAN_NAME(n2enthalpy) (&temp, &tv, _HTRANSUM, _HROTUM, _HVIBRUM, _HFORMUM, _HTOTALUM);

  CFreal rho = 0.0;
  for (CFint i = 0; i < _NS; ++i) { 
      dhe[1] += _Y[i]*_HTOTALUM[i];
      rho += _X[i]*_MOLARMASSP[i]*pressure/(_Rgas*temp);
  }

  dhe[0] = rho;
  dhe[2] = dhe[1] - pressure/rho;   	

  // additional component for Park multi-temperature model (N2 in the only molecule)
  if (_nbTvib == 1) {
      const CFuint idxN2 = 1;
      dhe[2+idxN2] = _Y[idxN2]*_HVIBRUM[idxN2];
  } 

  // extra data
  if (storeExtraData) {
     _extraData.dEdT = 0.0;
     _extraData.dHdT = 0.0;

     // frozen specific heats       
     FORTRAN_NAME(n2cv)(&temp, _CV);
     FORTRAN_NAME(n2cp)(&temp, _CP);

     for (CFint i = 0;i < _NS; ++i) {
       
        _extraData.energyTr[i] = _HTRANSUM[i] + _HROTUM[i] + _HFORMUM[i] - (_Rgas/_MOLARMASSP[i])*temp; 
        _extraData.enthalpyTt[i] = _HTOTALUM[i];
        _extraData.dRhoEdRhoi[i] = _extraData.energyTr[i];
        _extraData.energyVib[i] = _HVIBRUM[i];
        _extraData.eElec[i] = 0.0;
       
        _extraData.dHdT += _Y[i]*_CV[i];
	_extraData.dEdT += _Y[i]*_CP[i];
    
     }
  }

}

//////////////////////////////////////////////////////////////////////////////

CFdouble NitrogenNASALibrary::density(CFdouble& temp,
				   CFdouble& pressure,
				   CFreal* tVec)
{
  CFdouble ND = 0.0;
  CFdouble rho = 0.0;

  FORTRAN_NAME(numberd)(&pressure, &temp, &ND);
  FORTRAN_NAME(density)(_X, &ND, &rho);

  return rho;  
}

//////////////////////////////////////////////////////////////////////////////

CFdouble NitrogenNASALibrary::pressure(CFdouble& rho,
            CFdouble& temp,
            CFreal* tVec)
{
  CFdouble p = 0.0;
  FORTRAN_NAME(pressure)(_Y, &rho, &temp, &p); 

  return p; 
}

//////////////////////////////////////////////////////////////////////////////

CFdouble NitrogenNASALibrary::electronPressure(CFreal rhoE,
              CFreal tempE)
{
  return 0;
}

//////////////////////////////////////////////////////////////////////////////

CFdouble NitrogenNASALibrary::energy(CFdouble& temp,
          CFdouble& pressure)
{
  CFdouble ND = 0.0;
  CFdouble rho = 0.0;

  FORTRAN_NAME(numberd)(&pressure, &temp, &ND);
  FORTRAN_NAME(density)(_X, &ND, &rho);

  FORTRAN_NAME(n2enthalpy) (&temp, &temp, _HTRANSUM, _HROTUM, _HVIBRUM, _HFORMUM, _HTOTALUM);

  const CFreal pOvRho = pressure/rho;
  CFdouble intEnergy = 0.;
  for(CFint i = 0; i < _NS; ++i) {
    intEnergy += _Y[i]*(_HTOTALUM[i] - pOvRho);
  }

  return intEnergy;    
}

//////////////////////////////////////////////////////////////////////////////

CFdouble NitrogenNASALibrary::enthalpy(CFdouble& temp,
                                   CFdouble& pressure)
{
  FORTRAN_NAME(n2enthalpy) (&temp, &temp, _HTRANSUM, _HROTUM, _HVIBRUM, _HFORMUM, _HTOTALUM);

  CFdouble ND = 0.0;
  CFdouble rho = 0.0;

  // store the density
  FORTRAN_NAME(numberd)(&pressure, &temp, &ND);
  FORTRAN_NAME(density)(_X, &ND, &rho);

  CFdouble h = 0.0;
  for(CFint i = 0; i < _NS; ++i) {
    h += _Y[i]*_HTOTALUM[i];
  }

  return h;
}

//////////////////////////////////////////////////////////////////////////////

void NitrogenNASALibrary::setElemFractions(const RealVector& yn)
{
  
}

//////////////////////////////////////////////////////////////////////////////

void NitrogenNASALibrary::setElementXFromSpeciesY(const RealVector& ys)
{
  
}

//////////////////////////////////////////////////////////////////////////////

void NitrogenNASALibrary::setElectronFraction(RealVector& ys)
{
  
}

//////////////////////////////////////////////////////////////////////

void NitrogenNASALibrary::setSpeciesFractions(const RealVector& ys)
{
   for (CFint is = 0; is < _NS; ++is) {
   _Y[is] = ys[is];

    if (_Y[is] < 0.0) _Y[is] = 0.0;
    cf_assert(_Y[is] < 1.1);
  }

  // Fills X according to Y
  FORTRAN_NAME(specmasstomolfrac)(_Y, _X);
  
}

//////////////////////////////////////////////////////////////////////////////
void NitrogenNASALibrary::setSpeciesMolarFractions(const RealVector& xs)
{
  std::cout << "Function setSpeciesMolarFractions Not implemeted yet in Nitrogen NASA library" << endl;
  abort();
}

//////////////////////////////////////////////////////////////////////////////

void NitrogenNASALibrary::getSpeciesMolarFractions
(const RealVector& ys, RealVector& xs)
{
  CFreal massTot = 0.;
  for (CFint is = 0; is < _NS; ++is) {
    cf_always_assert (ys[is] > 1.1);

    const CFreal mm = ys[is]/_MOLARMASSP[is];
    massTot += mm;
    xs[is] = mm;
  }
  xs *= 1./massTot; 
}

//////////////////////////////////////////////////////////////////////////////

void NitrogenNASALibrary::getSpeciesMassFractions
(const RealVector& xs, RealVector& ys)
{
  CFreal massTot = 0.;

  for (CFint is = 0; is < _NS; ++is) {
    cf_assert (_X[is] < 1.1);
    const CFreal mm = _X[is]*_MOLARMASSP[is];
    massTot += mm;
    ys[is] = mm;
  }
  ys /= massTot;
}

//////////////////////////////////////////////////////////////////////////////

void NitrogenNASALibrary::getSpeciesMassFractions(RealVector& ys)
{
  for (CFint is = 0; is < _NS; ++is) {
    ys[is] = _Y[is];
  } 
}

//////////////////////////////////////////////////////////////////////////////

void NitrogenNASALibrary::getTransportCoefs(CFdouble& temp,
          CFdouble& pressure,
          CFdouble& lambda,
          CFdouble& lambdacor,
          RealVector& lambdael,
          RealMatrix& eldifcoef,
          RealVector& eltdifcoef)
{
  
}

//////////////////////////////////////////////////////////////////////////////

void NitrogenNASALibrary::getMassProductionTerm(CFdouble& temperature,
              RealVector& tVec,
              CFdouble& pressure,
              CFdouble& rho,
              const RealVector& ys,
              bool flagJac,
              RealVector& omega,
              RealMatrix& jacobian)
{

  for (CFint is = 0; is < _NS; ++is) {
      _Y[is] = ys[is];

      if (_Y[is] < 0.0) _Y[is] = 0.0;
      cf_assert(_Y[is] < 1.1);
      assert(_Y[is] <= 1.0);
    }

   // computation of mass productions terms
   CFdouble temp = temperature;
 
   // vibrational temperature
   CFdouble tv = 0.0;
   if (_nbTvib == 1) {
      tv = tVec[_nbTvib - 1];
   } else {
      tv = temp;
   } 
 
   FORTRAN_NAME(omegan2)(&rho, _Y, &temp, &tv, _OMEGA, _OMEGAVIB);

   // Returning the mass production terms
    for (CFint is = 0; is < _NS; ++is) {
      omega[is] = _OMEGA[is];
    }
  
}

//////////////////////////////////////////////////////////////////////////////

void NitrogenNASALibrary::getSource(CFdouble& temperature,
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
    
   for (CFint is = 0; is < _NS; ++is) {
      _Y[is] = ys[is];

      if (_Y[is] < 0.0) _Y[is] = 0.0;
      cf_assert(_Y[is] < 1.1);
      assert(_Y[is] <= 1.0);
    }

   // computation of mass productions terms
   CFdouble temp = temperature;
 
   // vibrational temperature
   CFdouble tv = 0.0;
   if (_nbTvib == 1) {
      tv = tVec[_nbTvib - 1];
   } else {
      tv = temp;
   } 

   // computation of mass productions terms
   FORTRAN_NAME(omegan2)(&rho, _Y, &temp, &tv, _OMEGA, _OMEGAVIB);

   // Returning the mass production terms
    for (CFint is = 0; is < _NS; ++is) {
      omega[is] = _OMEGA[is];
    }
  
   // source term for the vibrational energy conservation equation
   if (_nbTvib == 1) {
       assert(omegav.size() >= 1);
       omegav[_nbTvib - 1] = _OMEGAVIB[_nbTvib - 1];
   }

}
 
//////////////////////////////////////////////////////////////////////////////

void NitrogenNASALibrary::getRhoUdiff(CFdouble& temperature,
				      CFdouble& pressure,
				      RealVector& normConcGradients,
				      RealVector& normTempGradients,
				      CFreal* tVec,
				      RealVector& rhoUdiff,
				      bool fast)
{
  
}

//////////////////////////////////////////////////////////////////////////////

void NitrogenNASALibrary::getSpeciesTotEnthalpies(CFdouble& temp,
					       RealVector& tVec,
					       CFdouble& pressure,
					       RealVector& hsTot,
					       RealVector* hsVib,
					       RealVector* hsEl)
{
  // vibrational temperature
  CFdouble tv = 0.0;
  
  // vibrational temperature
  if (_nbTvib == 1) {
      tv = tVec[_nbTvib - 1];
  } else {
      tv = temp;
  } 
 
  FORTRAN_NAME(n2enthalpy) (&temp, &tv, _HTRANSUM, _HROTUM, _HVIBRUM, _HFORMUM, _HTOTALUM);

  // returning the total enthalpies per unit mass of species
  for(CFint i = 0; i < _NS; ++i) {
    hsTot[i] = _HTOTALUM[i];
  }
        
}

//////////////////////////////////////////////////////////////////////////////

void NitrogenNASALibrary::getSourceTermVT(CFdouble& temperature,
              RealVector& tVec,
              CFdouble& pressure,
              CFdouble& rho,
              RealVector& omegav,
              CFdouble& omegaRad)
{

}

//////////////////////////////////////////////////////////////////////////////

void NitrogenNASALibrary::getDij_fick(RealVector& dx,
				      CFdouble& pressure,
				      CFdouble& temperature,
				      RealMatrix& Dij,
				      RealVector& rhoUdiff)
{
  
}

//////////////////////////////////////////////////////////////////////////////

void NitrogenNASALibrary::getMolarMasses(RealVector& mm)
{
  assert(mm.size() == static_cast<CFuint>(_NS));

  // check the units
  for (CFint i = 0; i < _NS; ++i) {
    mm[i] = _MOLARMASSP[i];
  }  
}

//////////////////////////////////////////////////////////////////////////////

void NitrogenNASALibrary::testAndWriteProperties()
{
  
}

//////////////////////////////////////////////////////////////////////////////

void NitrogenNASALibrary::testAndWriteSourceTerms()
{
  
}

//////////////////////////////////////////////////////////////////////////////
      
void NitrogenNASALibrary::transportCoeffNEQ(CFreal& temperature, 
					    CFdouble& pressure,
					    CFreal* tVec, 
					    RealVector& normConcGradients,
					    RealVector& normTempGradients,
					    CFreal& eta,
					    CFreal& lambdaTrRo, 
					    RealVector& lambdaInt,
					    RealVector& rhoUdiff)
{
  
}
        
//////////////////////////////////////////////////////////////////////////////
   
void NitrogenNASALibrary::getSourceEE(CFdouble& temperature,
				   RealVector& tVec,
				   CFdouble& pressure,
				   CFdouble& rho,
				   const RealVector& ys,
				   bool flagJac,
				   CFdouble& omegaEE)
{
  
}
      
//////////////////////////////////////////////////////////////////////////////
      
} // namespace NitrogenNASA

} // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

