#include "MutationI/MutationLibrary.hh"
#include "MutationI/Mutation.hh"

#include "Common/CFLog.hh"
#include "Environment/ObjectProvider.hh"
#include "Environment/CFEnv.hh"
#include "Common/OSystem.hh"
#include "Common/BadValueException.hh"
#include "Environment/DirPaths.hh"
#include "Common/Stopwatch.hh"
#include "Common/PEFunctions.hh"
#include "Framework/MeshData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace Mutation {

//////////////////////////////////////////////////////////////////////////////

extern "C" {

  /**
   * Importing of the FORTRAN subroutine LENGTH
   */
  void FORTRAN_NAME(length)(FINT, FINT, FINT, FINT, FINT,
			    FINT, FINT, FINT, FINT, FINT);

  /**
   * Importing of the FORTRAN subroutine INITIALIZE
   */
  void FORTRAN_NAME(initialize)(FDOUBLE, FINT, FDOUBLE, FINT,
				FINT, FINT, FCHAR, FINT, FINT);

  /**
   * Importing of the FORTRAN subroutine NUCLEAR
   */
  void FORTRAN_NAME(nuclear)(FDOUBLE,FINT,FDOUBLE);

  /**
   * Importing of the FORTRAN subroutine EQUIGAMMA
   */
  void FORTRAN_NAME(equigamma)(FDOUBLE, FINT, FINT, FINT, FDOUBLE,
                               FDOUBLE, FDOUBLE, FDOUBLE, FDOUBLE,
                               FDOUBLE, FDOUBLE, FDOUBLE);

  /**
   * Importing of the FORTRAN subroutine FROZENGAMMA
   */
  void FORTRAN_NAME(frozengamma)(FDOUBLE, FINT, FINT, FINT, FDOUBLE,
				 FDOUBLE, FDOUBLE, FDOUBLE, FDOUBLE);

  /**
   * Importing of the FORTRAN subroutine FROZENGAMMANEQ
   */
  void FORTRAN_NAME(frozengammaneq)(FDOUBLE, FINT, FINT, FINT, FINT,
				    FDOUBLE, FDOUBLE, FDOUBLE, FDOUBLE,
				    FDOUBLE, FDOUBLE);
  /**
   * Importing of the FORTRAN subroutine FROZENGAMMA
   */
  void FORTRAN_NAME(frozengammafast)(FDOUBLE, FINT, FINT, FINT, FDOUBLE,
				 FDOUBLE, FDOUBLE, FDOUBLE, FDOUBLE);

  /**
   * Importing of the FORTRAN subroutine FROZENGAMMANEQ
   */
  void FORTRAN_NAME(frozengammaneqfast)(FDOUBLE, FINT, FINT, FINT, FINT,
				    FDOUBLE, FDOUBLE, FDOUBLE, FDOUBLE,
				    FDOUBLE, FDOUBLE, FDOUBLE);

  /**
   * Importing of the FORTRAN subroutine COMPOSITION
   */
  void FORTRAN_NAME(composition)(FDOUBLE, FINT, FINT, FINT, FDOUBLE,
                                 FDOUBLE, FDOUBLE, FDOUBLE, FDOUBLE);
  /**
   * Importing of the FORTRAN subroutine COMPFRACT
   */
  void FORTRAN_NAME(compfract)(FDOUBLE, FINT, FDOUBLE);

  /**
   * Importing of the FORTRAN subroutine NUMBERD
   */
  void FORTRAN_NAME(numberd)(FDOUBLE, FINT, FDOUBLE, FDOUBLE,
                             FDOUBLE, FDOUBLE, FDOUBLE);

  /**
   * Importing of the FORTRAN subroutine DENSITY
   */
  void FORTRAN_NAME(density)(FDOUBLE, FINT, FDOUBLE, FDOUBLE, FDOUBLE);

  /**
   * Importing of the FORTRAN subroutine MOLARMASS
   */
  void FORTRAN_NAME(molarmass)(FDOUBLE ,FINT, FDOUBLE, FDOUBLE, FDOUBLE);

  /**
   * Importing of the FORTRAN subroutine MASSCOMPOSITION
   */
  void FORTRAN_NAME(masscomposition)(FDOUBLE, FINT, FINT,FINT, FDOUBLE,
                                     FDOUBLE, FDOUBLE, FDOUBLE, FDOUBLE);

  /**
   * Importing of the FORTRAN subroutine ENTHALPYMASS
   */
  void FORTRAN_NAME(enthalpymass)(FDOUBLE, FINT, FINT ,FINT, FDOUBLE,
                                  FDOUBLE, FDOUBLE, FDOUBLE, FDOUBLE,
                                  FDOUBLE, FDOUBLE, FDOUBLE, FDOUBLE,
                                  FDOUBLE, FDOUBLE);

  /**
   * Importing of the FORTRAN subroutine ENTHALPYFORM
   */
  void FORTRAN_NAME(enthalpyform)(FDOUBLE, FINT, FDOUBLE);

  /**
   * Importing of the FORTRAN subroutine ENERGY
   */
  void FORTRAN_NAME(energy)(FDOUBLE, FINT, FINT, FINT, FDOUBLE,
                            FDOUBLE, FDOUBLE, FDOUBLE, FDOUBLE,
                            FDOUBLE, FDOUBLE, FDOUBLE, FDOUBLE,
                            FDOUBLE, FDOUBLE);
  /**
   * Importing of the FORTRAN subroutine ENERGYVIB
   */
  void FORTRAN_NAME(energyvib)(FDOUBLE, FINT, FINT, FINT, FINT, FDOUBLE,
                            FDOUBLE, FDOUBLE, FDOUBLE, FDOUBLE,
                            FDOUBLE, FDOUBLE, FDOUBLE, FDOUBLE,
                            FDOUBLE, FDOUBLE);
  /**
   * Importing of the FORTRAN subroutine ENERGYVIBIONS
   */
  void FORTRAN_NAME(energyvibions)(FDOUBLE, FINT, FINT, FINT, FINT, FDOUBLE,
                            FDOUBLE, FDOUBLE, FDOUBLE, FDOUBLE,
                            FDOUBLE, FDOUBLE, FDOUBLE, FDOUBLE,
                            FDOUBLE, FDOUBLE);

  /**
   * Importing of the FORTRAN subroutine ENTHALPY
   */
  void FORTRAN_NAME(enthalpy)(FDOUBLE, FINT, FINT, FINT, FDOUBLE,
                              FDOUBLE, FDOUBLE, FDOUBLE, FDOUBLE,
                              FDOUBLE, FDOUBLE, FDOUBLE, FDOUBLE,
                              FDOUBLE, FDOUBLE);
  /**
   * Importing of the FORTRAN subroutine ENTHALPYVIB
   */
  void FORTRAN_NAME(enthalpyvib)(FDOUBLE, FINT, FINT, FINT, FINT, FDOUBLE,
                              FDOUBLE, FDOUBLE, FDOUBLE, FDOUBLE,
                              FDOUBLE, FDOUBLE, FDOUBLE, FDOUBLE,
                              FDOUBLE, FDOUBLE);
  /**
   * Importing of the FORTRAN subroutine ENTHALPYVIBIONS
   */
  void FORTRAN_NAME(enthalpyvibions)(FDOUBLE, FINT, FINT, FINT, FINT, FDOUBLE,
                              FDOUBLE, FDOUBLE, FDOUBLE, FDOUBLE,
                              FDOUBLE, FDOUBLE, FDOUBLE, FDOUBLE,
                              FDOUBLE, FDOUBLE);

  /**
   * Importing of the FORTRAN subroutine COMPOTOL
   */
  void FORTRAN_NAME(compotol)(FDOUBLE, FDOUBLE, FDOUBLE);

  /**
   * Importing of the FORTRAN subroutine COMPOTOL2
   */
  void FORTRAN_NAME(compotol2)(FDOUBLE, FDOUBLE, FDOUBLE);

  /**
   * Importing of the FORTRAN subroutine COLLISION
   */
  void FORTRAN_NAME(collision)(FDOUBLE, FINT, FDOUBLE, FINT,
                               FDOUBLE, FDOUBLE, FDOUBLE, FDOUBLE);

  /**
   * Importing of the FORTRAN subroutine LAMBDAINT
   */
  void FORTRAN_NAME(lambdaint)(FDOUBLE, FINT, FDOUBLE, FINT,
                               FDOUBLE, FDOUBLE, FDOUBLE);


 /**
   * Importing of the FORTRAN subroutine LAMBDAVIB
   */
  void FORTRAN_NAME(lambdavib)(FDOUBLE, FINT, FDOUBLE, FINT,
                               FDOUBLE, FDOUBLE, FDOUBLE);

/**
   * Importing of the FORTRAN subroutine LAMBDAVIBIONS
   */
  void FORTRAN_NAME(lambdavibions)(FDOUBLE, FINT, FDOUBLE, FINT,
                               FDOUBLE, FDOUBLE, FDOUBLE);

  /**
   * Importing of the FORTRAN subroutine LAMBDACHID
   */
  void FORTRAN_NAME(lambdachid)(FDOUBLE, FINT, FDOUBLE, FINT,
                               FDOUBLE, FDOUBLE, FDOUBLE);

  /**
   * Importing of the FORTRAN subroutine LAMBDAREAD
   */
  void FORTRAN_NAME(lambdaread)(FINT, FINT, FDOUBLE, FINT, FDOUBLE, FINT,
                               FDOUBLE, FDOUBLE, FDOUBLE, FDOUBLE);

  /**
   * Importing of the FORTRAN subroutine LAMBDACHICG
   */
  void FORTRAN_NAME(lambdachicg)(FDOUBLE, FINT, FDOUBLE, FINT,
                               FDOUBLE, FDOUBLE, FDOUBLE);

  /**
   * Importing of the FORTRAN subroutine LAMBDAREACG
   */
  void FORTRAN_NAME(lambdareacg)(FINT, FINT, FDOUBLE, FINT, FDOUBLE, FINT,
                               FDOUBLE, FDOUBLE, FDOUBLE, FDOUBLE);

  /**
   * Importing of the FORTRAN subroutine LAMBDADAE
   */
  void FORTRAN_NAME(lambdae)(FDOUBLE, FINT, FDOUBLE, FINT,
                               FDOUBLE, FDOUBLE, FDOUBLE, FINT);

  /**
   * Importing of the FORTRAN subroutine LAMBDADAE
   */
  void FORTRAN_NAME(etad)(FDOUBLE, FINT, FDOUBLE, FINT,
                          FDOUBLE, FDOUBLE);

  /**
   * Importing of the FORTRAN subroutine LAMBDADACG
   */
  void FORTRAN_NAME(etacg)(FDOUBLE, FINT, FDOUBLE, FINT,
                          FDOUBLE, FDOUBLE);

   /**
   * Importing of the FORTRAN subroutine SIGMAE
   */
  void FORTRAN_NAME(sigmae)(FDOUBLE, FINT, FDOUBLE, FINT,
                          FDOUBLE, FDOUBLE, FDOUBLE, FINT);

  /**
   * Importing of the FORTRAN subroutine LTEVEFNEU
   */
  void FORTRAN_NAME(ltevefneu)(FINT, FINT, FDOUBLE, FINT, FDOUBLE,
                   FINT, FDOUBLE, FDOUBLE, FDOUBLE, FDOUBLE, FDOUBLE,
                   FDOUBLE, FDOUBLE, FDOUBLE, FDOUBLE, FDOUBLE);

  /**
   * Importing of the FORTRAN subroutine NUCMASSTOMOLFRAC
   */
  void FORTRAN_NAME(nucmasstomolfrac)(FDOUBLE, FINT, FINT, FINT, FDOUBLE, FDOUBLE);

  /**
   * Importing of the FORTRAN subroutine SPECMASSTOMOLFRAC
   */
  void FORTRAN_NAME(specmasstomolfrac)(FDOUBLE, FINT, FINT, FINT, FDOUBLE, FDOUBLE);

  /**
   * Importing of the FORTRAN subroutine SETSPECIESMOLARMASS
   */
  void FORTRAN_NAME(setspeciesmolarmass)(FDOUBLE, FINT, FDOUBLE);

  /**
   * Importing of the FORTRAN subroutine PRESSURE
   */
  void FORTRAN_NAME(pressure)(FDOUBLE, FINT, FDOUBLE, FDOUBLE, FDOUBLE, FDOUBLE, FDOUBLE);

  /**
   * Importing of the FORTRAN subroutine EPRESSURE
   */
  void FORTRAN_NAME(epressure)(FDOUBLE, FINT, FDOUBLE, FDOUBLE,
			       FDOUBLE, FDOUBLE, FDOUBLE, FDOUBLE, FDOUBLE);

  /**
   * Importing of the FORTRAN subroutine CORRECTION
   */
  void FORTRAN_NAME(correction)(FDOUBLE, FINT, FDOUBLE, FINT, FDOUBLE, FINT, FDOUBLE);

  /**
   * Importing of the FORTRAN subroutine ARRHENIUS
   */
  void FORTRAN_NAME(arrhenius)(FDOUBLE, FINT, FDOUBLE, FINT, FINT, FINT,
                               FDOUBLE, FDOUBLE, FDOUBLE, FDOUBLE, FDOUBLE,
                               FDOUBLE, FDOUBLE, FDOUBLE);
  /**
   * Importing of the FORTRAN subroutine ARRHENIUSVT
   */
  void FORTRAN_NAME(arrheniusvt)(FDOUBLE, FINT, FDOUBLE, FINT, FINT, FINT,
                               FDOUBLE, FDOUBLE, FDOUBLE, FDOUBLE, FDOUBLE,
                               FDOUBLE, FDOUBLE);
  /**
   * Importing of the FORTRAN subroutine ARRHENIUSVTIONS
   */
  void FORTRAN_NAME(arrheniusvtions)(FDOUBLE, FINT, FDOUBLE, FINT, FINT, FINT,
                               FDOUBLE, FDOUBLE, FDOUBLE, FDOUBLE, FDOUBLE,
                               FDOUBLE, FDOUBLE);
  /**
   * Importing of the FORTRAN subroutine ARRHENIUSPVTY
   */
  void FORTRAN_NAME(arrheniuspvty)(FDOUBLE, FINT, FDOUBLE, FINT, FINT, FINT,
                                   FDOUBLE, FDOUBLE, FDOUBLE, FDOUBLE, FINT,
                                   FDOUBLE, FDOUBLE, FDOUBLE, FDOUBLE);

  /**
   * Importing of the FORTRAN subroutine VTTRANSFER
   */
  void FORTRAN_NAME(vttransfer)(FINT, FINT, FDOUBLE, FDOUBLE, FDOUBLE, FDOUBLE,
                                FDOUBLE, FDOUBLE, FDOUBLE, FDOUBLE, FINT,
                                FDOUBLE, FDOUBLE, FDOUBLE, FDOUBLE, FDOUBLE);

  /**
   * Importing of the FORTRAN subroutine SMD
   */
  void FORTRAN_NAME(smd)(FDOUBLE, FINT, FDOUBLE, FINT, FINT, FINT,FDOUBLE, FDOUBLE,
			 FDOUBLE, FDOUBLE, FDOUBLE,FDOUBLE, FDOUBLE, FDOUBLE);

  /**
   * Importing of the FORTRAN subroutine SMNEUTCG
   */
  void FORTRAN_NAME(smneutcg)(FDOUBLE, FINT, FDOUBLE, FINT, FDOUBLE, FDOUBLE,
			      FDOUBLE, FDOUBLE, FDOUBLE);

  /**
   * Importing of the FORTRAN subroutine SMNEUTSU
   */
  void FORTRAN_NAME(smneutsut)(FDOUBLE, FINT, FDOUBLE, FINT, FDOUBLE, FDOUBLE,
			       FDOUBLE, FDOUBLE, FDOUBLE);

  /**
   * Importing of the FORTRAN subroutine SMNEUTD
   */
  void FORTRAN_NAME(smneutd)(FDOUBLE, FINT, FDOUBLE, FINT, FDOUBLE, FDOUBLE,
                             FDOUBLE, FDOUBLE, FDOUBLE);

  /**
   * Importing of the FORTRAN subroutine SOURCEVT
   */
  void FORTRAN_NAME(sourcevt) (FDOUBLE, FINT, FDOUBLE, FINT, FINT, FINT, FDOUBLE, FDOUBLE, FDOUBLE,
			       FDOUBLE, FDOUBLE, FINT, FDOUBLE, FDOUBLE, FDOUBLE);
  /**
   * Importing of the FORTRAN subroutine SOURCEVT
   */
  void FORTRAN_NAME(sourcevtn) (FDOUBLE, FINT, FDOUBLE, FINT, FINT, FINT, FDOUBLE, FDOUBLE, FDOUBLE,
			       FDOUBLE, FDOUBLE, FDOUBLE, FDOUBLE, FDOUBLE);
  /**
   * Importing of the FORTRAN subroutine SOURCEVT
   */
  void FORTRAN_NAME(sourcevtions) (FDOUBLE, FINT, FDOUBLE, FINT, FINT, FINT, FDOUBLE, FDOUBLE, FDOUBLE,
			       FDOUBLE, FDOUBLE, FINT, FDOUBLE, FDOUBLE, FDOUBLE);

 /**
   * Importing of the FORTRAN subroutine GETCHARGE
   */
  void FORTRAN_NAME(getcharge)(FINT, FINT, FINT);

  /**
   * Importing of the FORTRAN subroutine RGAS
   */
  void FORTRAN_NAME(rgas)(FDOUBLE, FINT, FDOUBLE);
}

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<MutationLibrary,
			    PhysicalPropertyLibrary,
			    MutationModule,
			    1>
mutationLibraryProvider("Mutation");

//////////////////////////////////////////////////////////////////////////////

void MutationLibrary::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFdouble >("deltaT","Delta temperature.");
  options.addConfigOption< bool >("useLookUpTable","Flag telling if to use the look up tables.");
  options.addConfigOption< bool >
    ("noElectronicEnergy","Flag telling to ignore the electronic energy.");
  options.addConfigOption< CFdouble >("Pmax","Maximum pressure in the table.");
  options.addConfigOption< std::vector<std::string> >("lookUpVars","Name of the vars to store in the table.");
  options.addConfigOption< CFdouble >("Tmax","Maximum temperature in the table.");
  options.addConfigOption< CFdouble >("Pmin","Minimum pressure in the table.");
  options.addConfigOption< int >("sonine","Sonine polynome order.");
  options.addConfigOption< int >("IMOD","Compute transport properties.");
  options.addConfigOption< CFdouble >("deltaP","Delta pressure.");
  options.addConfigOption< std::string >("mixtureName","Name of the mixture.");
  options.addConfigOption< std::string >("reactionName","Name of the reaction.");
  options.addConfigOption< CFdouble >("Tmin","Minimum temperature in the table.");
  options.addConfigOption< CFdouble >("Xlim","Species molar fractions limit.");
  options.addConfigOption< std::string >("thermCondAlgo","Algo to compute thermal conductivity.");
  options.addConfigOption< std::string >("dynViscAlgo","Algo to compute dynamic viscosity.");
  options.addConfigOption< CFdouble >("TminFix","Minimum allowed temperature for chemistry.");
  options.addConfigOption< CFdouble,
    Config::DynamicOption<> >("FactorOmega","Factor to reduce stiffness of chemical sorce terms.");
  options.addConfigOption< CFdouble >("GammaN","Factor of Catalycity of N.");
  options.addConfigOption< CFdouble >("GammaO","Factor of Catalycity of O.");
}

//////////////////////////////////////////////////////////////////////////////

MutationLibrary::MutationLibrary(const std::string& name)
  : PhysicalChemicalLibrary(name),
    _nameToIdxVar(),
    _lookUpTables()
{
  addConfigOptionsTo(this);
  
  _mixtureName = "";
  setParameter("mixtureName",&_mixtureName);
  
  //"empty" for LTE, no finite rate chemistry
  _reactionName = "empty";
  setParameter("reactionName",&_reactionName);

  _useLookUpTable = false;
  setParameter("useLookUpTable",&_useLookUpTable);

  _noElectEnergy = false;
  setParameter("noElectronicEnergy",&_noElectEnergy);

  _lkpVarNames = vector<std::string>();
  setParameter("lookUpVars",&_lkpVarNames);

  _Tmin = 100.0;
  setParameter("Tmin",&_Tmin);

  _Tmax = 2000.0;
  setParameter("Tmax",&_Tmax);

  _deltaT = 10.0;
  setParameter("deltaT",&_deltaT);
  
  _GammaN = 0.0;
  setParameter("GammaN",&_GammaN);
  
  _GammaO = 0.0;
  setParameter("GammaO",&_GammaO);

  _pmin = 10000.0;
  setParameter("Pmin",&_pmin);

  _pmax = 1000000.0;
  setParameter("Pmax",&_pmax);

  _deltaP = 1000.0;
  setParameter("deltaP",&_deltaP);
  // IT DOESN t work with SONINE =2 The problem is somewhere inside
  // CORRECTION
  // _sonine = 2;
  _sonine = 1;
  setParameter("sonine",&_sonine);
  // 1: zero order, 2: first order in Mutation
  cf_assert(_sonine < 3);

  // = 0: default, = 1: if transport properties will be computed when using
  // the library (you can use mutation without specifying the transport
  // collision integrals if not needed)
  _imod = 1;
  setParameter("IMOD",&_imod);

  _Xlim = 1.0e-15;
  setParameter("Xlim",&_Xlim);

  _lambdaAlgoStr = "CG";
  setParameter("thermCondAlgo",&_lambdaAlgoStr);

  _etaAlgoStr = "CG";
  setParameter("dynViscAlgo",&_etaAlgoStr);

  _TminFix = 800.;
  setParameter("TminFix",&_TminFix);
  
  _factorOmega = 1.0;
  setParameter("FactorOmega",&_factorOmega);
  
  EPS = 1e-6; // this value gives problems with air11 (composition, speed of sound) and LTE
  TOL = 1e-9; // recommended to replace with the use of Xlim and COMPOTOL2
}

//////////////////////////////////////////////////////////////////////////////

MutationLibrary::~MutationLibrary()
{
}

//////////////////////////////////////////////////////////////////////////////

void MutationLibrary::configure ( Config::ConfigArgs& args )
{
  PhysicalChemicalLibrary::configure(args);

  CFLog(NOTICE, "Tmax   = " << _Tmax << "\n");
  CFLog(NOTICE, "Tmin   = " << _Tmin << "\n");
  CFLog(NOTICE, "deltaT = " << _deltaT << "\n");
  CFLog(NOTICE, "Pmax   = " << _pmax << "\n");
  CFLog(NOTICE, "Pmin   = " << _pmin << "\n");
  CFLog(NOTICE, "deltaP = " << _deltaP << "\n");
  CFLog(NOTICE, "GammaN = " << _GammaN << "\n");
  CFLog(NOTICE, "GammaO = " << _GammaO << "\n");

  if (_lambdaAlgoStr == "CG") _lambdaAlgo = LAMBDACG;
  if (_lambdaAlgoStr == "Direct") _lambdaAlgo = LAMBDAD;
  if (_etaAlgoStr == "CG") _etaAlgo = ETACG;
  if (_etaAlgoStr == "Direct") _etaAlgo = ETAD;
}

//////////////////////////////////////////////////////////////////////////////

void MutationLibrary::setup()
{
  // if this is a parallel simulation, only ONE process at a time
  // sets the library
  setLibrarySequentially();
}
      
//////////////////////////////////////////////////////////////////////////////

void MutationLibrary:: setLibrarySequentially()
{
  // if nobody has configured the library path, it is done here
  // with default = basedir + "plugins/MutationI/data/mutation/"
  if (m_libPath == "") {
    CFLog(NOTICE, "MutationLibrary::libpath set to default" << "\n");
    std::string baseDir = Environment::DirPaths::getInstance().getBaseDir().string();
    m_libPath = baseDir + "/plugins/MutationI/data/mutation/";
  }

  ofstream fout("mutation.in");
  fout << _mixtureName << endl;
  fout << _reactionName << endl;

  std::string command1 = "cp -R " + m_libPath + "data .";
  Common::OSystem::getInstance().executeCommand(command1);

  CFout << "MutationLibrary::libpath     = "  << m_libPath << "\n";
  CFout << "MutationLibrary::mixtureName = "  << _mixtureName << "\n";
  CFout << "MutationLibrary::reactionName = " << _reactionName << "\n";

  try {
    if (_mixtureName == "") {
      throw Common::BadValueException (FromHere(), "Mixture name not set!!!");
    }
  }
  catch (Common::BadValueException& ex) {
    CFout << ex.what() << "\n";
    CFout << "Aborting ..." << "\n";
    abort();
  }

  // calculate the needed size for the work vectors
  FORTRAN_NAME(length)(&LWR1,&LWR2,&LWR3,&LWI,&LWC,&_NS,&NE,&_NC,&NREA,&IBINIJ);

  _atomicityCoeff.resize(_NS);

  // allocate the private data
  WR1 = new CFdouble[LWR1];
  WR2 = new CFdouble[LWR2];
  WR3 = new CFdouble[LWR3];
  WI  = new int[LWI];
  WC  = new char[LWC*4+1];
  FIJ = new CFdouble[_NS*(_NS-1)/2];

  Xn   =  new CFdouble[_NC];
  Yn   =  new CFdouble[_NC];
  Y   =  new CFdouble[_NS];
  X    =  new CFdouble[_NS];
  XTOL =  new CFdouble[_NS];
  Xini =  new CFdouble[_NS];
  Yini =  new CFdouble[_NS];
  LAMBDAEL =  new CFdouble[_NC];
  ELDIFCOEF = new CFdouble[_NC*_NC];
  ELTDIFCOEF = new CFdouble[_NC];
  OMEGA = new CFdouble[_NS];
  DWDP = new CFdouble[_NS];
  DWDT = new CFdouble[_NS];
  DWDYI = new CFdouble[_NS*_NS];
  MOLARMASSP = new CFdouble[_NS];
  DF =  new CFdouble[_NS];
  JDIF =  new CFdouble[_NS];
  DF =  new CFdouble[_NS];
  JDIF =  new CFdouble[_NS];

  for(int i = 0; i < _NS; ++i) {
    X[i] = 1.0;
    Y[i] = 1.0;
  }

  _HTOTAL =  new CFdouble[_NS];
  _HTRANS =  new CFdouble[_NS];
  _HELECT =  new CFdouble[_NS];
  _HROT   =  new CFdouble[_NS];
  _HVIBR  =  new CFdouble[_NS];
  _HFORM  =  new CFdouble[_NS];

  _HTOTALP =  new CFdouble[_NS];
  _HTRANSP =  new CFdouble[_NS];
  _HELECTP =  new CFdouble[_NS];
  _HROTP   =  new CFdouble[_NS];
  _HVIBRP  =  new CFdouble[_NS];
  _HFORMP  =  new CFdouble[_NS];

  _HMTOTAL =  new CFdouble[_NS];
  _HMTRANS =  new CFdouble[_NS];
  _HMELECT =  new CFdouble[_NS];
  _HMROT   =  new CFdouble[_NS];
  _HMVIBR  =  new CFdouble[_NS];
  _HMFORM  =  new CFdouble[_NS];

  _ETOTAL =  new CFdouble[_NS];
  _ETRANS =  new CFdouble[_NS];
  _EELECT =  new CFdouble[_NS];
  _EROT   =  new CFdouble[_NS];
  _EVIBR  =  new CFdouble[_NS];
  _EFORM  =  new CFdouble[_NS];

  /// @TODO remove this once the number of Tvibs is equal in coolfluid and
  /// in the allocated array inside Mutation
  _TVARRAY =  new CFdouble[2];
  _STVIB   =  new CFdouble[2];

  _CPE =  new CFdouble[_NS];
  _CPR =  new CFdouble[_NS];
  _CPV =  new CFdouble[_NS];
  _CPINT = new CFdouble[_NS];

  _CHIH = new CFdouble[_NS];
  _CHARGE = new int[_NS];

  _TVIBEPS1 = new CFdouble[2];
  _TVIBEPS  = new CFdouble[2];

  _LAMBDAVIB = new CFdouble[2];
  _CPVIB = new CFdouble[2];

  // species translational-rotational energies
  _extraData.energyTr.resize(_NS);
  _extraData.energyVib.resize(_NS);
  _extraData.enthalpyForm.resize(_NS);
  _extraData.dRhoEdRhoi.resize(_NS);
  _extraData.dRhoEvdRhoi.resize(_NS);
  _tmp.resize(_NS);

  // initialize some data in work vectors
  // Careful!!! WR2 is filled by function "collision"
  FORTRAN_NAME(initialize)(WR1,&LWR1,WR3,&LWR3,WI,&LWI,WC,&LWC,&_imod);

  // calculate nuclear fraction ratios;
  FORTRAN_NAME(nuclear)(WR1,&LWR1,Xn);

  PhysicalChemicalLibrary::setup();

  if (_lkpVarNames.size() > 0) {
    setLookUpTables();
    _useLookUpTable = true;
  }

  // set the species molar masses
  FORTRAN_NAME(setspeciesmolarmass)(WR1, &LWR1, MOLARMASSP);

  // get the charge of each species
  FORTRAN_NAME(getcharge)(WI, &LWI, _CHARGE);

  // constant of perfect gases
  FORTRAN_NAME(rgas)(WR1, &LWR1, &_Rgas);

  // formation enthalpy per unity mass
  FORTRAN_NAME(enthalpyform)(WR1, &LWR1, _HFORM);
  for (int i = 0; i < _NS; ++i) {
    _extraData.enthalpyForm[i] = _HFORM[i]/MOLARMASSP[i];
  }

  setMoleculesIDs(_moleculesIDs);
  _flagMoleculesIDs.resize(_NS,false);
  for (CFuint i = 0; i < _moleculesIDs.size(); ++i) {
    _flagMoleculesIDs[_moleculesIDs[i]] = true;
  }
  for (int i = 0; i < _NS; ++i) {
    _atomicityCoeff[i] = (_flagMoleculesIDs[i]) ? 2.5 : 1.5;
  }

  // fix this
  _molecule2EqIDs.resize(_nbTvib);
  if ((_NS == 5 || _NS == 11) && _nbTvib == 2) {
    _molecule2EqIDs[0] = 0;
    _molecule2EqIDs[1] = 2;
  }

  if (getNbTe() == 1) {
    _electrEnergyID = _nbTvib;
  }
  
  // set the flag telling if the mixture is ionized
  _hasElectrons = (NE == 1) ? true : false;
 
  std::string command2 = "rm -fr data/ mutation.in";
  Common::OSystem::getInstance().executeCommand(command2);


}

//////////////////////////////////////////////////////////////////////////////

void MutationLibrary::unsetup()
{
  if(isSetup()) {

    deletePtrArray(WR1);
    deletePtrArray(WR2);
    deletePtrArray(WR3);
    deletePtrArray(WC);
    deletePtrArray(WI);
    deletePtrArray(FIJ);

    deletePtrArray(Xn);
    deletePtrArray(Yn);
    deletePtrArray(Y);
    deletePtrArray(X);
    deletePtrArray(XTOL);
    deletePtrArray(Xini);
    deletePtrArray(Yini);
    deletePtrArray(LAMBDAEL);
    deletePtrArray(ELDIFCOEF);
    deletePtrArray(ELTDIFCOEF);
    deletePtrArray(OMEGA);
    deletePtrArray(DWDP);
    deletePtrArray(DWDT);
    deletePtrArray(DWDYI);
    deletePtrArray(MOLARMASSP);
    deletePtrArray(DF);
    deletePtrArray(JDIF);

    deletePtrArray(_HTOTAL);
    deletePtrArray(_HTRANS);
    deletePtrArray(_HELECT);
    deletePtrArray(_HROT);
    deletePtrArray(_HVIBR);
    deletePtrArray(_HFORM);

    deletePtrArray(_HTOTALP);
    deletePtrArray(_HTRANSP);
    deletePtrArray(_HELECTP);
    deletePtrArray(_HROTP);
    deletePtrArray(_HVIBRP);
    deletePtrArray(_HFORMP);

    deletePtrArray(_HMTOTAL);
    deletePtrArray(_HMTRANS);
    deletePtrArray(_HMELECT);
    deletePtrArray(_HMROT);
    deletePtrArray(_HMVIBR);
    deletePtrArray(_HMFORM);

    deletePtrArray(_ETOTAL);
    deletePtrArray(_ETRANS);
    deletePtrArray(_EELECT);
    deletePtrArray(_EROT);
    deletePtrArray(_EVIBR);
    deletePtrArray(_EFORM);

    /// @TODO remove this once the number of Tvibs is equal in coolfluid and
    /// in the allocated array inside Mutation
    deletePtrArray(_TVARRAY);
    deletePtrArray(_STVIB);
    deletePtrArray(_CPE);
    deletePtrArray(_CPR);
    deletePtrArray(_CPV);
    deletePtrArray(_CPINT);

    deletePtrArray(_CHIH);
    deletePtrArray(_CHARGE);
    deletePtrArray( _TVIBEPS1);
    deletePtrArray( _TVIBEPS);
    deletePtrArray(_LAMBDAVIB);
    deletePtrArray(_CPVIB);

    PhysicalChemicalLibrary::unsetup();
  }
}

//////////////////////////////////////////////////////////////////////////////

CFdouble MutationLibrary::lambdaD(CFdouble& temperature, //thermal conductivity by direct method
				  CFdouble& pressure)
{
  //composition must be called before!
  CFdouble lambdaTot = 0.0;
  CFdouble lambdaInth = 0.0;
  CFdouble lambdaRea = 0.0;
  CFdouble lambdaTh = 0.0;
  CFdouble lambdaTe = 0.0;
  CFdouble temp = temperature;
  if (temp < 100.) {
    // cout << " MutationLibrary::lambdaD -> Watch Out !! T <<<<< 100." << temp << endl;
    temp = 100.;
  }

  CFdouble ND = 0.0;
  FORTRAN_NAME(numberd)(WR1,&LWR1,&pressure,&temp,&temp,X,&ND); //ND needed for collision below

   //Thermodynamic properties

  // - Species specific heat per unit mole
  FORTRAN_NAME(enthalpy)(WR1, &LWR1, WI, &LWI, &temp, &temp, &temp, &temp, &pressure,
                         _HTOTAL, _HTRANS, _HELECT, _HROT, _HVIBR, _HFORM);
  CFdouble tempEps1 = temp*(1.0+EPS);
  CFdouble tempEps = temp*EPS;
  FORTRAN_NAME(enthalpy)(WR1, &LWR1, WI, &LWI, &tempEps1, &tempEps1, &tempEps1, &tempEps1,
                         &pressure, _HTOTALP, _HTRANSP, _HELECTP, _HROTP, _HVIBRP, _HFORMP);

  for(int is=0 ; is < _NS; is++) {
    _CPE[is]   = (_HELECTP[is] - _HELECT[is]) / (tempEps);
    _CPR[is]   = (_HROTP[is] - _HROT[is]) / (tempEps);
    _CPV[is]   = (_HVIBRP[is] - _HVIBR[is]) / (tempEps);
    _CPINT[is] = _CPE[is] + _CPR[is] + _CPV[is];
  }

  FORTRAN_NAME(compotol)(X, &TOL, XTOL);
  //Update of kinetic data
  FORTRAN_NAME(collision)(WR1, &LWR1, WR2, &LWR2, &temp, &temp, &ND, X);

  //A. Heavy particles properties

  //Eucken corrections for internal energy (no inelastic collisions)

  FORTRAN_NAME(lambdaint)(WR1, &LWR1, WR2, &LWR2, _CPINT, XTOL, &lambdaInth);
  FORTRAN_NAME(lambdachid)(WR1, &LWR1, WR2, &LWR2, XTOL, &lambdaTh, _CHIH);


  FORTRAN_NAME(lambdaread)(WI, &LWI, WR1, &LWR1, WR2, &LWR2, _HTOTAL, &temp, XTOL, &lambdaRea);

  //B. Electron gas properties
  if(presenceElectron()){
    int i=3;
    FORTRAN_NAME(lambdae)(WR1, &LWR1, WR2, &LWR2, X, &temp, &lambdaTe, &i);
  }

  lambdaTot = lambdaInth + lambdaTh + lambdaRea + lambdaTe;

  return lambdaTot;
}

//////////////////////////////////////////////////////////////////////////////

//dynamic viscosity by direct method
CFdouble MutationLibrary::etaD(CFdouble& temperature,
			       CFdouble& pressure,
			       CFreal* tVec)
{
  //composition must be called before!
  CFdouble eta = 0.0;
  CFdouble temp = temperature;
  if (temp < 100.) {
    //cout << " MutationLibrary::etaD -> Watch Out !! T <<<<< 100." << temp << endl;
    temp = 100.;
  }

  CFdouble Te = getTe(temp,tVec);

  //ND needed for collision below
  CFdouble ND = 0.0;
  FORTRAN_NAME(numberd)(WR1,&LWR1,&pressure,&temp,&Te,X,&ND);

  // FORTRAN_NAME(compotol2)(X, &_Xlim, XTOL);
  FORTRAN_NAME(compotol)(X, &TOL, XTOL);
  //Update of kinetic data
  FORTRAN_NAME(collision)(WR1, &LWR1, WR2, &LWR2, &temp, &temp, &ND, X);

  //A. Heavy particles properties

  //Eucken corrections for internal energy (no inelastic collisions)

   FORTRAN_NAME(etad)(WR1, &LWR1, WR2, &LWR2, XTOL, &eta);

   return eta;
}

//////////////////////////////////////////////////////////////////////////////

//thermal conductivity by conjugate gradient method
CFdouble MutationLibrary::lambdaCG(CFdouble& temperature,
                                   CFdouble& pressure)
{
  //composition must be called before!
  CFdouble lambdaInth = 0.0;
  CFdouble lambdaRea = 0.0;
  CFdouble lambdaTh = 0.0;
  CFdouble lambdaTe = 0.0;
  CFdouble ND = 0.0;

  CFdouble temp = temperature;
  if (temp < 100.) {
    //cout << " MutationLibrary::lambdaCG -> Watch Out !! T <<<<< 100." << temp << endl;
    temp = 100.;
  }

  FORTRAN_NAME(numberd)(WR1,&LWR1,&pressure,&temp,&temp,X,&ND); //ND needed for collision below

   //Thermodynamic properties

  // - Species specific heat per unit mole
  FORTRAN_NAME(enthalpy)(WR1, &LWR1, WI, &LWI, &temp, &temp, &temp, &temp, &pressure,
                         _HTOTAL, _HTRANS, _HELECT, _HROT, _HVIBR, _HFORM);
  CFdouble tempEps1 = temp*(1.0+EPS);
  CFdouble tempEps = temp*EPS;
  FORTRAN_NAME(enthalpy)(WR1, &LWR1, WI, &LWI, &tempEps1, &tempEps1, &tempEps1, &tempEps1,
                         &pressure, _HTOTALP, _HTRANSP, _HELECTP, _HROTP, _HVIBRP, _HFORMP);

  for(int is=0 ; is < _NS; is++) {
        _CPE[is]   = (_HELECTP[is] - _HELECT[is]) / (tempEps);
        _CPR[is]   = (_HROTP[is] - _HROT[is]) / (tempEps);
        _CPV[is]   = (_HVIBRP[is] - _HVIBR[is]) / (tempEps);
        _CPINT[is] = _CPE[is] + _CPR[is] + _CPV[is];
  }

  FORTRAN_NAME(compotol)(X, &TOL, XTOL);
  //Update of kinetic data
  FORTRAN_NAME(collision)(WR1, &LWR1, WR2, &LWR2, &temp, &temp, &ND, X);

  //A. Heavy particles properties

  //Eucken corrections for internal energy (no inelastic collisions)

  FORTRAN_NAME(lambdaint)(WR1, &LWR1, WR2, &LWR2, _CPINT, XTOL, &lambdaInth);
  FORTRAN_NAME(lambdachicg)(WR1, &LWR1, WR2, &LWR2, XTOL, &lambdaTh, _CHIH);


  FORTRAN_NAME(lambdareacg)(WI, &LWI, WR1, &LWR1, WR2, &LWR2, _HTOTAL,
			    &temp, XTOL, &lambdaRea);

  //B. Electron gas properties
  if(presenceElectron()){
    int i=3;
    FORTRAN_NAME(lambdae)(WR1, &LWR1, WR2, &LWR2, X, &temp, &lambdaTe, &i);
  }

  return lambdaInth + lambdaTh + lambdaRea + lambdaTe;
}

//////////////////////////////////////////////////////////////////////////////

CFdouble MutationLibrary::lambdaNEQ(CFdouble& temperature,
				    CFdouble& pressure)
{
//composition must be called before!
  CFdouble lambdaTot = 0.0;
  CFdouble lambdaInth = 0.0;
  CFdouble lambdaTh = 0.0;
  CFdouble lambdaTe = 0.0;
  CFdouble temp = temperature;

  if (temp < 100.) {
    //cout << " MutationLibrary::lambdaNEQ -> Watch Out !! T <<<<< 100." << temp << endl;
      temp = 100.;
  }

  CFdouble ND = 0.0;
  FORTRAN_NAME(numberd)(WR1,&LWR1,&pressure,&temp,&temp,X,&ND); //ND needed for collision below
  //  THIS ROUTINE NEEDS to BE REWRITTEN TO ACCOUNT FOR THE IONS
  // MARCO IONS
   //Thermodynamic properties

  // - Species specific heat per unit mole
  FORTRAN_NAME(enthalpy)(WR1, &LWR1, WI, &LWI, &temp, &temp, &temp, &temp, &pressure,
                         _HTOTAL, _HTRANS, _HELECT, _HROT, _HVIBR, _HFORM);
  CFdouble tempEps1 = temp*(1.0+EPS);
  CFdouble tempEps = temp*EPS;
  FORTRAN_NAME(enthalpy)(WR1, &LWR1, WI, &LWI, &tempEps1, &tempEps1, &tempEps1, &tempEps1,
                         &pressure, _HTOTALP, _HTRANSP, _HELECTP, _HROTP, _HVIBRP, _HFORMP);

  for(int is=0 ; is < _NS ; is++) {
    _CPE[is]   = (_HELECTP[is] - _HELECT[is]) / (tempEps);
    _CPR[is]   = (_HROTP[is] - _HROT[is]) / (tempEps);
    _CPV[is]   = (_HVIBRP[is] - _HVIBR[is]) / (tempEps);
    _CPINT[is] = _CPE[is] + _CPR[is] + _CPV[is];
  }

  FORTRAN_NAME(compotol)(X, &TOL, XTOL);
  //Update of kinetic data
  FORTRAN_NAME(collision)(WR1, &LWR1, WR2, &LWR2, &temp, &temp, &ND, X);

  //A. Heavy particles properties

  //Eucken corrections for internal energy (no inelastic collisions)
  FORTRAN_NAME(lambdaint)(WR1, &LWR1, WR2, &LWR2, _CPINT, XTOL, &lambdaInth);

  if (_lambdaAlgo == LAMBDACG) {
    FORTRAN_NAME(lambdachicg)(WR1, &LWR1, WR2, &LWR2, XTOL, &lambdaTh, _CHIH);
  }
  else {
    FORTRAN_NAME(lambdachid)(WR1, &LWR1, WR2, &LWR2, XTOL, &lambdaTh, _CHIH);
  }

  //B. Electron gas properties
  if(presenceElectron()){
    int i=3;
    FORTRAN_NAME(lambdae)(WR1, &LWR1, WR2, &LWR2, X, &temp, &lambdaTe, &i);
  }

  lambdaTot = lambdaInth + lambdaTh + lambdaTe;

  return lambdaTot;
}

//////////////////////////////////////////////////////////////////////////////

void MutationLibrary::lambdaVibNEQ(CFreal& temperature,
				   RealVector& tVec,
				   CFdouble& pressure,
				   CFreal& lambdaTrRo,
				   RealVector& lambdaVib)
{
  if (_NS >= 5) {
    CFdouble lambdaRot = 0.0;
    CFdouble lambdaTh  = 0.0;
    CFdouble lambdaE   = 0.0;
    CFdouble lambdaTe  = 0.0;
    CFdouble temp = temperature;
    if (temp < 100.) {
      //cout << " MutationLibrary::lambdaVibNEQ -> Watch Out !! T <<<<< 100." << temp << endl;
      temp = 100.;
    }

    for (CFuint i = 0; i < tVec.size(); ++i) {
      _TVARRAY[i] = tVec[i];
      _TVIBEPS1[i] = tVec[i] *(1.0+EPS);
      _TVIBEPS[i] = tVec[i]*EPS;
    }

    CFreal Te = getTe(temp,&tVec[0]);
    CFdouble tempEps1 = temp*(1.0+EPS);
    CFdouble tempEps = temp*EPS;
    CFdouble TeEps1 = Te*(1.0+EPS);
    CFdouble TeEps = Te*EPS;

    CFdouble ND = 0.0;
    FORTRAN_NAME(numberd)(WR1,&LWR1,&pressure,&temp,&Te,X,&ND);
    //ND needed for collision below

    //Thermodynamic properties
    if(presenceElectron()){
      // - Species specific heat per unit mole
      FORTRAN_NAME(enthalpyvibions)
	(WR1, &LWR1, WI, &LWI, &_nbTvib, &temp, &Te, &temp, _TVARRAY, &pressure,
	 _HTOTAL, _HTRANS, _HELECT, _HROT, _HVIBR, _HFORM);

      FORTRAN_NAME(enthalpyvibions)
	(WR1, &LWR1, WI, &LWI, &_nbTvib, &tempEps1, &TeEps1, &tempEps1, _TVIBEPS1,
	 &pressure, _HTOTALP, _HTRANSP, _HELECTP, _HROTP, _HVIBRP, _HFORMP);

      for(int is=0 ; is < _NS ; is++) {
	_CPE[is]   = (_HELECTP[is] - _HELECT[is])   / (TeEps);
	_CPR[is]   = (_HROTP[is] - _HROT[is])       / (tempEps);
      }

      if (_nbTvib == 2) {
	// rotational NO
	_CPR[4] += (_HVIBRP[4] - _HVIBR[4]) / (tempEps);
	// N2
	_CPVIB[0] = (_HVIBRP[3] - _HVIBR[3]) / (_TVIBEPS[0]);
	// O2
	_CPVIB[1] = (_HVIBRP[5] - _HVIBR[5]) / (_TVIBEPS[1]);
	// NO
	//_CPVIB[2]   = (_HVIBRP[3] - _HVIBR[3]) / (tempEps);
      }
      else if (_nbTvib == 1) {
	for (int is = 0; is <_NS; ++is) {
	  _CPINT[is] = (_HVIBRP[is] - _HVIBR[is]) / (_TVIBEPS[0]);
        }
      }
    }
    else {
      // - Species specific heat per unit mole
      FORTRAN_NAME(enthalpyvib)(WR1, &LWR1, WI, &LWI, &_nbTvib, &temp, &Te, &temp, _TVARRAY, &pressure,
				_HTOTAL, _HTRANS, _HELECT, _HROT, _HVIBR, _HFORM);

      FORTRAN_NAME(enthalpyvib)(WR1, &LWR1, WI, &LWI, &_nbTvib, &tempEps1, &TeEps, &tempEps1, _TVIBEPS1,
				&pressure, _HTOTALP, _HTRANSP, _HELECTP, _HROTP, _HVIBRP, _HFORMP);

      for(int is=0 ; is < _NS ; is++) {
	_CPE[is]   = (_HELECTP[is] - _HELECT[is])   / (TeEps);
	_CPR[is]   = (_HROTP[is] - _HROT[is])       / (tempEps);
      }

      if (_nbTvib == 2) {
	// rotational NO
	_CPR[3] += (_HVIBRP[3] - _HVIBR[3]) / (_TVIBEPS[0]);
	// N2
	_CPVIB[0] = (_HVIBRP[2] - _HVIBR[2]) / (_TVIBEPS[0]);
	// O2
	_CPVIB[1] = (_HVIBRP[4] - _HVIBR[4]) / (_TVIBEPS[1]);
	// NO
	//_CPVIB[2]   = (_HVIBRP[3] - _HVIBR[3]) / (tempEps);
      }
      else if (_nbTvib == 1) {
	for (int is = 0; is <_NS; ++is) {
	  _CPINT[is] = (_HVIBRP[is] - _HVIBR[is]) / (_TVIBEPS[0]);
        }
      }
    }

    FORTRAN_NAME(compotol)(X, &TOL, XTOL);
    //Update of kinetic data

    // here there was  Te = temp
    FORTRAN_NAME(collision)(WR1, &LWR1, WR2, &LWR2, &temp, &Te, &ND, X);

    //A. Heavy particles properties
    //Eucken corrections for rotational internal energy (no inelastic collisions)
    FORTRAN_NAME(lambdaint)(WR1, &LWR1, WR2, &LWR2, _CPR, XTOL, &lambdaRot);

    if (!_noElectEnergy) {
      FORTRAN_NAME(lambdaint)(WR1, &LWR1, WR2, &LWR2, _CPE, XTOL, &lambdaE);
    }

    if(presenceElectron()){
    if(_nbTvib == 1){
      CFdouble lambdaTmp = 0.0;
      FORTRAN_NAME(lambdaint)(WR1, &LWR1, WR2, &LWR2,_CPINT, XTOL,&lambdaTmp);
      _LAMBDAVIB[0] = lambdaTmp;
    }
      // for N2 and O2
    else {
      FORTRAN_NAME(lambdavibions)(WR1, &LWR1, WR2, &LWR2,_CPVIB, XTOL,_LAMBDAVIB);
       }
    }
    else {
    if(_nbTvib == 1){
      CFdouble lambdaTmp = 0.0;
      FORTRAN_NAME(lambdaint)(WR1, &LWR1, WR2, &LWR2,_CPINT, XTOL,&lambdaTmp);
    }
    else {
      FORTRAN_NAME(lambdavib)(WR1, &LWR1, WR2, &LWR2,_CPVIB, XTOL,_LAMBDAVIB);
      }
    }

    if (_lambdaAlgo == LAMBDACG) {
      FORTRAN_NAME(lambdachicg)(WR1, &LWR1, WR2, &LWR2, XTOL, &lambdaTh, _CHIH);
    }
    else {
      FORTRAN_NAME(lambdachid)(WR1, &LWR1, WR2, &LWR2, XTOL, &lambdaTh, _CHIH);
    }

    //B. Electron gas properties
    if(presenceElectron()){
      int i=3;
      FORTRAN_NAME(lambdae)(WR1, &LWR1, WR2, &LWR2, XTOL, &Te, &lambdaTe, &i);
      lambdaTrRo = lambdaRot + lambdaTh;
      _LAMBDAVIB[0] += lambdaTe + lambdaE;
    }
    else {
      lambdaTrRo = lambdaRot + lambdaTh;
    }
    for (int i = 0; i < _nbTvib; i++ ){
      lambdaVib[i] = _LAMBDAVIB[i];
    }
  }
  else {
    lambdaVibNEQN2(temperature, tVec, pressure,
		   lambdaTrRo, lambdaVib);
  }

/*  if (_TVARRAY[0] > 4000) {
  cout <<"lambdaTrRo" << lambdaTrRo<<endl;
  cout <<"lambdaVIB" << _LAMBDAVIB[0]<<endl;
  cout <<"lambdaTe" << _LAMBDAVIB[0]<<endl;
  cout <<"lambdaVIB" << _LAMBDAVIB[0]<<endl;
  //abort();
  } */
}

//////////////////////////////////////////////////////////////////////////////

void MutationLibrary::lambdaVibNEQN2(CFreal& temperature,
				   RealVector& tVec,
				   CFdouble& pressure,
				   CFreal& lambdaTrRo,
				   RealVector& lambdaVib)
{
  CFdouble lambdaRot = 0.0;
  CFdouble lambdaTh  = 0.0;
  CFdouble lambdaE   = 0.0;
  CFdouble lambdaVV   = 0.0;
  // CFdouble lambdaTe  = 0.0;
  CFdouble temp = temperature;

  if (temp < 100.) {
    //cout << " MutationLibrary::lambdaVibNEQN2 -> Watch Out !! T <<<<< 100." << temp << endl;
    temp = 100.;
  }

  for (CFuint i = 0; i < tVec.size(); ++i) {
    _TVARRAY[i] = tVec[i];
    _TVIBEPS1[i] = tVec[i] *(1.0+EPS);
    _TVIBEPS[i] = tVec[i]*EPS;
  }

  CFreal TvN2 = _TVARRAY[0];
  CFdouble tempEps1 = temp*(1.0+EPS);
  CFdouble tempEps = temp*EPS;
  CFdouble TvN2Eps1 = TvN2*(1.0+EPS);
  CFdouble TvN2Eps = TvN2*EPS;

  CFdouble ND = 0.0;
  FORTRAN_NAME(numberd)(WR1,&LWR1,&pressure,&temp,&TvN2,X,&ND); //ND needed for collision below

  //Thermodynamic properties
  FORTRAN_NAME(enthalpy)(WR1, &LWR1, WI, &LWI, &temp, &TvN2, &temp, &TvN2 , &pressure,
			 _HTOTAL, _HTRANS, _HELECT, _HROT, _HVIBR, _HFORM);

  FORTRAN_NAME(enthalpy)(WR1, &LWR1, WI, &LWI, &tempEps1, &TvN2Eps1, &tempEps1, &TvN2Eps1,
			 &pressure, _HTOTALP, _HTRANSP, _HELECTP, _HROTP, _HVIBRP, _HFORMP);

  for(int is=0 ; is < _NS ; is++) {
    _CPR[is]   = (_HROTP[is] - _HROT[is])/ (tempEps);
  }

  for(int is=0 ; is < _NS; is++) {
    _CPE[is]   = (_HELECTP[is] - _HELECT[is]) / (TvN2Eps);
    _CPR[is]   = (_HROTP[is] - _HROT[is]) / (tempEps);
    _CPV[is]   = (_HVIBRP[is] - _HVIBR[is]) / (TvN2Eps);
    _CPINT[is] = _CPE[is] +  _CPV[is];
  }

  FORTRAN_NAME(compotol)(X, &TOL, XTOL);

  FORTRAN_NAME(collision)(WR1, &LWR1, WR2, &LWR2, &temp, &TvN2, &ND, X);

  //A. Heavy particles properties
  //Eucken corrections for rotational internal energy (no inelastic collisions)
  FORTRAN_NAME(lambdaint)(WR1, &LWR1, WR2, &LWR2, _CPR, XTOL, &lambdaRot);

  if (!_noElectEnergy) {
    FORTRAN_NAME(lambdaint)(WR1, &LWR1, WR2, &LWR2, _CPE, XTOL, &lambdaE);
  }

  FORTRAN_NAME(lambdaint)(WR1, &LWR1, WR2, &LWR2, _CPV, XTOL, &lambdaVV);

  if (_lambdaAlgo == LAMBDACG) {
    FORTRAN_NAME(lambdachicg)(WR1, &LWR1, WR2, &LWR2, XTOL, &lambdaTh, _CHIH);
  }
  else {
    FORTRAN_NAME(lambdachid)(WR1, &LWR1, WR2, &LWR2, XTOL, &lambdaTh, _CHIH);
  }

  //B. Electron gas properties
  lambdaTrRo = lambdaRot + lambdaTh;
  lambdaVib[0] =  lambdaVV + lambdaE;
}

//////////////////////////////////////////////////////////////////////////////

CFdouble MutationLibrary::etaCG(CFdouble& temperature,
				CFdouble& pressure,
				CFreal* tVec)
{
  //composition must be called before!
  CFdouble eta = 0.0;
  CFdouble ND = 0.0;
  CFdouble temp = temperature;

  if (temp < 100.) {
    //cout << " MutationLibrary::etaCG -> Watch Out !! T <<<<< 100." << temp << endl;
    temp = 100.;
  }

  CFdouble Te = getTe(temp, tVec);
  FORTRAN_NAME(numberd)(WR1,&LWR1,&pressure,&temp,&Te,X,&ND); //ND needed for collision below

  FORTRAN_NAME(compotol)(X, &TOL, XTOL);

  //Update of kinetic data
  FORTRAN_NAME(collision)(WR1, &LWR1, WR2, &LWR2, &temp, &Te, &ND, X);

  //A. Heavy particles properties

  //Eucken corrections for internal energy (no inelastic collisions)
  FORTRAN_NAME(etacg)(WR1, &LWR1, WR2, &LWR2, XTOL, &eta);

/*  if (temp > 8000.0) {
    cout << "temp --> " << temp << endl;
    cout << "pressure --> " << pressure << endl;
    cout << "ETACG --> " << eta << endl;
//abort();
   } */
  
  return eta;
}

//////////////////////////////////////////////////////////////////////////////

CFdouble MutationLibrary::sigma(CFdouble& temp,   //electrical conductivity
                                CFdouble& pressure,
				CFreal* tVec)
{
  //composition must be called before!
  CFdouble sigma = 0.0;
  CFdouble ND = 0.0;
  CFdouble Te = getTe(temp, tVec);
  
  //  cout << "[T, Te, P] = [" << temp << ", " << Te  << ", " << pressure << "]\n";
  //   cout << "X = ";  
  //   CFreal sumX = 0.;
  //   for (int i = 0; i < _NS; ++i) {
  //     cout << X[i] << " ";
  //     sumX += X[i];
  //   }
  //   cout << "\n";
  //   cout << "sumX = " << sumX << endl;  
  
  FORTRAN_NAME(numberd)(WR1,&LWR1,&pressure,&temp,&Te,X,&ND); //ND needed for collision below
  
  FORTRAN_NAME(compotol)(X, &TOL, XTOL);
  //Update of kinetic data
  FORTRAN_NAME(collision)(WR1, &LWR1, WR2, &LWR2, &temp, &Te, &ND, X);

  int i=2;
  FORTRAN_NAME(sigmae)(WR1, &LWR1, WR2, &LWR2, XTOL, &temp, &sigma, &i);

  return sigma;
}

//////////////////////////////////////////////////////////////////////////////

void MutationLibrary::gammaAndSoundSpeed(CFdouble& temp,
					 CFdouble& pressure,
					 CFdouble& rho,
					 CFdouble& gamma,
					 CFdouble& soundSpeed)
{
  // compute the ratio of the mixture frozen specific heat in thermal equilibrium
  CFdouble drhodp = 0.0;
  CFdouble eps = 0.1;
  FORTRAN_NAME(equigamma)(WR1, &LWR1, WI, &LWI, &temp, &pressure,
			  &rho, Xn, X, &eps, &gamma, &drhodp);
  soundSpeed = sqrt(gamma/drhodp);
}

//////////////////////////////////////////////////////////////////////////////

void MutationLibrary::frozenGammaAndSoundSpeed(CFdouble& temp,
					       CFdouble& pressure,
					       CFdouble& rho,
					       CFdouble& gamma,
					       CFdouble& soundSpeed,
					       RealVector* tVec)
{
  // compute the ratio of the mixture frozen specific heat in thermal equilibrium
  
  if (tVec == CFNULL) {
      FORTRAN_NAME(frozengammafast)(WR1, &LWR1, WI, &LWI, &temp, &pressure,
  				  X, &EPS, &gamma);
    }
    else {
      for (CFuint i = 0; i < tVec->size(); ++i) {
        _TVARRAY[i] = (*tVec)[i];
      }
      CFdouble Te = _TVARRAY[0]; // AL : rough assumption here !!!
      FORTRAN_NAME(frozengammaneqfast)(WR1, &LWR1, WI, &LWI, &_nbTvib, &temp,
  				     &Te, _TVARRAY, &pressure, X, &EPS, &gamma);
    }
  
    soundSpeed = sqrt(gamma*pressure/rho);
  
 //  CFreal numBeta = 0.;
 //  CFreal denBeta = 0.;
 //  const CFuint start = (presenceElectron()) ? 1 : 0;
 //  for (int i = start; i < _NS; ++i) {
 //    const CFreal sigmai = Y[i]/MOLARMASSP[i];
 //    numBeta += sigmai;
 //    denBeta += sigmai*_atomicityCoeff[i];
 //  }
   
   
  // gamma = 1 + numBeta/denBeta;
  // soundSpeed = std::sqrt(gamma*pressure/rho);
}

//////////////////////////////////////////////////////////////////////////////

CFdouble MutationLibrary::soundSpeed(CFdouble& temp,
                                     CFdouble& pressure)
{
  if (!_useLookUpTable) {
    // compute the ratio of the mixture frozen specific heat in thermal equilibrium  c
    CFdouble gamma = 0.0;
    CFdouble drhodp = 0.0;
    CFdouble rho = density(temp, pressure, CFNULL);
    CFdouble eps = 0.1;
    FORTRAN_NAME(equigamma)(WR1, &LWR1, WI, &LWI, &temp, &pressure,
                            &rho, Xn, X, &eps, &gamma, &drhodp);
    return sqrt(gamma/drhodp);
  }
  else {
    return _lookUpTables.get(temp, pressure, _nameToIdxVar.find("a"));
  }
}

//////////////////////////////////////////////////////////////////////////////

void MutationLibrary::setComposition(CFdouble& temp,
				     CFdouble& pressure,
				     RealVector* x)
{
  if (!_useLookUpTable) {
    // initialization of the molar fractions
    for(int i = 0; i < _NS; ++i) {
      Xini[i] = 1.0;
    }

    // compute the molar fractions corresponding to the given temperature and
    // pressure for the mixture Xn
    FORTRAN_NAME(composition)(WR1,&LWR1,WI,&LWI,&temp,&pressure,Xn,Xini,X);

    if (x != CFNULL) {
      for(int i = 0; i < _NS; ++i) {
	(*x)[i] = static_cast<CFreal>(X[i]);
      }
    }
    
    // set mass fractions which will be used later
    CFreal massTot = 0.;
    for (int is = 0; is < _NS; ++is) {
      if (X[is] > 1.00000000001) {
	cout << "X[" << is << "] = " << X[is] << endl;
	abort();
      }
      const CFreal mm = X[is]*MOLARMASSP[is];
      massTot += mm;
      Y[is] = mm;
    }
    
    for (int is = 0; is < _NS; ++is) {
      Y[is] /= massTot;
    }    
  }
}

//////////////////////////////////////////////////////////////////////////////

void MutationLibrary::setDensityEnthalpyEnergy(CFdouble& temp,
					       CFdouble& pressure,
                                               RealVector& dhe)
{
  //setComposition must be called before
  if (!_useLookUpTable) {
    CFdouble ND = 0.0;
    CFdouble rho = 0.0;
    CFdouble MMass = 0.0;

    FORTRAN_NAME(numberd)(WR1,&LWR1,&pressure,&temp,&temp,X,&ND);       //ND

    FORTRAN_NAME(density)(WR1,&LWR1,X,&ND,&rho);        //rho

    FORTRAN_NAME(molarmass)(WR1,&LWR1,&rho,&ND,&MMass);         //MMass

    FORTRAN_NAME(enthalpy)(WR1,&LWR1,WI,&LWI,&temp,&temp,&temp,&temp,&pressure,
                           _HTOTAL,_HTRANS,_HELECT,_HROT,_HVIBR,_HFORM);

    const CFreal pOvRhoMM = pressure/rho*MMass;
    // sum up all the internal energies for each species
    dhe[1] = 0.0;
    dhe[2] = 0.0;
    for(int i = 0; i < _NS; ++i) {
      dhe[1] += X[i]*_HTOTAL[i];
      dhe[2] += X[i]*(_HTOTAL[i] - pOvRhoMM);
    }
    
    dhe[0] = rho;
    dhe[1] /= MMass;
    dhe[2] /= MMass;
  }
  else {
    dhe[0] = _lookUpTables.get(temp, pressure, _nameToIdxVar.find("d"));
    dhe[1] = _lookUpTables.get(temp, pressure, _nameToIdxVar.find("h"));
    dhe[2] = _lookUpTables.get(temp, pressure, _nameToIdxVar.find("e"));
  }
}

//////////////////////////////////////////////////////////////////////////////

void MutationLibrary::setDensityEnthalpyEnergy(CFdouble& temp,
					       RealVector& tVec,
                                               CFdouble& pressure,
                                               RealVector& dhe,
					       bool storeExtraData)
{
  CFdouble ND = 0.0;
  CFdouble rho = 0.0;
  CFdouble MMass = 0.0;

  FORTRAN_NAME(numberd)(WR1,&LWR1,&pressure,&temp,&temp,X,&ND);       //ND

  FORTRAN_NAME(density)(WR1,&LWR1,X,&ND,&rho);        //rho

  FORTRAN_NAME(molarmass)(WR1,&LWR1,&rho,&ND,&MMass);         //MMass
  const CFreal pOvRhoMM = pressure/rho*MMass;
  const CFreal ovMMass = 1./MMass;
  for (CFuint i = 0; i < tVec.size(); ++i) {
    _TVARRAY[i] = tVec[i];
  }

  if (presenceElectron()) {
    CFdouble TEL = _TVARRAY[0];  // I impose that the electronic temperature is equal
    // to the vibration of N2

    FORTRAN_NAME(enthalpyvibions)(WR1,&LWR1,WI,&LWI, &_nbTvib, &temp, &TEL,
				  &temp,_TVARRAY,&pressure,
				  _HTOTAL,_HTRANS,_HELECT,_HROT,_HVIBR,_HFORM);
    // sum up all the internal energies for each species
    dhe = 0.0;

    for(int i = 0; i < _NS; ++i) {
      dhe[1] += X[i]*_HTOTAL[i];
      dhe[2] += X[i]*(_HTOTAL[i] - pOvRhoMM);
      if (_nbTvib == 1) {
	dhe[3] += X[i]*_HVIBR[i]*ovMMass;
      }
    }

    if (_nbTvib == 2) {
      dhe[3] = X[3]*_HVIBR[3]/MOLARMASSP[3];
      dhe[4] = X[5]*_HVIBR[5]/MOLARMASSP[5];
    }

    dhe[0] = rho;
    dhe[1] *= ovMMass;
    dhe[2] *= ovMMass;
  }
  else {
    //
     CFdouble TEL = _TVARRAY[0];  // MP: electronic temperature is equal to the vibration of N2

    if (_NS >= 5) {
      FORTRAN_NAME(enthalpyvib)(WR1,&LWR1,WI,&LWI, &_nbTvib, &temp, &TEL,
				&temp,_TVARRAY,&pressure, _HTOTAL,_HTRANS,
				_HELECT,_HROT,_HVIBR,_HFORM);
    }
    else {
      FORTRAN_NAME(enthalpy)(WR1,&LWR1,WI,&LWI, &temp, &TEL, &temp, &TEL,
			     &pressure, _HTOTAL,_HTRANS,_HELECT,_HROT,
			     _HVIBR,_HFORM);
    }

    // sum up all the internal energies for each species
    dhe = 0.0;

    for(int i = 0; i < _NS; ++i) {
      dhe[1] += X[i]*_HTOTAL[i];
      dhe[2] += X[i]*(_HTOTAL[i] - pOvRhoMM);

      // to be tested
      if (_noElectEnergy) {
	dhe[1] -= X[i]*_HELECT[i];
	dhe[2] -= X[i]*_HELECT[i];
      }

      if (_nbTvib == 1) {
	dhe[3] += X[i]*_HVIBR[i]*ovMMass;
      }
    }

    if (_nbTvib == 2) {
      dhe[3] = X[2]*_HVIBR[2]/MOLARMASSP[2];
      dhe[4] = X[4]*_HVIBR[4]/MOLARMASSP[4];
    }

    dhe[0] = rho;
    dhe[1] *= ovMMass;
    dhe[2] *= ovMMass;

    // store the species translational-rotational energies
    if (storeExtraData) {
      cf_assert(_nbTvib == 1);
      _extraData.dEdT = 0.0;
      // store the roto-translational energy (plus formation enthalpy)
      for(int i = 0; i < _NS; ++i) {
	_extraData.energyTr[i] = (_HTRANS[i] - pOvRhoMM + _HROT[i] + _HFORM[i])/MOLARMASSP[i];
	_extraData.dRhoEdRhoi[i] = _extraData.energyTr[i];
	_extraData.energyVib[i] = _HVIBR[i]/MOLARMASSP[i];
	
	if (_flagMoleculesIDs[i]) {
	  _extraData.dRhoEvdRhoi[i] = _HVIBR[i]/MOLARMASSP[i];
	  _extraData.dEdT += 2.5*X[i];
	}
	else {
	  _extraData.dEdT += 1.5*X[i];
	}
      }

      _extraData.dEdT *= _Rgas*ovMMass;

      // perturb the vibrational temperature
      CFdouble startTv = _TVARRAY[0];
      _TVARRAY[0] = startTv*(1.+EPS);
      CFdouble epsTv = startTv*EPS;
      for (int is = 0; is < _NS; ++is) {
	_tmp[is] = _HVIBR[is];
      }

      if (_NS >= 5) {
	FORTRAN_NAME(enthalpyvib)(WR1,&LWR1,WI,&LWI, &_nbTvib, &temp, &TEL,
				  &temp,_TVARRAY,&pressure, _HTOTAL,_HTRANS,
				  _HELECT,_HROT,_HVIBR,_HFORM);
      }
      else {
	FORTRAN_NAME(enthalpy)(WR1,&LWR1,WI,&LWI, &temp, &TEL, &temp, &_TVARRAY[0],
			       &pressure, _HTOTAL,_HTRANS,_HELECT,_HROT,
			       _HVIBR,_HFORM);
      }

      _extraData.dEvTv = 0.0;
      for(int is = 0; is < _NS; ++is) {
	_extraData.dEvTv += X[is]*(_HVIBR[is] - _tmp[is]);
      }
      _extraData.dEvTv /= (epsTv*MMass);
    }
  }

  // Checked OK Compared with MUTATION 2
/*   if (temp > 8000.0) {
  cout << "temp --> " << temp << endl;
  cout << "pressure --> " << pressure << endl;
  for(int i = 0; i < _NS; ++i) {
    cout << "HTOTA --> " << _HTOTAL[i] << endl;
  }
  for(int i = 0; i < _NS; ++i) {
    cout << "HTRANS --> " << _HTRANS[i] << endl;
  }
  for(int i = 0; i < _NS; ++i) {
    cout << "HELECT --> " << _HELECT[i] << endl;
  }
  for(int i = 0; i < _NS; ++i) {
    cout << "HVIBR --> " << _HVIBR[i] << endl;
  }
  for(int i = 0; i < _NS; ++i) {
    cout << "HFORM --> " << _HFORM[i] << endl;
    cout << "MOLARMASS --> " << MOLARMASSP[i] << endl;
  }

  cout << "dhe = " << dhe << endl;
  //abort();
  } */
}

//////////////////////////////////////////////////////////////////////////////

CFdouble MutationLibrary::density(CFdouble& temp,
                                  CFdouble& pressure,
				  CFreal* tVec)
{
  //  CF_DEBUG_POINT;
  if (!_useLookUpTable) {
    CFdouble ND = 0.0;
    CFdouble rho = 0.0;

    CFdouble Te = getTe(temp,tVec);
    FORTRAN_NAME(numberd)(WR1,&LWR1,&pressure,&temp,&Te,X,&ND);
    FORTRAN_NAME(density)(WR1,&LWR1,X,&ND,&rho);
    return rho;
  }
  return _lookUpTables.get(temp, pressure, _nameToIdxVar.find("d"));
}
//////////////////////////////////////////////////////////////////////////////

void MutationLibrary::setSpeciesMolarFractions(const RealVector& xs)
{
  for (int is = 0; is < _NS; ++is) {
    X[is] = xs[is];
     if (X[is] < 0.0) X[is] = 0.0;
    if (X[is] > 1.00000000001) {
       cout << "X[is] = " << X[is] << endl;
      // abort();
    }
  }
    
}

//////////////////////////////////////////////////////////////////////////////

CFdouble MutationLibrary::pressure(CFdouble& rho,
				   CFdouble& temp,
				   CFreal* tVec)
{
  CFdouble p = 0.0;
  if (NE == 0) {
  // set the mixture pressure
    FORTRAN_NAME(pressure)(WR1, &LWR1, &rho, &temp, MOLARMASSP, Y, &p);
  }
  else {
    _electronPress = 0.0;
    // if the electronic temperature is available, it is
    // the last entry in the array tVec (whose size is _nbTvib+1)
    CFdouble Te = getTe(temp, tVec);
    FORTRAN_NAME(epressure)(WR1, &LWR1, &rho, &temp, &Te, MOLARMASSP, Y,
			    &p, &_electronPress);
  }

  return p;
}

//////////////////////////////////////////////////////////////////////////////

CFdouble MutationLibrary::electronPressure(CFreal rhoE,
					   CFreal tempE)
{
  return rhoE*tempE*_Rgas/MOLARMASSP[0];
}

//////////////////////////////////////////////////////////////////////////////

CFdouble MutationLibrary::energy(CFdouble& temp,
                                 CFdouble& pressure)
{
  if (!_useLookUpTable) {
    FORTRAN_NAME(energy)(WR1,&LWR1,WI,&LWI,&temp,&temp,&temp,&temp,&pressure,
                         _ETOTAL,_ETRANS,_EELECT,_EROT,_EVIBR,_EFORM);

    CFdouble ND = 0.0;
    CFdouble rho = 0.0;
    CFdouble MMass = 0.0;

    // store the density
    FORTRAN_NAME(numberd)(WR1,&LWR1,&pressure,&temp,&temp,X,&ND);
    FORTRAN_NAME(density)(WR1,&LWR1,X,&ND,&rho);
    FORTRAN_NAME(molarmass)(WR1,&LWR1,&rho,&ND,&MMass);

    // sum up all the internal energies for each species
    CFdouble intEnergy = 0.0;
    for(int i = 0; i < _NS; ++i) {
      intEnergy += X[i]*_ETOTAL[i];
    }
    return intEnergy /= MMass;
  }
  else {
    return _lookUpTables.get(temp, pressure, _nameToIdxVar.find("e"));
  }
}

//////////////////////////////////////////////////////////////////////////////

CFdouble MutationLibrary::enthalpy(CFdouble& temp,
                                   CFdouble& pressure)
{
  if (!_useLookUpTable) {
    FORTRAN_NAME(enthalpy)(WR1,&LWR1,WI,&LWI,&temp,&temp,&temp,&temp,&pressure,
                           _HTOTAL,_HTRANS,_HELECT,_HROT,_HVIBR,_HFORM);

    CFdouble ND = 0.0;
    CFdouble rho = 0.0;
    CFdouble MMass = 0.0;

    // store the density
    FORTRAN_NAME(numberd)(WR1,&LWR1,&pressure,&temp,&temp,X,&ND);
    FORTRAN_NAME(density)(WR1,&LWR1,X,&ND,&rho);
    FORTRAN_NAME(molarmass)(WR1,&LWR1,&rho,&ND,&MMass);

    // sum up all the internal energies for each species
    CFdouble h = 0.0;
    for(int i = 0; i < _NS; ++i) {
      h += X[i]*_HTOTAL[i];
    }
    return h /= MMass;
  }
  else {
    return _lookUpTables.get(temp, pressure, _nameToIdxVar.find("e"));
  }
}

//////////////////////////////////////////////////////////////////////////////

void MutationLibrary::setElemFractions(const RealVector& yn)
{
  for (int ic = 0; ic < _NC; ++ic) {
    Yn[ic] = yn[ic];

    if (!(Yn[ic] >= 0.0 && Yn[ic] <= 1.0)) {
      // cout << "Yn[ic] = " << Yn[ic] << endl;
      // abort();
    }

    cf_assert(Yn[ic] >= 0.0);
    cf_assert(Yn[ic] <= 1.0);
  }

  // Fills Xn according to Yn
  FORTRAN_NAME(nucmasstomolfrac)(WR1, &LWR1, WI, &LWI, Yn, Xn);
}

//////////////////////////////////////////////////////////////////////////////

void MutationLibrary::setElementXFromSpeciesY(const RealVector& ys)
{
  cf_assert(ys.size() == static_cast<CFuint>(_NS));
  cf_assert(_NC == 2);

  CFreal massTot = 0.;
  for (int is = 0; is < _NS; ++is) {
    if (ys[is] > 1.00000000001) {
      // cout << "ys[is] = " << ys[is] << endl;
      // abort();
    }

    const CFreal mm = ys[is]/MOLARMASSP[is];
    massTot += mm;
    X[is] = mm;
  }

  massTot = 1./massTot;
  for (int is = 0; is < _NS; ++is) {
    X[is] *= massTot;
  }

  FORTRAN_NAME(compfract)(X,&NE,Xn);
}

//////////////////////////////////////////////////////////////////////////////

void MutationLibrary::setElectronFraction(RealVector& ys)
{
  // get molar mass of species
  FORTRAN_NAME(setspeciesmolarmass)(WR1, &LWR1, MOLARMASSP);

  // charge neutrality: xEl = sum(xIon)
  CFdouble yEl = 0.0;
  for (int is = 0; is < _NS; ++is) {
    if (_CHARGE[is] > 0) {
      yEl += ys[is] / MOLARMASSP[is];
    }
  }

  yEl *= MOLARMASSP[0]; // 1st species: electron
  ys[0] = yEl; // overwrite electron mass fraction
}

//////////////////////////////////////////////////////////////////////

void MutationLibrary::setSpeciesFractions(const RealVector& ys)
{
  // test this !!!!
  if (presenceElectron()) {
    setElectronFraction(const_cast<RealVector&>(ys));
  }

  for (int is = 0; is < _NS; ++is) {
    Y[is] = ys[is];

    if (Y[is] < 0.0) Y[is] = 0.0;
    if (Y[is] > 1.00000000001) {
      // cout << "Y[is] = " << Y[is] << endl;
      // abort();
    }
  }
    
  // Fills X according to Y
  FORTRAN_NAME(specmasstomolfrac)(WR1, &LWR1, WI, &LWI, Y, X);
  
  
}

//////////////////////////////////////////////////////////////////////////////

void MutationLibrary::getSpeciesMolarFractions
(const RealVector& ys, RealVector& xs)
{
  CFreal massTot = 0.;
  for (int is = 0; is < _NS; ++is) {
    if (ys[is] > 1.00000000001) {
      //cout << "ys[is] = " << ys[is] << endl;
      //abort();
    }

    const CFreal mm = ys[is]/MOLARMASSP[is];
    massTot += mm;
    xs[is] = mm;
  }
  xs *= 1./massTot;

  //   for (int is = 0; is < _NS; ++is) {
  //     Y[is] = ys[is];

  //     if (Y[is] < 0.0) Y[is] = 0.0;
  //     if (Y[is] > 1.0) {
  //       cout << "Y[is] = " << Y[is] << endl;
  //       abort();
  //     }
  //   }

  //   // Fills X according to Y
  //   FORTRAN_NAME(specmasstomolfrac)(WR1, &LWR1, WI, &LWI, Y, X);

  //   for (int is = 0; is < _NS; ++is) {
  //     xs[is] = X[is];
  //   }
}

//////////////////////////////////////////////////////////////////////////////

void MutationLibrary::getSpeciesMassFractions
(const RealVector& xs, RealVector& ys)
{
  CFreal massTot = 0.;
  for (int is = 0; is < _NS; ++is) {
    if (xs[is] > 1.00000000001) {
      cout << "xs[" << is << "] = " << xs[is] << endl;
      abort();
    }
    
    const CFreal mm = xs[is]*MOLARMASSP[is];
    massTot += mm;
    ys[is] = mm;
  }

  ys /= massTot;
}

//////////////////////////////////////////////////////////////////////////////

void MutationLibrary::getSpeciesMassFractions(RealVector& ys)
{
  for (int is = 0; is < _NS; ++is) {
    ys[is] = Y[is];
  }
}

//////////////////////////////////////////////////////////////////////////////

void MutationLibrary::getTransportCoefs(CFdouble& temp,
					CFdouble& pressure,
					CFdouble& lambda,
					CFdouble& lambdacor,
					RealVector& lambdael,
					RealMatrix& eldifcoef,
					RealVector& eltdifcoef)
{
// NOTE: X must be set before the calling of this function!
  CFdouble ND = 0.0;
  // Cumpute enthalpies
  FORTRAN_NAME(enthalpy)(WR1, &LWR1, WI, &LWI, &temp, &temp, &temp,
			 &temp, &pressure, _HTOTAL, _HTRANS, _HELECT,
			 _HROT, _HVIBR, _HFORM);
  // Compute number density
  FORTRAN_NAME(numberd)(WR1, &LWR1, &pressure, &temp, &temp, X, &ND);

  // Avoid error in case of vanishing molar fractions of species
  //FORTRAN_NAME(compotol2)(X, &_Xlim, XTOL);
  FORTRAN_NAME(compotol)(X, &TOL, XTOL);

  // Fill WR2 vector, necessary if temperature is changed,
  // or composition in case of electron presence is changed ??? Here ???
  FORTRAN_NAME(collision)(WR1, &LWR1, WR2, &LWR2, &temp, &temp, &ND, X);
  // Compute FIJ
  FORTRAN_NAME(correction)(WR1, &LWR1, WR2, &LWR2, XTOL, &_sonine, FIJ);
  // Compute transport properties
  FORTRAN_NAME(ltevefneu)(WI, &LWI, WR1, &LWR1, WR2, &LWR2, _HTOTAL,
			  &temp, X, &ND, FIJ, &lambda, &lambdacor,
			  LAMBDAEL, ELDIFCOEF, ELTDIFCOEF);

  // Copying non scalar data
  for (int ic = 0; ic < _NC; ++ic) {
    lambdael[ic] = LAMBDAEL[ic];
    eltdifcoef[ic] = ELTDIFCOEF[ic];
    for (int jc = 0; jc < _NC; ++jc) {
      // Careful because C++ and Fortran label matrices differently!
      eldifcoef(ic,jc) = ELDIFCOEF[jc*_NC+ic];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void MutationLibrary::getMassProductionTerm(CFdouble& temperature,
					    RealVector& tVec,
					    CFdouble& pressure,
					    CFdouble& rho,
					    const RealVector& ys,
					    bool flagJac,
					    RealVector& omega,
					    RealMatrix& jacobian)
{
  if (!_freezeChemistry) {
    // Limit the mass composition of species !!!
    CFdouble YLIM = 1.0e-15;
    // Flag to Fortran
    // int getJacobian = flagJac ? 1 : 0;
    CFdouble temp = max(_TminFix, temperature);

    // Fill species mass fractions
    for (int is = 0; is < _NS; ++is) {
      Y[is] = ys[is];

      if (Y[is] < 0.0) Y[is] = 0.0;
      if (Y[is] > 1.01) {
      	cout << "Y[is] = " << Y[is] << endl;
	// abort();
      }
      //  cf_assert(Y[is] <= 1.0);
    }

    // Computing the reaction rates (using one temperature now!)
    // FORTRAN_NAME(arrheniuspvty)(WR1, &LWR1, WR3, &LWR3, WI, &LWI,
    //                               &pressure, &temp, Y, &YLIM, &getJacobian,
    //                               OMEGA, DWDP, DWDT, DWDYI);

    for (CFuint i = 0; i < tVec.size(); ++i) {
      _TVARRAY[i] = tVec[i];
    }
    CFdouble TV = tVec[0];

    // This is wrong we should specify that if it is CNEQ then we should call arrhenius
    // on the other hand if it is TCNEQ we should call arrheniusvt because could be the
    // case that Tvib = 1 but I have TCNEQ ! Ask me if you don t understand Marco
    if (tVec.size() == 1) {
      FORTRAN_NAME(arrhenius)(WR1, &LWR1, WR3, &LWR3, WI, &LWI,
			      Y, &YLIM, &pressure, &temp, &TV,
			      &TV, &rho, OMEGA);
    }
    else {
      if (presenceElectron()) {
        FORTRAN_NAME(arrheniusvtions)(WR1, &LWR1, WR3, &LWR3, WI, &LWI,
				      Y, &YLIM, &pressure, &temp, _TVARRAY,
				      &rho, OMEGA);
      }
      else {
        FORTRAN_NAME(arrheniusvt)(WR1, &LWR1, WR3, &LWR3, WI, &LWI,
				  Y, &YLIM, &pressure, &temp, _TVARRAY,
				  &rho, OMEGA);
      }
    }

    // Returning the reaction rates
    for (int is = 0; is < _NS; ++is) {
      omega[is] = _factorOmega*OMEGA[is];
    }
    
    // Assembling Jacobian matrix of the source term corresponding to
    // primitive variables: p, u, v,( w,) T, Ys
    if (flagJac) {
      jacobian = 0.0;
      CFuint NbEqu = jacobian.nbRows();
      cf_assert ( jacobian.nbCols() == NbEqu );
      CFuint NbEuler = NbEqu - _NS;
      for (int is = 0; is < _NS; ++is) {
	// Pressure contribution
	jacobian(NbEuler+is,0) = DWDP[is];
	// Temperature contribution
	jacobian(NbEuler+is,NbEuler-1) = DWDT[is];
	// Species mass fractions contribution
	for (int js = 0; js < _NS; ++js) {
	  // Careful because C++ and Fortran label matrices differently!
	  jacobian(NbEuler+is,NbEuler+js) = DWDYI[js*_NS+is];
	}
      }
    }
  }
  else {
    // Returning the reaction rates
    for (int is = 0; is < _NS; ++is) {
      omega[is] = 0.0;
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void MutationLibrary::getSource(CFdouble& temp,
			        RealVector& tempVib,
			        CFdouble& pressure,
			        CFdouble& rho,
			        const RealVector& ys,
			        bool flagJac,
			        RealVector& omega,
                                RealVector& omegav,
                                CFdouble& omegaRad,
                                RealMatrix& jacobian)
{
} 

//////////////////////////////////////////////////////////////////////////////

void MutationLibrary::getRhoUdiff(CFdouble& temperature,
				  CFdouble& pressure,
                                  RealVector& normConcGradients,
                                  CFreal* tVec,
				  RealVector& rhoUdiff,
				  bool fast)
{
  // NOTE: X must be set before the calling of this function!
  CFdouble ND = 0.0;
  CFdouble rho = 0.0;
  CFdouble MMass = 0.0;
  CFdouble temp = max(_TminFix, temperature);
  CFdouble Te = getTe(temp,tVec);

  // Compute number density
  FORTRAN_NAME(numberd)(WR1, &LWR1, &pressure, &temp, &Te, X, &ND);
  // Avoid error in case of vanishing molar fractions of species


  /// @TODO AL: check this for ionized cases !!!!!!!!!!!!!!!!!!!!!!!!!!

  //  FORTRAN_NAME(compotol2)(X, &_Xlim, XTOL);
  //  if (presenceElectron()) {
  //     FORTRAN_NAME(compotol)(X, &TOL, XTOL);
  //  }
  //  else {
  FORTRAN_NAME(compotol2)(X, &_Xlim, XTOL);
  //  }

  if (!fast) {
    // Fill WR2 vector, necessary if temperature is changed,
    // or composition in case of electron presence is changed ??? Here ???
    FORTRAN_NAME(collision)(WR1, &LWR1, WR2, &LWR2, &temp, &Te, &ND, X);
    // Compute FIJ
    FORTRAN_NAME(correction)(WR1, &LWR1, WR2, &LWR2, XTOL, &_sonine, FIJ);
  }
  
  // density of the mixture
  FORTRAN_NAME(density)(WR1,&LWR1,X,&ND,&rho);
  // set the mixture molar mass
  FORTRAN_NAME(molarmass)(WR1,&LWR1,&rho,&ND,&MMass);

  // Set driving forces as gradients of molar fractions
  CFreal normMMassGradient = 0.0;
  for (int is = 0; is < _NS; ++is) {
    normMMassGradient += normConcGradients[is] / MOLARMASSP[is];
  }
  normMMassGradient *= -MMass*MMass;

  for (int is = 0; is < _NS; ++is) {
    DF[is] = (MMass*normConcGradients[is] + Y[is]*normMMassGradient) /
      MOLARMASSP[is];
  }

  CFdouble eamb = 0.;
  if (presenceElectron()) {
    FORTRAN_NAME(smd)(WR1, &LWR1, WR2, &LWR2, WI, &LWI, XTOL,
		      &temp, &Te, &ND, DF, FIJ, JDIF, &eamb);
  }
  else {
    FORTRAN_NAME(smneutd)(WR1, &LWR1, WR2, &LWR2, XTOL, &ND, DF, FIJ, JDIF);
    //FORTRAN_NAME(smneutsut)(WR1, &LWR1, WR2, &LWR2, XTOL, &ND, DF, FIJ, JDIF);
  }

  for (int is = 0; is < _NS; ++is) {
    rhoUdiff[is] = JDIF[is];
  }

  // if (Te > 4000.) {
//     cout << "Te = " << Te << endl;
//     cout << "temp = " << temp << endl;
//     cout << "rhoU = " << rhoUdiff << endl;
//     abort();
//   }
}

//////////////////////////////////////////////////////////////////////////////

void MutationLibrary::getDij_fick(RealVector& dx,
				  CFdouble& pressure,
				  CFdouble& temperature,
				  RealMatrix& Dij,
                                  RealVector& rhoUdiff)
{
  // NOTE: X must be set before the calling of this function!
  CFdouble ND = 0.0;
  CFdouble MMass = 0.0;
  CFdouble rho = 0.0;
  CFdouble yi = 0.0;

  // Compute number density and density 
  FORTRAN_NAME(numberd)(WR1, &LWR1, &pressure, &temperature, &temperature, X, &ND);
  FORTRAN_NAME(density)(WR1,&LWR1,X,&ND,&rho);
  
   // Fill WR2 vector, necessary if temperature is changed
   FORTRAN_NAME(collision)(WR1, &LWR1, WR2, &LWR2, &temperature, &temperature, &ND, X);
 
  // First we fill the matrix of the binary diffusion coefficients
  CFuint ij;
  for (int is = 1; is < _NS+1; ++is) {
    for (int js = 1; js < _NS+1; ++js) {
          ij = ((is-1)*(2*_NS-is)+2*js)/2;
          Dij(is-1,js-1) = WR2[IBINIJ+ij-1] /ND;
          
    }
  }

  // set the mixture molar mass
  FORTRAN_NAME(molarmass)(WR1,&LWR1,&rho,&ND,&MMass);
  
  // Compute the Diffusion coefficient
  // Since we use Fick law with a constant coefficient for all spicies
  // It is necessary to choose a spicies that is present
    CFreal sum = 0.0;
    if (_NS == 5){
      if (X[0] != 0.0){
	for (int js = 1; js < _NS; ++js)
	  {
	    sum += (X[js]/Dij(0,js));
	  }
	sum *= X[0];
	yi = X[0]*MOLARMASSP[0]/MMass;
      }
      
      else if (X[1] != 0.0){
	for (int js = 0; js < _NS; ++js)
	  { if (js != 1)
	      sum += (X[js]/Dij(1,js));
	  }
	sum *= X[1];
	yi = X[1]*MOLARMASSP[1]/MMass;
	
      }
      
      else if (X[2] != 0.0){
	for (int js = 0; js < _NS; ++js)
	  {	if (js != 2)
	      sum += (X[js]/Dij(2,js));
	  }
	sum *= X[2];
	yi = X[2]*MOLARMASSP[2]/MMass;
	
      }
    
      else if(X[3] != 0.0){
	for (int js = 0; js < _NS; ++js){
	  if (js != 2)
	    sum += (X[js]/Dij(3,js));
	}
	sum *= X[3];
	yi = X[3]*MOLARMASSP[3]/MMass;
	
      }
	
      else if(X[4] != 0.0){
	for (int js = 0; js < _NS; ++js){
	  sum += (X[js]/Dij(4,js));
	}
	sum *= X[4];
	yi = X[4]*MOLARMASSP[4]/MMass;
	
      }
      

	else std::cout<<"There is a bug\n";
    }
    
    else if (_NS==11 ){ 
if (X[1] != 0.0){
	for (int js = 0; js < _NS; ++js)
	  { if (js != 1)
	      sum += (X[js]/Dij(1,js));
	  }
	sum *= X[1];
	yi = X[1]*MOLARMASSP[1]/MMass;
	
      }
      
      else if (X[2] != 0.0){
	for (int js = 0; js < _NS; ++js)
	  {	if (js != 2)
	      sum += (X[js]/Dij(2,js));
	  }
	sum *= X[2];
	yi = X[2]*MOLARMASSP[2]/MMass;
	
      }
    
      else if(X[3] != 0.0){
	for (int js = 0; js < _NS; ++js){
	  if (js != 3)
	    sum += (X[js]/Dij(3,js));
	}
	sum *= X[3];
	yi = X[3]*MOLARMASSP[3]/MMass;
	
      }
	
      else if(X[4] != 0.0){
	for (int js = 0; js < _NS; ++js){
	   if (js != 4)
	     sum += (X[js]/Dij(4,js));
	}
	sum *= X[4];
	yi = X[4]*MOLARMASSP[4]/MMass;
	
      }
      else if(X[5] != 0.0){
	for (int js = 0; js < _NS; ++js){
	  if (js != 5)
	    sum += (X[js]/Dij(5,js));
	}
	sum *= X[5];
	yi = X[5]*MOLARMASSP[5]/MMass;
	
      }
      // else if(X[6] != 0.0){
      // 	for (int js = 0; js < _NS; ++js){
      // 	  if (js != 6)
      // 	  sum += (X[js]/Dij(6,js));
      // 	}
      // 	sum *= X[6];
      // 	yi = X[6]*MOLARMASSP[6]/MMass;
	
      // }
      // else if(X[7] != 0.0){
      // 	for (int js = 0; js < _NS; ++js){
      // 	  if (js != 7)
      // 	  sum += (X[js]/Dij(7,js));
      // 	}
      // 	sum *= X[7];
      // 	yi = X[7]*MOLARMASSP[7]/MMass;
	
      // }
    //   else if(X[8] != 0.0){
  // 	for (int js = 0; js < _NS; ++js){
  // 	  if (js != 8)
  // 	  sum += (X[js]/Dij(8,js));
  // 	}
  // 	sum *= X[8];
  // 	yi = X[8]*MOLARMASSP[8]/MMass;
	
  //     }
  //     else if(X[9] != 0.0){
  // 	for (int js = 0; js < _NS; ++js){
  // 	  if (js != 9)
  // 	  sum += (X[js]/Dij(9,js));
  // 	}
  // 	sum *= X[9];
  // 	yi = X[9]*MOLARMASSP[9]/MMass;
	
  //     }
  // else if(X[10] != 0.0){
  // 	for (int js = 0; js < _NS; ++js){
  // 	  if (js != 10)
  // 	  sum += (X[js]/Dij(10,js));
  // 	}
  // 	sum *= X[10];
  // 	yi = X[10]*MOLARMASSP[10]/MMass;
	
  //     }

	else std::cout<<"There is a bug\n";
      
    }
   
    //    Diff_coeff = rho*(1.0 - yi)/sum;
    
    for (int is = 0; is < _NS; ++is){
      rhoUdiff[is] = -dx[is]*rho*(1.0 - yi)/sum;
    }
}

//////////////////////////////////////////////////////////////////////////////

void MutationLibrary::getSpeciesTotEnthalpies(CFdouble& temp,
					      RealVector& tVec,
                                              CFdouble& pressure,
                                              RealVector& hsTot,
					      RealVector* hsVib,
					      RealVector* hsEl)
{
  if (tVec.size() == 0) {
    FORTRAN_NAME(enthalpy)(WR1,&LWR1,WI,&LWI,&temp,&temp,&temp,&temp,&pressure,
                           _HTOTAL,_HTRANS,_HELECT,_HROT,_HVIBR,_HFORM);
  }
  else {
    for (CFuint i = 0; i < tVec.size(); ++i){
      _TVARRAY[i] = tVec[i];
    }

    // MARCO ION Here we have to check that ENTHALPY VIB!!! CHANGE IT
    // THIS IS FOR TCNEQ
    int sizeTvib = tVec.size();
    if (_NS >= 5) {
      if (presenceElectron()) {
	CFdouble TEL = _TVARRAY[0];  // I impose that the electronic temperature is equal
	// to the vibration of N2
	FORTRAN_NAME (enthalpyvibions)(WR1,&LWR1,WI,&LWI,&sizeTvib,&temp,&TEL,&temp,_TVARRAY,&pressure,
				       _HTOTAL,_HTRANS,_HELECT,_HROT,_HVIBR,_HFORM);

	if (_nbTvib == 2) {
	  (*hsVib)[0] = _HVIBR[3] / MOLARMASSP[3];
	  (*hsVib)[1] = _HVIBR[5] / MOLARMASSP[5];
	}

	if (_nbTvib == 1) {
	  for (CFuint i = 0; i < _moleculesIDs.size(); ++i) {
	    const CFuint molID = _moleculesIDs[i];
	    (*hsVib)[i] = _HVIBR[molID] / MOLARMASSP[molID];
	  }
	}
      }
      else {
	CFdouble TEL = _TVARRAY[0];  // I impose that the electronic temperature is equal
	// to the vibration of N2
	FORTRAN_NAME (enthalpyvib)(WR1,&LWR1,WI,&LWI,&sizeTvib,&temp,&TEL,&temp,_TVARRAY,&pressure,
				   _HTOTAL,_HTRANS,_HELECT,_HROT,_HVIBR,_HFORM);

	if (_nbTvib == 2) {
	  (*hsVib)[0] = _HVIBR[2] / MOLARMASSP[2];
	  (*hsVib)[1] = _HVIBR[4] / MOLARMASSP[4];
        }

        if (_nbTvib == 1) {
	  for (CFuint i = 0; i < _moleculesIDs.size(); ++i) {
	    const CFuint molID = _moleculesIDs[i];
	    (*hsVib)[i] = _HVIBR[molID] / MOLARMASSP[molID];
	  }
	}
      }
    }
    else {
      CFdouble TEL = _TVARRAY[0];  // I impose that the electronic temperature is equal
      // to the vibration of N2

      FORTRAN_NAME (enthalpy)(WR1,&LWR1,WI,&LWI,&temp,&TEL,&temp,&TEL,&pressure,
			      _HTOTAL,_HTRANS,_HELECT,_HROT,_HVIBR,_HFORM);

      (*hsVib)[0] = _HVIBR[1] / MOLARMASSP[1];
    }
  }

  // returning the total enthalpies per unit mass of species
  for(int i = 0; i < _NS; ++i) {
    hsTot[i] = (!_noElectEnergy) ? _HTOTAL[i] / MOLARMASSP[i] :
      (_HTOTAL[i] -_HELECT[i])/MOLARMASSP[i];
  }
}

//////////////////////////////////////////////////////////////////////////////

void MutationLibrary::setLookUpTables()
{
  vector<ComputeQuantity> varComputeVec;
  varComputeVec.reserve(_lkpVarNames.size());

  // store the pointers to member functions for
  // computing physical quantities
  for (CFuint i = 0; i < _lkpVarNames.size(); ++i) {

    if (_lkpVarNames[i] == "e") {
      varComputeVec.push_back(&MutationLibrary::energy);
    }
    else if (_lkpVarNames[i] == "h") {
      varComputeVec.push_back(&MutationLibrary::enthalpy);
    }
    else if (_lkpVarNames[i] == "a") {
      varComputeVec.push_back(&MutationLibrary::soundSpeed);
    }
    else if (_lkpVarNames[i] == "d") {
      varComputeVec.push_back(&MutationLibrary::density);
    }
    else {
      throw Common::NoSuchValueException (FromHere(), "Variable name not found");
    }
  }

  // set the look up tables
  setTables(varComputeVec);
}

//////////////////////////////////////////////////////////////////////////////

void MutationLibrary::setTables(vector<ComputeQuantity>& varComputeVec)
{
  Common::Stopwatch<Common::WallTime> stp;
  stp.start();

  const CFuint nbT = static_cast<CFuint>
    ((_Tmax - _Tmin)/_deltaT);
  const CFuint nbP = static_cast<CFuint>
    ((_pmax - _pmin)/_deltaP);
  const CFuint nbLookUpVars = _lkpVarNames.size();

  // create mapping variable names with correspondind idxs
  _nameToIdxVar.reserve(nbLookUpVars);
  for (CFuint iVar = 0; iVar < nbLookUpVars; ++iVar) {
    _nameToIdxVar.insert(_lkpVarNames[iVar], iVar);
  }
  _nameToIdxVar.sortKeys();

  // compute the arrays containing all the keys for T and P
  vector<CFdouble> vT(nbT);
  vector<CFdouble> vP(nbP);
  for(CFuint i = 0; i < nbT; ++i) {
    vT[i] = _Tmin + i*_deltaT;
  }
  for(CFuint j = 0; j < nbP; ++j) {
    vP[j] =_pmin + j*_deltaP;
  }

  // initialize the table
  _lookUpTables.initialize(vT, vP, nbLookUpVars);

  // insert the values
  for(CFuint i = 0; i < nbT; ++i) {
    for(CFuint j = 0; j < nbP; ++j) {
      // set the composition at first
      setComposition(vT[i],vP[j], CFNULL);

      for (CFuint iVar = 0; iVar < nbLookUpVars; ++iVar) {
        CFdouble result = (this->*varComputeVec[iVar])(vT[i],vP[j]);
        _lookUpTables.insert(vT[i],vP[j],
                             _nameToIdxVar.find(_lkpVarNames[iVar]),
                             result);
      }
    }
  }

  CFLog(VERBOSE, "MutationLibrary::setTables() took " << stp << "s\n");
}

//////////////////////////////////////////////////////////////////////////////

void MutationLibrary::getSourceTermVT(CFdouble& temperature,
				      RealVector& tVec,
				      CFdouble& pressure,
				      CFdouble& rho,
				      RealVector& omegav,
	                              CFdouble& omegaRad)
{
  cf_assert(omegav.size() >= 1);

  CFdouble temp = max(_TminFix,temperature);
  CFdouble YTOL = 1e-20;
  CFdouble Predissociation = 0.3*8.314;

  for (CFuint i = 0; i < tVec.size(); ++i) {
   _TVARRAY[i] = tVec[i];
  }


  if (presenceElectron()) {
     FORTRAN_NAME(sourcevtions)(WR1, &LWR1, WR3, &LWR3, WI, &LWI,
   			 Y, &YTOL, &temp, _TVARRAY, &rho, &_nbTvib,
  			 &pressure, _STVIB, MOLARMASSP);

    for (CFuint i = 0; i < omegav.size(); ++i) {
     omegav[i] = _STVIB[i];
    }

    if (_nbTvib == 1) {
      //@TODO this is only valid for air5
/*      omegav[0] += Predissociation *(OMEGA[3]*113200./MOLARMASSP[3] +
				     OMEGA[4]*75500./MOLARMASSP[4] +
				     OMEGA[5]*59900./MOLARMASSP[5])+ OMEGA[6]*_HVIBR[6]/MOLARMASSP[6]
	+ OMEGA[7]*_HVIBR[7]/MOLARMASSP[7] + OMEGA[9]*_HVIBR[9]/MOLARMASSP[9];*/

    omegav[0] += OMEGA[3]*_HVIBR[3]/MOLARMASSP[3]+
    OMEGA[4]*_HVIBR[4]/MOLARMASSP[4]+OMEGA[5]*_HVIBR[5]/MOLARMASSP[5]+ OMEGA[6]*_HVIBR[6]/MOLARMASSP[6]
	+ OMEGA[7]*_HVIBR[7]/MOLARMASSP[7] + OMEGA[9]*_HVIBR[9]/MOLARMASSP[9];
/*           if (_TVARRAY[0] > 4000.0 ) {
       cout << "OMEGA " << OMEGA[0] <<endl;
       cout << "OMEGA " << OMEGA[1] <<endl;
       cout << "OMEGA " << OMEGA[2] <<endl;
       cout << "OMEGA " << OMEGA[3] <<endl;
       cout << "OMEGA " << OMEGA[4] <<endl;
       cout << "OMEGA " << OMEGA[5] <<endl;
       cout << "OMEGA " << OMEGA[6] <<endl;
       cout << "OMEGA " << OMEGA[7] <<endl;
       cout << "OMEGA " << OMEGA[8] <<endl;
       cout << "OMEGA " << OMEGA[9] <<endl;
       cout << "OMEGA " << OMEGA[10] <<endl;
       cout << "OMEGAVT " << omegav[0] <<endl;
       cout << "OMEGACV " <<OMEGA[3]*_HVIBR[3]/MOLARMASSP[3]+
    OMEGA[4]*_HVIBR[4]/MOLARMASSP[4]+OMEGA[5]*_HVIBR[5]/MOLARMASSP[5]+ OMEGA[6]*_HVIBR[6]/MOLARMASSP[6]
	+ OMEGA[7]*_HVIBR[7]/MOLARMASSP[7] + OMEGA[9]*_HVIBR[9]/MOLARMASSP[9] <<endl;
       abort();
    }*/

    }
    else {
      cf_assert(_nbTvib == 2);

      // Please remove it from here ASAP and do it more generally
      omegav[0] += OMEGA[3]*_HVIBR[3]/MOLARMASSP[3];  // This is for Nitrogen molecules
      omegav[1] += OMEGA[5]*_HVIBR[5]/MOLARMASSP[5];  // This is for Oxygen molecules
      // NOTICE THAT THE MOLARMASS IS NEEDED DUE TO THE FACT THAT
      // RT is per unit  mole.
      //      omegav[0] += Predissociation*OMEGA[3]*113200./MOLARMASSP[3];  // This is for Nitrogen molecules
      //      omegav[1] += Predissociation*OMEGA[5]*59900./MOLARMASSP[5];  // This is for Oxygen molecules
    }
  }
  else {
      //  TEMP TO BE FIXED !!!!
    if (_NS < 5) {
      FORTRAN_NAME(sourcevtn)(WR1, &LWR1, WR3, &LWR3, WI, &LWI,
			      Y, &YTOL, &temp, _TVARRAY, &rho,
			      &pressure, _STVIB, MOLARMASSP);

      omegav[0] = _STVIB[0];
      //@TODO this is only valid for NITROGEN2
      // omegav[0] += OMEGA[1]*_HVIBR[1]/MOLARMASSP[1];
      omegav[0] += Predissociation*(OMEGA[1]*113200./MOLARMASSP[1]);
    }
    else {
      FORTRAN_NAME(sourcevt)(WR1, &LWR1, WR3, &LWR3, WI, &LWI,
			     Y, &YTOL, &temp, _TVARRAY, &rho, &_nbTvib,
			     &pressure, _STVIB, MOLARMASSP);

      for (CFuint i = 0; i < omegav.size(); ++i) {
	omegav[i] = _STVIB[i];
      }

      if (_nbTvib == 1) {
	//@TODO this is only valid for air5
	omegav[0] += OMEGA[2]*_HVIBR[2]/MOLARMASSP[2]+
    OMEGA[3]*_HVIBR[3]/MOLARMASSP[3]+OMEGA[4]*_HVIBR[4]/MOLARMASSP[4];

	//      omegav[0] += Predissociation*(OMEGA[2]*113200./MOLARMASSP[2] +
	//      OMEGA[3]*75500./MOLARMASSP[3] +
	//      OMEGA[4]*59900./MOLARMASSP[4]);
      }
      else {
	cf_assert(_nbTvib == 2);
	//  omegav[0] += Predissociation*OMEGA[2]*113200./MOLARMASSP[2];  // This is for Nitrogen molecules
	omegav[0] += OMEGA[2]*_HVIBR[2]/MOLARMASSP[2];  // This is for Nitrogen molecules
	//omegav[1] += Predissociation*OMEGA[4]*59900./MOLARMASSP[4];  // This is for Oxygen molecules
	omegav[1] += OMEGA[4]*_HVIBR[4]/MOLARMASSP[4];  // This is for Oxygen molecules
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void MutationLibrary::getMolarMasses(RealVector& mm)
{
  cf_assert(mm.size() == static_cast<CFuint>(_NS));

  // check the units
  for (int i = 0; i < _NS; ++i) {
    mm[i] = MOLARMASSP[i];
  }
}


//////////////////////////////////////////////////////////////////////////////

void MutationLibrary::getGammaN(CFreal& m_GN)
{
  m_GN= _GammaN;
}

//////////////////////////////////////////////////////////////////////////////

void MutationLibrary::getGammaO(CFreal& m_GO)
{
  m_GO= _GammaO;
}

//////////////////////////////////////////////////////////////////////////////

void MutationLibrary::setMoleculesIDs(std::vector<CFuint>& v)
{
  if (_NS < 5) {
    v.resize(1);
    // @TODO AL: gory fix !!! only N2 work
      v[0] = 1;
  }
  else if (_NS >= 5) {
    if(presenceElectron()){
      v.resize(6);
      // @TODO MP: gory fix !!! only air 11 will work
      v[0] = 3;
      v[1] = 4;
      v[2] = 5;
      v[3] = 6;
      v[4] = 7;
      v[5] = 9;
    }
    else { // @TODO AL: gory fix !!! only air 5 will work
      v.resize(3);
      v[0] = 2;
      v[1] = 3;
      v[2] = 4;
    }
  }
}

//////////////////////////////////////////////////////////////////////////////



} // namespace Mutation

} // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

