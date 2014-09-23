#include "Mutation2.0.0I/MutationLibrary2OLD.hh"
#include "Mutation2.0.0I/Mutation2OLD.hh"
#include "Common/CFLog.hh"
#include "Environment/ObjectProvider.hh"
#include "Common/BadValueException.hh"
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
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace Mutation2OLD {

//////////////////////////////////////////////////////////////////////////////

extern "C" {

  /**
   * Importing of the FORTRAN subroutine LENGTHCF
   */
  void FORTRAN_NAME(lengthcf)(FINT, FINT, FINT, FINT, FINT,
                              FINT, FINT, FINT, FINT, FINT, FINT, FINT, FINT, FINT);

  /**
   * Importing of the FORTRAN subroutine MUTINITCF
   */
  void FORTRAN_NAME(mutinitcf)(FINT, FDOUBLE, FINT, FDOUBLE, FINT, FDOUBLE, FINT,
                             FDOUBLE, FINT, FINT, FINT, FCHAR, FINT, FINT, FINT,
           FINT, FINT, FINT, FINT);

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
  void FORTRAN_NAME(frozengammaneq)(FDOUBLE, FINT, FINT, FINT,
            FDOUBLE, FDOUBLE, FDOUBLE,
            FDOUBLE, FDOUBLE, FDOUBLE);

  /**
   * Importing of the FORTRAN subroutine FROZENGAMMA
   */
  void FORTRAN_NAME(frozengammafast)(FDOUBLE, FINT, FINT, FINT, FDOUBLE,
             FDOUBLE, FDOUBLE, FDOUBLE, FDOUBLE);

  /**
   * Importing of the FORTRAN subroutine FROZENGAMMANEQ
   */
  void FORTRAN_NAME(frozengammaneqfast)(FDOUBLE, FINT, FINT, FINT,
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
   * Importing of the FORTRAN subroutine _MOLARMASS
   */
  void FORTRAN_NAME(molarmass)(FDOUBLE ,FINT, FDOUBLE, FDOUBLE, FDOUBLE);

 /**
   * Importing of the FORTRAN subroutine ENTHALPY
   */
  void FORTRAN_NAME(enthalpy)(FDOUBLE, FINT, FINT ,FINT, FDOUBLE,
            FDOUBLE, FDOUBLE, FDOUBLE, FDOUBLE,
            FDOUBLE, FDOUBLE, FDOUBLE, FDOUBLE,
            FDOUBLE, FDOUBLE);
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
  void FORTRAN_NAME(lambdavib)(FDOUBLE, FINT, FDOUBLE, FINT, FINT,FINT,
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
   * Importing of the FORTRAN subroutine SETSPECIES_MOLARMASS
   */
  void FORTRAN_NAME(setspeciesmolarmass)(FDOUBLE, FINT, FDOUBLE);

  /**
   * Importing of the FORTRAN subroutine PRESSURE
   */
  void FORTRAN_NAME(pressure)(FDOUBLE, FINT, FDOUBLE, FDOUBLE,
            FDOUBLE, FDOUBLE, FDOUBLE);

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
  void FORTRAN_NAME(arrhenius)(FDOUBLE, FINT, FDOUBLE, FINT, FDOUBLE, FINT, 
                               FDOUBLE, FINT, FINT, FINT,
                               FDOUBLE, FDOUBLE, FDOUBLE, FDOUBLE, FDOUBLE,
                               FDOUBLE, FDOUBLE, FDOUBLE, FDOUBLE, FDOUBLE, 
			       FINT, FINT);
  /**
   * Importing of the FORTRAN subroutine ARRHENIUSRAD
   */
  void FORTRAN_NAME(arrheniusrad)(FDOUBLE,FINT, FDOUBLE, FINT, FINT, FINT, FDOUBLE, FDOUBLE, FDOUBLE,
                               FDOUBLE, FDOUBLE, FDOUBLE, FDOUBLE);

  /**
   * Importing of the FORTRAN subroutine OMEGAVTRANSFER
   */
  void FORTRAN_NAME(omegavtransfer)(FINT, FDOUBLE, FINT, FDOUBLE, FINT, FDOUBLE, FINT, FDOUBLE,
            FINT, FINT, FINT, FDOUBLE, FDOUBLE, FDOUBLE, FDOUBLE,
            FDOUBLE, FDOUBLE, FDOUBLE, FDOUBLE, FDOUBLE, FDOUBLE,
            FDOUBLE, FDOUBLE);

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
   * Importing of the FORTRAN subroutine GETCHARGE
   */
  void FORTRAN_NAME(getcharge)(FINT, FINT, FINT);

  /**
   * Importing of the FORTRAN subroutine RGAS
   */
  void FORTRAN_NAME(rgas)(FDOUBLE, FINT, FDOUBLE);

   /**
   * Importing of the FORTRAN subroutine LOADCR
   */
  void FORTRAN_NAME(loadcr)();
}

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<MutationLibrary2OLD,
			    PhysicalPropertyLibrary,
			    Mutation2OLDModule,
			    1>
mutation2OLDLibraryProvider("Mutation2OLD");

//////////////////////////////////////////////////////////////////////////////

void MutationLibrary2OLD::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< bool >
    ("noElectronicEnergy","Flag telling to ignore the electronic energy.");
  options.addConfigOption< bool >("includeElectronicEnergy","Include the electronic energy.");
  options.addConfigOption< CFint >("sonine","Sonine polynome order.");
  options.addConfigOption< CFint >("IMOD","Compute transport properties.");
  options.addConfigOption< std::string >("mixtureName","Name of the mixture.");
  options.addConfigOption< std::string >("reactionName","Name of the reaction.");
  options.addConfigOption< std::string >("transfName","Name of the transfer file.");
  options.addConfigOption< std::string >("thermoDir","Name of the directory where thermodynamic data are.");
  options.addConfigOption< CFdouble >("Xlim","Species molar fractions limit.");
  options.addConfigOption< std::string >("thermCondAlgo","Algo to compute thermal conductivity.");
  options.addConfigOption< std::string >("dynViscAlgo","Algo to compute dynamic viscosity.");
  options.addConfigOption< CFdouble >("TminFix","Minimum allowed temperature for chemistry.");
  options.addConfigOption< CFdouble >("GammaN","Factor of Catalycity of N.");
  options.addConfigOption< CFdouble >("GammaO","Factor of Catalycity of O.");
  options.addConfigOption< CFdouble >("Escape","Escape factor.");
  options.addConfigOption< CFint >("MoleculeID","Molecule ID for EV exchanges.");
  options.addConfigOption< CFuint >("CVModel","Model for the CV coupling.");
  options.addConfigOption< vector<CFreal> >
    ("PrefDissFactor","Preferential dissociation factors.");
  options.addConfigOption< vector<CFuint> >
    ("BoltzmannIDs","IDs of the species whose electronic energy has to be considered.");

  options.addConfigOption< CFuint >
    ("MolecTvID","ID to identify the molecule that gets all Tv contribution (ex. N2).");

  options.addConfigOption< CFdouble,
    Config::DynamicOption<> >("FactorOmega","Factor to reduce stiffness of chemical sorce terms.");
}

//////////////////////////////////////////////////////////////////////////////

MutationLibrary2OLD::MutationLibrary2OLD(const std::string& name)
  : PhysicalChemicalLibrary(name)
{
  addConfigOptionsTo(this);
  
  _mixtureName = "";
  setParameter("mixtureName",&_mixtureName);

  //"empty" for LTE, no finite rate chemistry
  _reactionName = "empty";
  setParameter("reactionName",&_reactionName);

  _transfName = "empty";
  setParameter("transfName",&_transfName);

  _thermoDir = "thermo";
  setParameter("thermoDir",&_thermoDir);

  _noElectEnergy = false;
  setParameter("noElectronicEnergy",&_noElectEnergy);

  _includeElectronicEnergy = false;
  setParameter("includeElectronicEnergy",&_includeElectronicEnergy);

  // IT DOESN t work with SONINE =2 The problem is somewhere inside
  // CORRECTION
  // _sonine = 2;
  _sonine = 1;
  setParameter("sonine",&_sonine);
  // 1: zero order, 2: first order in Mutation
  assert(_sonine < 3);

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

  _TminFix =100.;
  setParameter("TminFix",&_TminFix);

  _GammaN = 0.0;
  setParameter("GammaN",&_GammaN);
  
  _GammaO = 0.0;
  setParameter("GammaO",&_GammaO);
  
  _escape = 0.0;
  setParameter("Escape",&_escape);

  _iimol = -1;
  setParameter("MoleculeID",&_iimol);

  _cvModel = 0;
  setParameter("CVModel",&_cvModel);

  _prefDissFactor = vector<CFreal>();
  setParameter("PrefDissFactor",&_prefDissFactor);

  _boltzmannIDs = vector<CFuint>();
  setParameter("BoltzmannIDs",&_boltzmannIDs);

  _molecTvID = 2;
  setParameter("MolecTvID",&_molecTvID);

  _factorOmega = 1.0;
  setParameter("FactorOmega",&_factorOmega);

  _EPS = 1e-6; // this value gives problems with air11 (composition, speed of sound) and LTE
  _TOL = 1e-9; // recommended to replace with the use of Xlim and COMPOTOL2
}

//////////////////////////////////////////////////////////////////////////////

MutationLibrary2OLD::~MutationLibrary2OLD()
{
}

//////////////////////////////////////////////////////////////////////////////

void MutationLibrary2OLD::configure ( Config::ConfigArgs& args )
{
  PhysicalChemicalLibrary::configure(args);
  
  CFLog(NOTICE, "GammaN = " << _GammaN << "\n");
  CFLog(NOTICE, "GammaO = " << _GammaO << "\n");
  
  if (_lambdaAlgoStr == "CG") _lambdaAlgo = LAMBDACG;
  if (_lambdaAlgoStr == "Direct") _lambdaAlgo = LAMBDAD;
  if (_etaAlgoStr == "CG") _etaAlgo = ETACG;
  if (_etaAlgoStr == "Direct") _etaAlgo = ETAD;
}

//////////////////////////////////////////////////////////////////////////////

void MutationLibrary2OLD::setup()
{  
  // if this is a parallel simulation, only ONE process at a time
  // sets the library
  if (PE::GetPE().IsParallel()) {

    PE::GetPE().setBarrier();

    if (PE::GetPE().GetRank() == 0) {
       copyDataFiles();
     }
    
    PE::GetPE().setBarrier();

    for (CFuint i = 0; i < PE::GetPE().GetProcessorCount(); ++i) {

     if (i == PE::GetPE().GetRank ()) {
        setLibrarySequentially();
      }

      PE::GetPE().setBarrier();
    }
  
    if (PE::GetPE().GetRank() == 0) {
      deleteDataFiles();
    }
  }
  else {
    copyDataFiles();
    setLibrarySequentially();
    deleteDataFiles();
  }
}

//////////////////////////////////////////////////////////////////////////////

void MutationLibrary2OLD::copyDataFiles()
{
  // if nobody has configured the library path, it is done here
  // with default = basedir + "plugins/MutationI/data/mutation/"
  if (m_libPath == "") {
    CFLog(NOTICE, "MutationLibrary::libpath set to default" << "\n");
    std::string baseDir = Environment::DirPaths::getInstance().getBaseDir().string();
    m_libPath = baseDir + "/plugins/Mutation2.0.0I/";
  }
  
  cout << "m_libPath  =  " << m_libPath << endl;

  ofstream fout("mutation.in");
  fout << _mixtureName << endl;
  fout << _reactionName << endl;
  fout << _transfName << endl;
  fout << _thermoDir << endl;

  std::string command1 = "cp -R " + m_libPath + "data .";
  Common::OSystem::getInstance().executeCommand(command1);

  std::string command2 = "echo $PWD";
  Common::OSystem::getInstance().executeCommand(command2);

  CFout << "MutationLibrary::executing " << command1 << "\n";
  CFout << "MutationLibrary::libpath     = "  << m_libPath << "\n";
  CFout << "MutationLibrary::mixtureName = "  << _mixtureName << "\n";
  CFout << "MutationLibrary::reactionName = " << _reactionName << "\n";
  CFout << "MutationLibrary::transfName = " << _transfName << "\n";
  CFout << "MutationLibrary::thermoDir = " << _thermoDir << "\n";

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
}

//////////////////////////////////////////////////////////////////////////////

void MutationLibrary2OLD::deleteDataFiles()
{
  std::string command2 = "rm -fr data/ mutation.in";
  Common::OSystem::getInstance().executeCommand(command2);
}

//////////////////////////////////////////////////////////////////////////////

void MutationLibrary2OLD::setLibrarySequentially()
{
  // calculate the needed size for the work vectors
  FORTRAN_NAME(lengthcf)(&_LWR1,&_LWR2,&_LWR3,&_LWR4,&_LWI,&_LWC,&_NS,&_NE,&_NC,
                         &_NREA, &_NV, &_nbTvib, &_nbTe, &_IBINIJ);
  
  // _extraData vector inizialization
 /* _extraData.enthalpyTt.resize(_NS);
  _extraData.energyTr.resize(_NS);
  _extraData.energyVib.resize(_NS);
  _extraData.cpVib.resize(_nbTvib);
  _extraData.cpElec.resize(_NS);
  _extraData.eElec.resize(_NS);
  _extraData.enthalpyForm.resize(_NS);
  _extraData.dRhoEdRhoi.resize(_NS);
  _extraData.dRhoEvdRhoi.resize(_NS);  
  _atomicityCoeff.resize(_NS); */
 
  // allocate the private data
  _WR1 = new CFdouble[_LWR1];
  _WR2 = new CFdouble[_LWR2];
  _WR3 = new CFdouble[_LWR3];
  _WR4 = new CFdouble[_LWR4];
  _WI  = new CFint[_LWI];
  _WC  = new char[_LWC*4+1];
  _FIJ = new CFdouble[_NS*(_NS-1)/2];

  _X    =  new CFdouble[_NS];
  _XTOL =  new CFdouble[_NS];
  _Xini =  new CFdouble[_NS];
  _Xn   =  new CFdouble[_NC];
  _Yn   =  new CFdouble[_NC];
  _Y    =  new CFdouble[_NS];
  _Yini =  new CFdouble[_NS];

  // initialization
  for(CFint i = 0; i < _NS; ++i) {
    _X[i] = 1.0;
    _Y[i] = 1.0;
  }

  _LAMBDAEL =  new CFdouble[_NC];
  _ELDIFCOEF = new CFdouble[_NC*_NC];
  _ELTDIFCOEF = new CFdouble[_NC];
  _OMEGA = new CFdouble[_NS];
  _MOLARMASSP = new CFdouble[_NS];
  _DF =  new CFdouble[_NS];
  _JDIF =  new CFdouble[_NS];

  const CFuint sizeJ = _NS + _nbTvib + _nbTe + 1;
  _OMEGAJACOB = new CFdouble[_NS*sizeJ];
  _OMEGAVTJACOB = new CFdouble[(_nbTvib + _nbTe)*(sizeJ-1)];
  
  _HMTOTAL =  new CFdouble[_NS];
  _HMTRANS =  new CFdouble[_NS];
  _HMELECT =  new CFdouble[_NS];
  _HMROT   =  new CFdouble[_NS];
  _HMVIBR  =  new CFdouble[_NS];
  _HMFORM  =  new CFdouble[_NS];

  _HT  =  new CFdouble[_NS];
  _HTT =  new CFdouble[_NS];
  _HTE =  new CFdouble[_NS];
  _HTR =  new CFdouble[_NS];
  _HTV =  new CFdouble[_NS];
  _HTF =  new CFdouble[_NS];

  _HE  =  new CFdouble[_NS];
  _HET =  new CFdouble[_NS];
  _HEE =  new CFdouble[_NS];
  _HER =  new CFdouble[_NS];
  _HEV =  new CFdouble[_NS];
  _HEF =  new CFdouble[_NS];

  _HTOTAL =  new CFdouble[_NS];
  _HTRANS =  new CFdouble[_NS];
  _HELECT =  new CFdouble[_NS];
  _HROT   =  new CFdouble[_NS];
  _HVIBR  =  new CFdouble[_NS];
  _HFORM  =  new CFdouble[_NS];

  _HTOTALUM =  new CFdouble[_NS];
  _HTRANSUM =  new CFdouble[_NS];
  _HELECTUM =  new CFdouble[_NS];
  _HROTUM   =  new CFdouble[_NS];
  _HVIBRUM  =  new CFdouble[_NS];
  _HFORMUM  =  new CFdouble[_NS];

  _HTOTALP =  new CFdouble[_NS];
  _HTRANSP =  new CFdouble[_NS];
  _HELECTP =  new CFdouble[_NS];
  _HROTP   =  new CFdouble[_NS];
  _HVIBRP  =  new CFdouble[_NS];
  _HFORMP  =  new CFdouble[_NS];

  _HINT =  new CFdouble[_NS];

  _ETOTAL =  new CFdouble[_NS];
  _ETRANS =  new CFdouble[_NS];
  _EELECT =  new CFdouble[_NS];
  _EROT   =  new CFdouble[_NS];
  _EVIBR  =  new CFdouble[_NS];
  _EFORM  =  new CFdouble[_NS];

  _CPE =  new CFdouble[_NS];
  _CPR =  new CFdouble[_NS];
  _CPV =  new CFdouble[_NS];
  _CPINT = new CFdouble[_NS];
  _CPVIB = new CFdouble[_NS];

  _CHIH = new CFdouble[_NS];
  _CHARGE = new CFint[_NS];

  // initialize some data in work vectors
  // Careful!!! WR2 is filled by function "collision"
  FORTRAN_NAME(mutinitcf) (&_LWR1,_WR1,&_LWR2,_WR2,&_LWR3,_WR3,&_LWR4,_WR4,&_LWI,_WI,&_LWC,_WC,
         &_imod,&_IVIBTEMPI,&_IVIBSPEI,&_INVIBSPEI,&_IVIBI,&_IATOMI,&_IDIS);

  // calculate nuclear fraction ratios;
  FORTRAN_NAME(nuclear)(_WR1,&_LWR1,_Xn);

  const CFuint sizeTV = _nbTvib+1;
  _TVARRAY = new CFdouble[_NV];
  _STVIB   =  new CFdouble[sizeTV];
  _TVIBEPS1 = new CFdouble[_NV];
  _TVIBEPS  = new CFdouble[_NV];
  _LAMBDAVIB = new CFdouble[_nbTvib];
  _OMEGACV = new CFdouble[_nbTvib];
  _OMEGAVE = new CFdouble[_nbTvib];
  _OMEGAVV = new CFdouble[_nbTvib];
  _TVEC = new CFdouble[_nbTvib+2];
  _TT = new CFdouble[_NV];
  _TE = new CFdouble[_NV]; 

  _extraData.enthalpyTt.resize(_NS);
  _extraData.energyTr.resize(_NS);
  _extraData.energyVib.resize(_NS);
  _extraData.cpVib.resize(_nbTvib);
  _extraData.cpElec.resize(_NS);
  _extraData.eElec.resize(_NS);
  _extraData.enthalpyForm.resize(_NS);
  _extraData.dRhoEdRhoi.resize(_NS);
  _extraData.dRhoEvdRhoi.resize(_NS);
  _atomicityCoeff.resize(_NS);

  PhysicalChemicalLibrary::setup();

  // set the species molar masses
  FORTRAN_NAME(setspeciesmolarmass)(_WR1, &_LWR1, _MOLARMASSP);

  // get the charge of each species
  FORTRAN_NAME(getcharge)(_WI, &_LWI, _CHARGE);

  // constant of perfect gases
  FORTRAN_NAME(rgas)(_WR1, &_LWR1, &_Rgas);

  // index of N2 in the mixture (in FORTRAN style, i.e. starting from 1)
  // ask Marco ...This is just valid for air  or air11
  if (_iimol < 0) {
    _iimol = (presenceElectron()) ? 3 : 2;
  }
  
  setMoleculesIDs(_molecIDs);
  _flagMoleculesIDs.resize(_NS,false);
  for (CFuint i = 0; i < _molecIDs.size(); ++i) {
    _flagMoleculesIDs[_molecIDs[i]] = true;
  }
  for (CFint i = 0; i < _NS; ++i) {
    _atomicityCoeff[i] = (_flagMoleculesIDs[i]) ? 2.5 : 1.5;
  }
  
  // formation enthalpy per unity mass
  FORTRAN_NAME(enthalpyform)(_WR1, &_LWR1, _HFORM);
  for (CFuint i = 0; i < (CFuint) _NS; ++i) {
    _extraData.enthalpyForm[i] = _HFORM[i]/_MOLARMASSP[i];
  }

  _molecule2EqIDs.resize(_nbTvib);    // this vector now is no longer used
  if ((_NS == 5 || _NS == 11) && _nbTvib > 1) {
    _molecule2EqIDs[0] = 0;
    _molecule2EqIDs[1] = 1;
    _molecule2EqIDs[2] = 2;
}

  if (getNbTe() == 1) {
    _electrEnergyID = _nbTvib;
  }

  // Set by default all species IDs if none is specified for the
  // electronic energy to be considered
  if (_boltzmannIDs.size() == 0) {
    _boltzmannIDs.resize(_NS);
    for (CFuint i = 0; i < (CFuint) _NS; ++i) {
      _boltzmannIDs[i] = i;
    }
  }
  if (_boltzmannIDs.size() != static_cast<CFuint>(_NS)) {
    FORTRAN_NAME(loadcr)();
  }
  
  // set the flag telling if the mixture is ionized
  _hasElectrons = (_NE == 1) ? true : false;
  
   if (_mixtureName == "air5") {
   _molecTvID = 2;
   }
   if (_mixtureName == "air11") {
   _molecTvID = 3;
   }
}

//////////////////////////////////////////////////////////////////////////////

void MutationLibrary2OLD::unsetup()
{
  if(isSetup()) {

    deletePtrArray(_WR1);
    deletePtrArray(_WR2);
    deletePtrArray(_WR3);
    deletePtrArray(_WR4);
    deletePtrArray(_WC);
    deletePtrArray(_WI);
    deletePtrArray(_FIJ);

    deletePtrArray(_X);
    deletePtrArray(_XTOL);
    deletePtrArray(_Xini);
    deletePtrArray(_Xn);
    deletePtrArray(_Yn);
    deletePtrArray(_Y);
    deletePtrArray(_Yini);

    deletePtrArray(_LAMBDAEL);
    deletePtrArray(_ELDIFCOEF);
    deletePtrArray(_ELTDIFCOEF);
    deletePtrArray(_OMEGA);
    deletePtrArray(_MOLARMASSP);
    deletePtrArray(_DF);
    deletePtrArray(_JDIF);
    deletePtrArray(_OMEGAJACOB);

    deletePtrArray(_HMTOTAL);
    deletePtrArray(_HMTRANS);
    deletePtrArray(_HMELECT);
    deletePtrArray(_HMROT);
    deletePtrArray(_HMVIBR);
    deletePtrArray(_HMFORM);

    deletePtrArray(_HT);
    deletePtrArray(_HTT);
    deletePtrArray(_HTE);
    deletePtrArray(_HTR);
    deletePtrArray(_HTV);
    deletePtrArray(_HTF);

    deletePtrArray(_HE);
    deletePtrArray(_HET);
    deletePtrArray(_HEE);
    deletePtrArray(_HER);
    deletePtrArray(_HEV);
    deletePtrArray(_HEF);

    deletePtrArray(_HTOTAL);
    deletePtrArray(_HTRANS);
    deletePtrArray(_HELECT);
    deletePtrArray(_HROT);
    deletePtrArray(_HVIBR);
    deletePtrArray(_HFORM);

    deletePtrArray(_HTOTALUM);
    deletePtrArray(_HTRANSUM);
    deletePtrArray(_HELECTUM);
    deletePtrArray(_HROTUM);
    deletePtrArray(_HVIBRUM);
    deletePtrArray(_HFORMUM);

    deletePtrArray(_HTOTALP);
    deletePtrArray(_HTRANSP);
    deletePtrArray(_HELECTP);
    deletePtrArray(_HROTP);
    deletePtrArray(_HVIBRP);
    deletePtrArray(_HFORMP);

    deletePtrArray(_HINT);

    deletePtrArray(_ETOTAL);
    deletePtrArray(_ETRANS);
    deletePtrArray(_EELECT);
    deletePtrArray(_EROT);
    deletePtrArray(_EVIBR);
    deletePtrArray(_EFORM);

    deletePtrArray(_CPE);
    deletePtrArray(_CPR);
    deletePtrArray(_CPV);
    deletePtrArray(_CPINT);
    deletePtrArray(_CPVIB);
    deletePtrArray(_CHIH);
    deletePtrArray(_CHARGE);

    deletePtrArray(_TVARRAY);
    deletePtrArray(_STVIB);
    deletePtrArray(_TVIBEPS1);
    deletePtrArray(_TVIBEPS);
    deletePtrArray(_LAMBDAVIB);
    deletePtrArray(_OMEGACV);
    deletePtrArray(_OMEGAVE);
    deletePtrArray(_OMEGAVV);
    deletePtrArray(_TVEC);
    deletePtrArray(_TT);
    deletePtrArray(_TE);

    PhysicalChemicalLibrary::unsetup();
  }
}

//////////////////////////////////////////////////////////////////////////////

//thermal conductivity by direct method
CFdouble MutationLibrary2OLD::lambdaD(CFdouble& temperature,
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
    //cout << " MutationLibrary2OLD::lambdaD() => Watch Out !! T <<<<< 100. " << temp << endl;
    temp = 100.;
  }

  CFdouble ND = 0.0;
  FORTRAN_NAME(numberd)(_WR1,&_LWR1,&pressure,&temp,&temp,_X,&ND); //ND needed for collision below

  //Thermodynamic properties
  setDefaultTVarray(temp);

  // - Species specific heat per unit mole
  FORTRAN_NAME(enthalpy)(_WR1, &_LWR1, _WI, &_LWI, &temp, &temp, &temp, _TVARRAY, &pressure,
       _HTOTAL, _HTRANS, _HELECT, _HROT, _HVIBR, _HFORM);

  CFdouble tempEps1 = temp*(1.0+_EPS);
  CFdouble tempEps = temp*_EPS;

  setDefaultTVarray(tempEps1);

  FORTRAN_NAME(enthalpy)(_WR1, &_LWR1, _WI, &_LWI, &tempEps1, &tempEps1, &tempEps1, _TVARRAY,
       &pressure, _HTOTALP, _HTRANSP, _HELECTP, _HROTP, _HVIBRP, _HFORMP);

  for(int is=0 ; is < _NS; is++) {
    _CPE[is]   = (_HELECTP[is] - _HELECT[is]) / (tempEps);
    _CPR[is]   = (_HROTP[is] - _HROT[is]) / (tempEps);
    _CPV[is]   = (_HVIBRP[is] - _HVIBR[is]) / (tempEps);
    _CPINT[is] = _CPE[is] + _CPR[is] + _CPV[is];
  }

  FORTRAN_NAME(compotol)(_X, &_TOL, _XTOL);
  //Update of kinetic data
  FORTRAN_NAME(collision)(_WR1, &_LWR1, _WR2, &_LWR2, &temp, &temp, &ND, _X);

  //A. Heavy particles properties

  //Eucken corrections for internal energy (no inelastic collisions)
  FORTRAN_NAME(lambdaint)(_WR1, &_LWR1, _WR2, &_LWR2, _CPINT, _XTOL, &lambdaInth);
  FORTRAN_NAME(lambdachid)(_WR1, &_LWR1, _WR2, &_LWR2, _XTOL, &lambdaTh, _CHIH);
  FORTRAN_NAME(lambdaread)(_WI, &_LWI, _WR1, &_LWR1, _WR2, &_LWR2, _HTOTAL,
         &temp, _XTOL, &lambdaRea);

  //B. Electron gas properties
  if(presenceElectron()){
    int i=3;
    FORTRAN_NAME(lambdae)(_WR1, &_LWR1, _WR2, &_LWR2, _X, &temp, &lambdaTe, &i);
  }

  lambdaTot = lambdaInth + lambdaTh + lambdaRea + lambdaTe;

  return lambdaTot;
}

//////////////////////////////////////////////////////////////////////////////

//dynamic viscosity by direct method
CFdouble MutationLibrary2OLD::etaD(CFdouble& temperature,
        CFdouble& pressure,
        CFreal* tVec)
{
  //composition must be called before!
  CFdouble eta = 0.0;
  CFdouble ND = 0.0;
  CFdouble temp = temperature;
  if (temp < 100.) {
    //cout << " MutationLibrary2OLD::etaD() => Watch Out !! T <<<<< 100. " << temp << endl;
    temp = 100.;
  }

  CFdouble Te = getTe(temp,tVec);
  FORTRAN_NAME(numberd)(_WR1,&_LWR1,&pressure,&temp,&Te,_X,&ND);

  // FORTRAN_NAME(compotol2)(_X, &_Xlim, _XTOL);
  FORTRAN_NAME(compotol)(_X, &_TOL, _XTOL);
  //Update of kinetic data
  FORTRAN_NAME(collision)(_WR1, &_LWR1, _WR2, &_LWR2, &temp, &Te, &ND, _X);

  //A. Heavy particles properties
  //Eucken corrections for internal energy (no inelastic collisions)
  FORTRAN_NAME(etad)(_WR1, &_LWR1, _WR2, &_LWR2, _XTOL, &eta);

  return eta;
}

//////////////////////////////////////////////////////////////////////////////

//thermal conductivity by conjugate gradient method
CFdouble MutationLibrary2OLD::lambdaCG(CFdouble& temperature,
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
    //cout << " MutationLibrary2OLD::lambdaCG() => Watch Out !! T <<<<< 100. " << temp << endl;
    temp = 100.;
  }

  FORTRAN_NAME(numberd)(_WR1,&_LWR1,&pressure,&temp,&temp,_X,&ND); //ND needed for collision below

  //Thermodynamic properties

  // - Species specific heat per unit mole
  setDefaultTVarray(temp);

  FORTRAN_NAME(enthalpy)(_WR1, &_LWR1, _WI, &_LWI, &temp, &temp, &temp, _TVARRAY, &pressure,
       _HTOTAL, _HTRANS, _HELECT, _HROT, _HVIBR, _HFORM);

  CFdouble tempEps1 = temp*(1.0+_EPS);
  CFdouble tempEps = temp*_EPS;

  setDefaultTVarray(tempEps1);

  FORTRAN_NAME(enthalpy)(_WR1, &_LWR1, _WI, &_LWI, &tempEps1, &tempEps1, &tempEps1, _TVARRAY,
       &pressure, _HTOTALP, _HTRANSP, _HELECTP, _HROTP, _HVIBRP, _HFORMP);

  for(int is=0 ; is < _NS; is++) {
    _CPE[is]   = (_HELECTP[is] - _HELECT[is]) / (tempEps);
    _CPR[is]   = (_HROTP[is] - _HROT[is]) / (tempEps);
    _CPV[is]   = (_HVIBRP[is] - _HVIBR[is]) / (tempEps);
    _CPINT[is] = _CPE[is] + _CPR[is] + _CPV[is];
  }

  FORTRAN_NAME(compotol)(_X, &_TOL, _XTOL);
  //Update of kinetic data
  FORTRAN_NAME(collision)(_WR1, &_LWR1, _WR2, &_LWR2, &temp, &temp, &ND, _X);

  //A. Heavy particles properties

  //Eucken corrections for internal energy (no inelastic collisions)

  FORTRAN_NAME(lambdaint)(_WR1, &_LWR1, _WR2, &_LWR2, _CPINT, _XTOL, &lambdaInth);
  FORTRAN_NAME(lambdachicg)(_WR1, &_LWR1, _WR2, &_LWR2, _XTOL, &lambdaTh, _CHIH);
  FORTRAN_NAME(lambdareacg)(_WI, &_LWI, _WR1, &_LWR1, _WR2, &_LWR2, _HTOTAL,
          &temp, _XTOL, &lambdaRea);

  //B. Electron gas properties
  if(presenceElectron()){
    int i=3;
    FORTRAN_NAME(lambdae)(_WR1, &_LWR1, _WR2, &_LWR2, _X, &temp, &lambdaTe, &i);
  }

  return lambdaInth + lambdaTh + lambdaRea + lambdaTe;

}

//////////////////////////////////////////////////////////////////////////////

CFdouble MutationLibrary2OLD::lambdaNEQ(CFdouble& temperature,
             CFdouble& pressure)
{
  //composition must be called before!
  CFdouble lambdaTot = 0.0;
  CFdouble lambdaInth = 0.0;
  CFdouble lambdaTh = 0.0;
  CFdouble lambdaTe = 0.0;
  CFdouble temp = temperature;

  if (temp < 100.) {
    //    cout << " MutationLibrary2OLD::lambdaNEQ() => Watch Out !! T <<<<< 100. " << temp << endl;
    temp = 100.;
  }

  CFdouble ND = 0.0;
  FORTRAN_NAME(numberd)(_WR1,&_LWR1,&pressure,&temp,&temp,_X,&ND); //ND needed for collision below
  //Thermodynamic properties

  // - Species specific heat per unit mole
  setDefaultTVarray(temp);

  FORTRAN_NAME(enthalpy)(_WR1, &_LWR1, _WI, &_LWI, &temp, &temp, &temp, _TVARRAY, &pressure,
       _HTOTAL, _HTRANS, _HELECT, _HROT, _HVIBR, _HFORM);
  CFdouble tempEps1 = temp*(1.0+_EPS);
  CFdouble tempEps = temp*_EPS;

  setDefaultTVarray(tempEps1);

  FORTRAN_NAME(enthalpy)(_WR1, &_LWR1, _WI, &_LWI, &tempEps1, &tempEps1, &tempEps1, _TVARRAY,
       &pressure, _HTOTALP, _HTRANSP, _HELECTP, _HROTP, _HVIBRP, _HFORMP);

  for(int is=0 ; is < _NS ; is++) {
    _CPE[is]   = (_HELECTP[is] - _HELECT[is]) / (tempEps);
    _CPR[is]   = (_HROTP[is] - _HROT[is]) / (tempEps);
    _CPV[is]   = (_HVIBRP[is] - _HVIBR[is]) / (tempEps);
    _CPINT[is] = _CPE[is] + _CPR[is] + _CPV[is];
  }

  FORTRAN_NAME(compotol)(_X, &_TOL, _XTOL);
  //Update of kinetic data
  FORTRAN_NAME(collision)(_WR1, &_LWR1, _WR2, &_LWR2, &temp, &temp, &ND, _X);

  //A. Heavy particles properties

  //Eucken corrections for internal energy (no inelastic collisions)
  FORTRAN_NAME(lambdaint)(_WR1, &_LWR1, _WR2, &_LWR2, _CPINT, _XTOL, &lambdaInth);

  if (_lambdaAlgo == LAMBDACG) {
    FORTRAN_NAME(lambdachicg)(_WR1, &_LWR1, _WR2, &_LWR2, _XTOL, &lambdaTh, _CHIH);
  }
  else {
    FORTRAN_NAME(lambdachid)(_WR1, &_LWR1, _WR2, &_LWR2, _XTOL, &lambdaTh, _CHIH);
   }

  //B. Electron gas properties
  if(presenceElectron()){
    int i=3;
    FORTRAN_NAME(lambdae)(_WR1, &_LWR1, _WR2, &_LWR2, _X, &temp, &lambdaTe, &i);
  }

  lambdaTot = lambdaInth + lambdaTh + lambdaTe;

  return lambdaTot;
}

//////////////////////////////////////////////////////////////////////////////

void MutationLibrary2OLD::lambdaVibNEQ(CFreal& temperature,
            RealVector& tVec,
            CFdouble& pressure,
            CFreal& lambdaTrRo,
            RealVector& lambdaInt)
{
  CFdouble lambdaRot = 0.0;
  CFdouble lambdaTh  = 0.0;
  CFdouble lambdaE   = 0.0;
  CFdouble lambdaTe  = 0.0;
  CFdouble temp = temperature;

  if (temp < 100.) {
    //cout << "MutationLibrary2OLD::lambdaNEQ() => Watch Out !! T <<<<< 100. " << temp << endl;
    temp = 100.;
  }

  CFuint ic = 0;
  for (CFint i = 0; i < _NS; ++i) {
    const CFuint nvibMode = _WI[_IVIBI+i-1];
    for (CFuint j = 0; j < nvibMode; ++j, ++ic) {
      const CFuint ivib = _WI[_IVIBTEMPI+i-1]-2;
      cf_assert(ivib < tVec.size());
      _TVARRAY[ic] = tVec[ivib];
      _TVIBEPS1[ic] = tVec[ivib] *(1.0+_EPS);
      _TVIBEPS[ic] = tVec[ivib]*_EPS;
    }
  }

  const CFuint TeID = tVec.size()-1;
  CFreal Te = (_nbTe == 0)? tVec[0] : tVec[TeID];


  CFdouble TeEps = Te*_EPS;
  CFdouble TeEps1 = Te*(1.0+_EPS);
  CFdouble tempEps1 = temp*(1.0+_EPS);
  CFdouble tempEps = temp*_EPS;
  CFdouble ND = 0.0;
  FORTRAN_NAME(numberd)(_WR1,&_LWR1,&pressure,&temp,&Te,_X,&ND); //ND needed for collision below

  //Thermodynamic properties
  // - Species specific heat per unit mole
  FORTRAN_NAME(enthalpy)(_WR1, &_LWR1, _WI, &_LWI, &temp, &Te, &temp, _TVARRAY, &pressure,
       _HTOTAL, _HTRANS, _HELECT, _HROT, _HVIBR, _HFORM);

  FORTRAN_NAME(enthalpy)(_WR1, &_LWR1, _WI, &_LWI, &tempEps1, &TeEps1, &tempEps1, _TVIBEPS1,
       &pressure, _HTOTALP, _HTRANSP, _HELECTP, _HROTP, _HVIBRP, _HFORMP);

  for(CFint is = 0; is < _NS ; is++) {
    // _CPE[is]   = (_HELECTP[is] - _HELECT[is])/TeEps;
    _CPR[is]   = (_HROTP[is] - _HROT[is])/tempEps;

    // to generalize this !!!
    if (_WI[_IVIBI+is-1] == 1) {
      _CPVIB[is] = (_HVIBRP[is] - _HVIBR[is])/(tVec[_WI[_IVIBTEMPI+is-1]-2] *_EPS);
    }
    else {
      _CPVIB[is] = 0.0;
    }
  }
  
  // to test ...............
  const CFuint bsize = _boltzmannIDs.size();
  for(CFuint is = 0; is < bsize; is++) {
    const CFuint elEnergyID = _boltzmannIDs[is];
    _CPE[elEnergyID]   = (_HELECTP[elEnergyID] - _HELECT[elEnergyID])/TeEps;
  }

  FORTRAN_NAME(compotol)(_X, &_TOL, _XTOL);
  //Update of kinetic data
  FORTRAN_NAME(collision)(_WR1, &_LWR1, _WR2, &_LWR2, &temp, &Te, &ND, _X);

  //A. Heavy particles properties
  //Eucken corrections for rotational internal energy (no inelastic collisions)
  FORTRAN_NAME(lambdaint)(_WR1, &_LWR1, _WR2, &_LWR2, _CPR, _XTOL, &lambdaRot);

  FORTRAN_NAME(lambdaint)(_WR1, &_LWR1, _WR2, &_LWR2, _CPE, _XTOL, &lambdaE);

  FORTRAN_NAME(lambdavib)(_WR1, &_LWR1, _WR2, &_LWR2, _WI, &_LWI, _CPVIB, _XTOL,_LAMBDAVIB);

  if (_lambdaAlgo == LAMBDACG) {
    FORTRAN_NAME(lambdachicg)(_WR1, &_LWR1, _WR2, &_LWR2, _XTOL, &lambdaTh, _CHIH);
  }
  else {
    FORTRAN_NAME(lambdachid)(_WR1, &_LWR1, _WR2, &_LWR2, _XTOL, &lambdaTh, _CHIH);
  }

  //B. Electron gas properties
  int i=3;
  FORTRAN_NAME(lambdae)(_WR1, &_LWR1, _WR2, &_LWR2, _XTOL, &Te, &lambdaTe, &i);
  lambdaTrRo = lambdaRot + lambdaTh;

  for (CFint i = 0; i < _nbTvib; i++ ){
    lambdaInt[i] = _LAMBDAVIB[i];
  }

  if(_nbTe == 1) {
    lambdaInt[_nbTvib] = lambdaTe + lambdaE;
  }
  else {
    if (presenceElectron()) {
      lambdaInt[max(0,_electrEnergyID)] +=  lambdaE + lambdaTe;
    }
    else {
      if (_includeElectronicEnergy) {
  lambdaInt[max(0,_electrEnergyID)] +=  lambdaE;
      }
      lambdaTe = 0.0;
    }
  }
  /*  if (Te > 4000) {
  cout <<"lambdaTrRo" << lambdaTrRo<<endl;
  cout <<"lambdaVIB" << _LAMBDAVIB[0]<<endl;
//  abort();
  } */

}

//////////////////////////////////////////////////////////////////////////////

CFdouble MutationLibrary2OLD::etaCG(CFdouble& temperature,
         CFdouble& pressure,
         CFreal* tVec)
{
  //composition must be called before!
  CFdouble eta = 0.0;
  CFdouble ND = 0.0;
  CFdouble temp = temperature;
  if (temp < 100.) {
    // cout << " MutationLibrary2OLD::etaCG() => Watch Out !! T <<<<< 100. " << temp << endl;
    temp = 100.;
  }

  CFdouble Te = getTe(temp,tVec);
  FORTRAN_NAME(numberd)(_WR1,&_LWR1,&pressure,&temp,&Te,_X,&ND); //ND needed for collision below

  FORTRAN_NAME(compotol)(_X, &_TOL, _XTOL);

  //Update of kinetic data
  FORTRAN_NAME(collision)(_WR1, &_LWR1, _WR2, &_LWR2, &temp, &Te, &ND, _X);

  //A. Heavy particles properties
  //Eucken corrections for internal energy (no inelastic collisions)
  FORTRAN_NAME(etacg)(_WR1, &_LWR1, _WR2, &_LWR2, _XTOL, &eta);

  return eta;
}

//////////////////////////////////////////////////////////////////////////////

CFdouble MutationLibrary2OLD::sigma(CFdouble& temp,   //electrical conductivity
				    CFdouble& pressure,
				    CFreal* tVec)
{
//composition must be called before!
  CFdouble sigma = 0.0;
  CFdouble ND = 0.0;
  CFdouble Te = getTe(temp,tVec);
  
  FORTRAN_NAME(numberd)(_WR1,&_LWR1,&pressure,&temp,&Te,_X,&ND); //ND needed for collision below
  
  FORTRAN_NAME(compotol)(_X, &_TOL, _XTOL);
  //Update of kinetic data
  FORTRAN_NAME(collision)(_WR1, &_LWR1, _WR2, &_LWR2, &temp, &Te, &ND,_X);
  
  int i=2;
  FORTRAN_NAME(sigmae)(_WR1, &_LWR1, _WR2, &_LWR2, _XTOL, &temp, &sigma, &i);

  return sigma;
}

//////////////////////////////////////////////////////////////////////////////

void MutationLibrary2OLD::gammaAndSoundSpeed(CFdouble& temp,
           CFdouble& pressure,
           CFdouble& rho,
           CFdouble& gamma,
           CFdouble& soundSpeed)
{
  // compute the ratio of the mixture frozen specific heat in thermal equilibrium
  CFdouble drhodp = 0.0;
  CFdouble eps = 0.1;
  FORTRAN_NAME(equigamma)(_WR1, &_LWR1, _WI, &_LWI, &temp, &pressure,
        &rho, _Xn, _X, &eps, &gamma, &drhodp);
  soundSpeed = sqrt(gamma/drhodp);
}

//////////////////////////////////////////////////////////////////////////////

void MutationLibrary2OLD::frozenGammaAndSoundSpeed(CFdouble& temp,
                 CFdouble& pressure,
                 CFdouble& rho,
                 CFdouble& gamma,
                 CFdouble& soundSpeed,
                 RealVector* tVec)
{
  // compute the ratio of the mixture frozen specific heat in thermal equilibrium
  // if (tVec == CFNULL) {
//     FORTRAN_NAME(frozengammafast)(_WR1, &_LWR1, _WI, &_LWI, &temp, &pressure,
//    				  _X, &_EPS, &gamma);
//   }
//   else { 
//      setTVarray(*tVec);
//      CFdouble Te = getTe(temp,&(*tVec)[0]);
  
//      FORTRAN_NAME(frozengammaneqfast)(_WR1, &_LWR1, _WI, &_LWI, &temp,
//  				     &Te, _TVARRAY, &pressure, _X, &_EPS, &gamma);
//    }
//    soundSpeed = sqrt(gamma*pressure/rho); 
  
  // cout << "(1) a,g = " <<  soundSpeed << ", "  << gamma << endl; 
  
  if (getNbTe() == 0) {
    CFreal numBeta = 0.;
    CFreal denBeta = 0.;
    const CFint start = (presenceElectron()) ? 1 : 0;
    for (CFint i = start; i < _NS; ++i) {
      const CFreal sigmai = _Y[i]/_MOLARMASSP[i];
      numBeta += sigmai;
      denBeta += sigmai*_atomicityCoeff[i];
    }
    
    gamma = 1 + numBeta/denBeta;
  }
  else {
    setTVarray(*tVec);
    CFdouble Te = getTe(temp,&(*tVec)[0]);
    
    FORTRAN_NAME(frozengammaneqfast)(_WR1, &_LWR1, _WI, &_LWI, &temp,
				     &Te, _TVARRAY, &pressure, _X, &_EPS, &gamma);
  }
  
  soundSpeed = std::sqrt(gamma*pressure/rho); 
  
  //cout << "(2) a,g = " <<  soundSpeed << ", "  << gamma << endl; 
}

//////////////////////////////////////////////////////////////////////////////

CFdouble MutationLibrary2OLD::soundSpeed(CFdouble& temp,
                                      CFdouble& pressure)
{
  // compute the ratio of the mixture frozen specific heat in thermal equilibrium  c
  CFdouble gamma = 0.0;
  CFdouble drhodp = 0.0;
  CFdouble rho = density(temp, pressure, CFNULL);
  CFdouble eps = 0.1;
  FORTRAN_NAME(equigamma)(_WR1, &_LWR1, _WI, &_LWI, &temp, &pressure,
        &rho, _Xn, _X, &eps, &gamma, &drhodp);
  return sqrt(gamma/drhodp);
}

//////////////////////////////////////////////////////////////////////////////

void MutationLibrary2OLD::setComposition(CFdouble& temp,
					 CFdouble& pressure,
					 RealVector* x)
{
  // initialization of the molar fractions
  for(CFint i = 0; i < _NS; ++i) {
    _Xini[i] = 1.0;
  }
  
  // compute the molar fractions corresponding to the given temperature and
  // pressure for the mixture Xn
  FORTRAN_NAME(composition)(_WR1,&_LWR1,_WI,&_LWI,&temp,&pressure,_Xn,_Xini,_X);
  
  if (x != CFNULL) {
    for(CFint i = 0; i < _NS; ++i) {
      (*x)[i] = static_cast<CFreal>(_X[i]);
    }
  } 
  
  // set mass fractions which will be used later
  CFreal massTot = 0.;
  for (CFint is = 0; is < _NS; ++is) {
    if (_X[is] > 1.00000000001) {
      cout << "X[" << is << "] = " << _X[is] << endl;
      abort();
    }
    const CFreal mm = _X[is]*_MOLARMASSP[is];
    massTot += mm;
    _Y[is] = mm;
  }
  
  for (CFint is = 0; is < _NS; ++is) {
    _Y[is] /= massTot;
  }
}

//////////////////////////////////////////////////////////////////////////////

void MutationLibrary2OLD::setDensityEnthalpyEnergy(CFdouble& temp,
						   CFdouble& pressure,
						   RealVector& dhe)
{
  //setComposition must be called before
  CFdouble rho = 0.0;

  setDefaultTVarray(temp);

  FORTRAN_NAME(enthalpymass)(_WR1,&_LWR1,_WI,&_LWI,&temp,&temp,&temp,_TVARRAY,
       &pressure,_HTOTALUM,_HTRANSUM,_HELECTUM,_HROTUM, _HVIBRUM,_HFORMUM);

  dhe = 0.0;
  
  if (_noElectEnergy) {
    for(CFint i = 0; i < _NS; ++i) {
      dhe[1] += _Y[i]*(_HTOTALUM[i] - _HELECTUM[i]);
      rho += _X[i]*pressure/((_Rgas/_MOLARMASSP[i])*temp);
    }
  } 
  else {
    for(CFint i = 0; i < _NS; ++i) {
      dhe[1] += _Y[i]*_HTOTALUM[i];
      rho += _X[i]*pressure/((_Rgas/_MOLARMASSP[i])*temp);
    }
  }
  
  dhe[0] = rho;
  dhe[2] = dhe[1] - pressure/rho; 
}

//////////////////////////////////////////////////////////////////////////////

void MutationLibrary2OLD::setDensityEnthalpyEnergy(CFdouble& temp,
            RealVector& tVec,
            CFdouble& pressure,
            RealVector& dhe,
            bool storeExtraData)
{
  CFdouble rho = 0.0;
  CFreal Te = 0.0;
  
  if (_nbTvib != 0) {
    assert(tVec.size() == static_cast<CFuint>(_nbTvib + _nbTe)); 
    const CFuint TeID = tVec.size()-1;
    Te = (_nbTe == 0)? tVec[0] : tVec[TeID];
    setTVarray(tVec);
  } 
  else {
    Te = temp;
    setDefaultTVarray(temp);    
  }
  
  // set the working _TVEC array
  setTVec(temp,tVec);
  
  FORTRAN_NAME(enthalpymass)(_WR1,&_LWR1,_WI,&_LWI, &temp, &Te, &temp, _TVARRAY,
			     &pressure, _HTOTALUM,_HTRANSUM,_HELECTUM,_HROTUM,
			     _HVIBRUM,_HFORMUM);  
  
  // Checked OK Compared with MUTATION, Checked HTOT
  // sum up all the internal energies for each species
  dhe = 0.0;
  
  for(CFint i = 0; i < _NS; ++i) {  
    dhe[1] += _Y[i]*_HTOTALUM[i];
    if (presenceElectron() && i == 0) {
      rho += _X[i]*pressure/((_Rgas/_MOLARMASSP[i])*Te);
    }  
    else {
      rho += _X[i]*pressure/((_Rgas/_MOLARMASSP[i])*temp); 
    }
    // to be tested
    if (_noElectEnergy) {
      dhe[1] -= _Y[i]*_HELECTUM[i];
    }
  } 
  
  CFuint icc=0;
  for (CFint i = 0; i < _nbTvib; ++i) {
    const CFuint nvSize = _WI[_INVIBSPEI+i-1];
    for (CFuint j = 0; j < nvSize; ++j, ++icc) {
      const CFuint idx = _WI[_IVIBSPEI+icc-1]-1; 
      dhe[3+i] += _Y[idx]*_HVIBRUM[idx]; 
    }
  }
  
  //      FOR CR IT WE NEED TO SPLIT THAT WAY !!
  if (getNbTe() == 0 && presenceElectron()) {
    //    TO BE CHANGED MARCO
    //    dhe[3+max(0,_electrEnergyID)] +=_Y[0]*_HTRANSUM[0];
    if (_includeElectronicEnergy) {
      dhe[3+max(0,_electrEnergyID)] +=_Y[0]*_HTRANSUM[0]; 
      for (CFuint i = 0; i < _boltzmannIDs.size(); ++i) {
            const CFuint elEnergyID = _boltzmannIDs[i];
            dhe[3+max(0,_electrEnergyID)] += _Y[elEnergyID]*_HELECTUM[elEnergyID];
      }
    }
  } 
  
  // Checked _WI[_IVIBSPEI+icc-1] vector
  // electronic energy
  if (getNbTe() == 1) {
    cf_assert(static_cast<CFuint>(3+_nbTvib) < dhe.size()); 
    dhe[3+_nbTvib] = _Y[0]*_HTRANSUM[0];
    if (_includeElectronicEnergy) {
      for (CFuint i = 0; i < _boltzmannIDs.size(); ++i) {
	const CFuint elEnergyID = _boltzmannIDs[i];
	dhe[3+_nbTvib] += _Y[elEnergyID]*_HELECTUM[elEnergyID];
      }
    }
  }
  
  dhe[0] = rho;
  dhe[2] = dhe[1] - pressure/rho;
  
  if (storeExtraData) {
    _extraData.dEdT = 0.0;
    _extraData.dHdT = 0.0;
    
    for(CFint i = 0; i < _NS; ++i) {
      const CFreal Ti = (presenceElectron() && (i == 0)) ? Te : temp; 
      _extraData.energyTr[i] = _HTRANSUM[i] + _HROTUM[i] + _HFORMUM[i] - (_Rgas/_MOLARMASSP[i])*Ti;
      _extraData.enthalpyTt[i] = (_noElectEnergy) ? _HTOTALUM[i] - _HELECTUM[i] : _HTOTALUM[i];
      _extraData.dRhoEdRhoi[i] = _extraData.energyTr[i];
      _extraData.energyVib[i] = _HVIBRUM[i];
      
      if (_WI[_IATOMI+i-1] == 0) {      
        _extraData.eElec[i] = _HTOTALUM[i];    //   free electrons 
      } else {
        _extraData.eElec[i] = _HELECTUM[i];    //   other species
      } 
      
      const CFreal yReli = _Y[i]*_Rgas/_MOLARMASSP[i]; 
      if (_flagMoleculesIDs[i]) {
	_extraData.dRhoEvdRhoi[i] = _HVIBRUM[i];     
	_extraData.dHdT += 3.5*yReli;
	_extraData.dEdT += 2.5*yReli;
      }
      else {
	_extraData.dHdT += 2.5*yReli;
	_extraData.dEdT += 1.5*yReli;
      }
    }
    
    // Thermochemical nonequilibrium 
    if (_nbTvib > 0) {
      if (_nbTvib == 1) {
	CFdouble startTv = _TVARRAY[0]; 
	CFdouble epsTv = startTv*_EPS;
	_TVARRAY[0] = startTv + epsTv;
	CFdouble TePert = Te*(1. + _EPS);
	
	FORTRAN_NAME(enthalpymass)
	  (_WR1,&_LWR1,_WI,&_LWI, &temp, &TePert, &temp, _TVARRAY,
	   &pressure, _HTOTALUM,_HTRANSUM,_HELECTUM,_HROTUM, _HVIBRUM,_HFORMUM);
	
        _extraData.cpVib = 0.0;
        _extraData.dEvTv = 0.0;
	for(CFint is = 0; is < _NS; ++is) {
	  _extraData.cpVib[_nbTvib-1] += _Y[is]*(_HVIBRUM[is] - _extraData.energyVib[is]);   
	}
        _extraData.cpVib[_nbTvib-1] /= epsTv; 
        _extraData.dEvTv = _extraData.cpVib[_nbTvib-1];  
      }
      else {
	CFdouble TePert = Te*(1.0+_EPS);
	for (CFint v = 0; v < _NV; v++) {
	  _TVARRAY[v] += _TVARRAY[v]*_EPS;      
	}  
	
	FORTRAN_NAME(enthalpymass)
	  (_WR1,&_LWR1,_WI,&_LWI, &temp, &TePert, &temp, _TVARRAY,
	   &pressure, _HTOTALUM,_HTRANSUM,_HELECTUM,_HROTUM, _HVIBRUM,_HFORMUM); 
	
	CFuint icc = 0;
	_extraData.list.clear(); // there was a mistake here! you cannot push_back continuously !!!  
	_extraData.cpVib = 0.;
	
	for (CFint v = 0; v < _nbTvib; v++) {
	  _extraData.list.push_back(_WI[_INVIBSPEI+v-1]);
	  for (CFint j = 0; j < _WI[_INVIBSPEI+v-1]; j++, icc++) {
	    CFuint ivt = _WI[_IVIBSPEI+icc-1]-1; 
	    _extraData.list.push_back(ivt); 
	    _extraData.cpVib[v] += 
	      _Y[ivt]*(_HVIBRUM[ivt] - _extraData.energyVib[ivt])/(_TVEC[v+1]*_EPS);  
	  }              
	}  
      }
      
      for (CFint is = 0; is < _NS; is++) {
	if ((_WI[_IATOMI+is-1]) == 0) { 
	  _extraData.cpElec[is] = 2.5*_Rgas/_MOLARMASSP[is];
	} else {
	  _extraData.cpElec[is] = (_HELECTUM[is] - _extraData.eElec[is])/(_TVEC[_nbTvib+1]*_EPS);
	}               
      }
      
    }
    // Chemical nonequilibrium
    else {
      CFdouble tempPert = temp*(1.0+_EPS);
      setDefaultTVarray(tempPert);
      FORTRAN_NAME(enthalpy)
	(_WR1,&_LWR1,_WI,&_LWI, &tempPert, &tempPert, &tempPert, _TVARRAY,
	 &pressure, _HTOTAL,_HTRANS,_HELECT,_HROT, _HVIBR,_HFORM);
      
      setDefaultTVarray(tempPert);
      FORTRAN_NAME(enthalpymass)
	(_WR1,&_LWR1,_WI,&_LWI, &tempPert, &tempPert, &tempPert, _TVARRAY,
	 &pressure, _HTOTALUM,_HTRANSUM,_HELECTUM,_HROTUM, _HVIBRUM,_HFORMUM); 
      
      
      if (presenceElectron()) {           // chemical nonequilibrium with free electrons
	const CFreal yRel = _Y[0]*_Rgas/_MOLARMASSP[0];
	_extraData.dHdT +=  2.5*yRel;
	_extraData.dEdT +=  1.5*yRel; 
      }
      
      if (!_noElectEnergy) { 
	for (CFint is = 0; is < _NS; ++is) {
	  const CFreal dEdT = 
	    _Y[is]*(_HVIBRUM[is] + _HELECTUM[is] - (_extraData.energyVib[is]
						    + _extraData.eElec[is]))/(_EPS*temp);
	  _extraData.dEdT += dEdT;
	  _extraData.dHdT += dEdT;
	}
      } else {
	for (CFint is = 0; is < _NS; ++is) {
	  const CFreal dEdT = _Y[is]*(_HVIBRUM[is] - _extraData.energyVib[is])/(_EPS*temp);   
	  _extraData.dEdT += dEdT;
	  _extraData.dHdT += dEdT;
	  
	}   
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

CFdouble MutationLibrary2OLD::density(CFdouble& temp,
				   CFdouble& pressure,
				   CFreal* tVec)
{
  CFdouble ND = 0.0;
  CFdouble rho = 0.0;
  CFdouble Te = getTe(temp,tVec);

  FORTRAN_NAME(numberd)(_WR1,&_LWR1,&pressure,&temp,&Te,_X,&ND);
  FORTRAN_NAME(density)(_WR1,&_LWR1,_X,&ND,&rho);
  return rho;
}

//////////////////////////////////////////////////////////////////////////////

CFdouble MutationLibrary2OLD::pressure(CFdouble& rho,
            CFdouble& temp,
            CFreal* tVec)
{
  CFdouble p = 0.0;
  if (_NE == 0) {
  // set the mixture pressure
    FORTRAN_NAME(pressure)(_WR1, &_LWR1, &rho, &temp, _MOLARMASSP, _Y, &p);
  }
  else {
    _electronPress = 0.0;   
    // if the electronic temperature is available, it is
    // the last entry in the array tVec (whose size is _nbTvib+1)  !! ask Andrea
    CFdouble Te = getTe(temp, tVec);
    FORTRAN_NAME(epressure)(_WR1, &_LWR1, &rho, &temp, &Te, _MOLARMASSP, _Y,
          &p, &_electronPress);
  }

  return p;
}

//////////////////////////////////////////////////////////////////////////////

CFdouble MutationLibrary2OLD::electronPressure(CFreal rhoE,
              CFreal tempE)
{
  return rhoE*tempE*_Rgas/_MOLARMASSP[0];
}

//////////////////////////////////////////////////////////////////////////////

CFdouble MutationLibrary2OLD::energy(CFdouble& temp,
          CFdouble& pressure)
{
  CFdouble ND = 0.0;
  CFdouble rho = 0.0;
  CFdouble MMass = 0.0;

  FORTRAN_NAME(numberd)(_WR1,&_LWR1,&pressure,&temp,&temp,_X,&ND);
  FORTRAN_NAME(density)(_WR1,&_LWR1,_X,&ND,&rho);
  FORTRAN_NAME(molarmass)(_WR1,&_LWR1,&rho,&ND,&MMass);

  setDefaultTVarray(temp);

  FORTRAN_NAME(enthalpy)(_WR1,&_LWR1,_WI,&_LWI,&temp,&temp,&temp,_TVARRAY,&pressure,
                         _HTOTAL,_HTRANS,_HELECT,_HROT,_HVIBR,_HFORM);

  const CFreal pOvRhoMM = pressure/rho*MMass;
  CFdouble intEnergy = 0.;
  for(CFint i = 0; i < _NS; ++i) {
    intEnergy += _X[i]*(_HTOTAL[i] - pOvRhoMM);
  }

  return intEnergy /= MMass;
}

//////////////////////////////////////////////////////////////////////////////

CFdouble MutationLibrary2OLD::enthalpy(CFdouble& temp,
                                   CFdouble& pressure)
{
  setDefaultTVarray(temp);

  FORTRAN_NAME(enthalpy)(_WR1,&_LWR1,_WI,&_LWI,&temp,&temp,&temp,_TVARRAY,&pressure,
       _HTOTAL,_HTRANS,_HELECT,_HROT,_HVIBR,_HFORM);

  CFdouble ND = 0.0;
  CFdouble rho = 0.0;
  CFdouble MMass = 0.0;

  // store the density
  FORTRAN_NAME(numberd)(_WR1,&_LWR1,&pressure,&temp,&temp,_X,&ND);
  FORTRAN_NAME(density)(_WR1,&_LWR1,_X,&ND,&rho);
  FORTRAN_NAME(molarmass)(_WR1,&_LWR1,&rho,&ND,&MMass);

  // sum up all the internal energies for each species
  CFdouble h = 0.0;
  for(CFint i = 0; i < _NS; ++i) {
    h += _X[i]*_HTOTAL[i];
  }
  return h /= MMass;
}

//////////////////////////////////////////////////////////////////////////////

void MutationLibrary2OLD::setElemFractions(const RealVector& yn)
{
  for (CFint ic = 0; ic < _NC; ++ic) {
    _Yn[ic] = yn[ic];

    if (!(_Yn[ic] >= 0.0 && _Yn[ic] <= 1.0)) {
      // cout << "Yn[ic] = " << Yn[ic] << endl;
      // abort();
    }

    assert(_Yn[ic] >= 0.0);
    assert(_Yn[ic] <= 1.0);
  }

  // Fills Xn according to Yn
  FORTRAN_NAME(nucmasstomolfrac)(_WR1, &_LWR1, _WI, &_LWI, _Yn, _Xn);
}

//////////////////////////////////////////////////////////////////////////////

void MutationLibrary2OLD::setElementXFromSpeciesY(const RealVector& ys)
{
  cf_assert(ys.size() == static_cast<CFuint>(_NS));
  
  CFreal massTot = 0.;
  for (CFint is = 0; is < _NS; ++is) {
    cf_assert (ys[is] < 1.1);
    
    const CFreal mm = ys[is]/_MOLARMASSP[is];
    massTot += mm;
    _X[is] = mm;
  }
  
  massTot = 1./massTot;
  for (CFint is = 0; is < _NS; ++is) {
    _X[is] *= massTot;
  }

  FORTRAN_NAME(compfract)(_X,&_NE,_Xn);
}

//////////////////////////////////////////////////////////////////////////////

void MutationLibrary2OLD::setElectronFraction(RealVector& ys)
{
  // get molar mass of species
  FORTRAN_NAME(setspeciesmolarmass)(_WR1, &_LWR1, _MOLARMASSP);

  // charge neutrality: xEl = sum(xIon)
  CFdouble yEl = 0.0;
  for (CFint is = 0; is < _NS; ++is) {
    if (_CHARGE[is] > 0) {
      yEl += ys[is] / _MOLARMASSP[is];
    }
  }

  yEl *= _MOLARMASSP[0]; // 1st species: electron
  ys[0] = yEl; // overwrite electron mass fraction
}

//////////////////////////////////////////////////////////////////////

void MutationLibrary2OLD::setSpeciesFractions(const RealVector& ys)
{
  // test this !!!!
  
  if (presenceElectron()) {
    setElectronFraction(const_cast<RealVector&>(ys));
  }
  
  for (CFint is = 0; is < _NS; ++is) {
    _Y[is] = ys[is];

    if (_Y[is] < 0.0) _Y[is] = 0.0;
    cf_assert(_Y[is] < 1.1);
  }

  // Fills X according to Y
  FORTRAN_NAME(specmasstomolfrac)(_WR1, &_LWR1, _WI, &_LWI, _Y, _X);
}

//////////////////////////////////////////////////////////////////////////////
      
void MutationLibrary2OLD::getSpeciesMolarFractions
(const RealVector& ys, RealVector& xs)
{
  CFreal massTot = 0.;
  for (CFint is = 0; is < _NS; ++is) {
    cf_always_assert(ys[is] < 1.1);
    
    const CFreal mm = ys[is]/_MOLARMASSP[is];
    massTot += mm;
    xs[is] = mm;
  }
  xs *= 1./massTot;
}
      
//////////////////////////////////////////////////////////////////////////////

void MutationLibrary2OLD::getSpeciesMassFractions
(const RealVector& xs, RealVector& ys)
{
  CFreal massTot = 0.;
  for (CFint is = 0; is < _NS; ++is) {
    cf_assert (xs[is] < 1.1);
    const CFreal mm = xs[is]*_MOLARMASSP[is];
    massTot += mm;
    ys[is] = mm;
  }
  ys /= massTot;
}
      
//////////////////////////////////////////////////////////////////////////////

void MutationLibrary2OLD::getSpeciesMassFractions(RealVector& ys)
{
  for (CFint is = 0; is < _NS; ++is) {
    ys[is] = _Y[is];
  }
}

//////////////////////////////////////////////////////////////////////////////

void MutationLibrary2OLD::getTransportCoefs(CFdouble& temp,
          CFdouble& pressure,
          CFdouble& lambda,
          CFdouble& lambdacor,
          RealVector& lambdael,
          RealMatrix& eldifcoef,
          RealVector& eltdifcoef)
{
  // NOTE: X must be set before the calling of this function!

  CFdouble ND = 0.0;

  setDefaultTVarray(temp);

  // Compute species enthalpies
  FORTRAN_NAME(enthalpy)(_WR1, &_LWR1, _WI, &_LWI, &temp, &temp, &temp,
       _TVARRAY, &pressure, _HTOTAL, _HTRANS, _HELECT,
       _HROT, _HVIBR, _HFORM);

  // Compute number density
  FORTRAN_NAME(numberd)(_WR1, &_LWR1, &pressure, &temp, &temp, _X, &ND);

  // Avoid error in case of vanishing molar fractions of species
  //FORTRAN_NAME(compotol2)(_X, &_Xlim, _XTOL);
  FORTRAN_NAME(compotol)(_X, &_TOL, _XTOL);

  // Fill WR2 vector, necessary if temperature is changed,
  // or composition in case of electron presence is changed ??? Here ???
  FORTRAN_NAME(collision)(_WR1, &_LWR1, _WR2, &_LWR2, &temp, &temp, &ND, _X);
  // Compute FIJ
  FORTRAN_NAME(correction)(_WR1, &_LWR1, _WR2, &_LWR2, _XTOL, &_sonine, _FIJ);
  // Compute transport properties
  FORTRAN_NAME(ltevefneu)(_WI, &_LWI, _WR1, &_LWR1, _WR2, &_LWR2, _HTOTAL,
        &temp, _X, &ND, _FIJ, &lambda, &lambdacor,
        _LAMBDAEL, _ELDIFCOEF, _ELTDIFCOEF);

  // Copying non scalar data
  for (CFint ic = 0; ic < _NC; ++ic) {
    lambdael[ic] = _LAMBDAEL[ic];
    eltdifcoef[ic] = _ELTDIFCOEF[ic];
    for (CFint jc = 0; jc < _NC; ++jc) {
      // Careful because C++ and Fortran label matrices differently!
      eldifcoef(ic,jc) = _ELDIFCOEF[jc*_NC+ic];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void MutationLibrary2OLD::getMassProductionTerm(CFdouble& temperature,
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
    CFdouble YLIM = 1.0e-20;
    // Flag to Fortran
    // CFint getJacobian = flagJac ? 1 : 0;
    CFdouble temp = max(_TminFix, temperature);

    // Fill species mass fractions
    for (CFint is = 0; is < _NS; ++is) {
      _Y[is] = ys[is];

      if (_Y[is] < 0.0) _Y[is] = 0.0;
      cf_assert(_Y[is] < 1.1);
      assert(_Y[is] <= 1.0);
    }

    // set the working _TVEC array
    setTVec(temp,tVec);
    
    CFint flagJ = (flagJac == true) ? 1 : 0;
    
    CFint iimolF77 = _iimol+1;
    
    FORTRAN_NAME(arrhenius) (_WR1, &_LWR1, _WR2, &_LWR2, _WR3, &_LWR3, _WR4, &_LWR4,
                             _WI, &_LWI, _Y, &YLIM, &pressure,
			     _TVEC, &rho, &_escape, _OMEGA, _OMEGACV, &_OMEGAI,
			     _OMEGAJACOB, &iimolF77, &flagJ);

    // Returning the reaction rates
    for (CFint is = 0; is < _NS; ++is) {
      omega[is] = _factorOmega*_OMEGA[is];
    }
    // Assembling Jacobian matrix of the source term corresponding to
    // primitive variables: rho_i, u, v,( w,) T, (Tv, Te)

    // suppose N = _NS-1 =>
    // jacobian (nbEq*nbEq) must contain the following derivatives in this order:

    // dOm(0)/dRho_s dOm(0)/dV dOm(0)/dT dOm(0)/dTv_m dOm(0)/dTe
    // ...
    // dOm(N)/dRho_s dOm(N)/dV dOm(N)/dT dOm(N)/dTv_m dOm(N)/dTe

    // -----_NS----- ---DIM--- ----1---- --_nbTvib--- --_nbTe---
    // with        dOm(s)/dV = 0


    //  _OMEGAJACOB is a 1-D array storing all the entries dOm(s)/dW (except dOm(s)/dV) in
    // in a row wise manner

    // this is meant for the case [rho_i v T Tv Te]
    if (flagJac) {
      const CFuint dim = PhysicalModelStack::getActive()->getDim();
      const CFuint nbEqs = jacobian.nbRows();
      assert(dim == nbEqs - _NS - _nbTvib - _nbTe - 1);
      const CFuint sizeJ = nbEqs - dim;
      const CFuint startT = _NS + dim;

      // this shouldn't be needed;
      jacobian = 0.0;

      const CFuint nns = (CFuint) _NS;
      for (CFuint is = 0; is < nns; ++is) {
        for (CFuint js = 0; js < nns; ++js) {
          jacobian(is,js) = _OMEGAJACOB[is*sizeJ+js];
        }
        for (CFuint js = startT; js < nbEqs; ++js) {
          const CFuint idx = js - dim;
          jacobian(is,js) = _OMEGAJACOB[is*sizeJ+idx];
        }
      }
      
      const CFuint nbTvTe = _nbTvib+_nbTe;
      for (CFuint v = 0; v < nbTvTe; ++v) {
        for (CFuint js = 0; js < nns; ++js) {
          jacobian(nns+dim+v+1,js) = _OMEGAJACOB[(nns+v)*sizeJ+js];
        }
        for (CFuint js = startT; js < nbEqs; ++js) {
          const CFuint idx = js - dim;
          jacobian(nns+dim+v+1,js) = _OMEGAJACOB[(nns+v)*sizeJ+idx];
        }
      }
    }
  }
  else {
    // Returning the reaction rates
    for (CFint is = 0; is < _NS; ++is) {
      omega[is] = 0.0;
    }
  } 
   
 }

//////////////////////////////////////////////////////////////////////////////

void MutationLibrary2OLD::getRhoUdiff(CFdouble& temperature,
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
  FORTRAN_NAME(numberd)(_WR1, &_LWR1, &pressure, &temp, &Te, _X, &ND);
  // Avoid error in case of vanishing molar fractions of species
  /// @TO DO AL: check this for ionized cases !!!!!!!!!!!!!!!!!!!!!!!!!!
  //  FORTRAN_NAME(compotol2)(_X, &_Xlim, _XTOL);
  //  if (presenceElectron()) {
  //     FORTRAN_NAME(compotol)(_X, &TOL, _XTOL);
  //  }
  //  else {
  FORTRAN_NAME(compotol2)(_X, &_Xlim, _XTOL);
  //  }

  // Fill WR2 vector, necessary if temperature is changed,
  // or composition in case of electron presence is changed ??? Here ???
  FORTRAN_NAME(collision)(_WR1, &_LWR1, _WR2, &_LWR2, &temp, &Te, &ND, _X);
  // Compute FIJ
  FORTRAN_NAME(correction)(_WR1, &_LWR1, _WR2, &_LWR2, _XTOL, &_sonine, _FIJ);

  // density of the mixture
  FORTRAN_NAME(density)(_WR1,&_LWR1,_X,&ND,&rho);
  // set the mixture molar mass
  FORTRAN_NAME(molarmass)(_WR1,&_LWR1,&rho,&ND,&MMass);

  // Set driving forces as gradients of molar fractions
  CFreal normMMassGradient = 0.0;
  for (CFint is = 0; is < _NS; ++is) {
    normMMassGradient += normConcGradients[is] / _MOLARMASSP[is];
  }
  normMMassGradient *= -MMass*MMass;

  for (CFint is = 0; is < _NS; ++is) {
    _DF[is] = (MMass*normConcGradients[is] + _Y[is]*normMMassGradient) /
      _MOLARMASSP[is];
  }

  // Compute the mass diffusion fluxes
  // FORTRAN_NAME(smneutsut)(_WR1, &_LWR1, _WR2, &_LWR2, _XTOL, &ND, _DF, _FIJ, _JDIF);

  CFdouble eamb = 0.;
  if (presenceElectron()) {
    FORTRAN_NAME(smd)(_WR1, &_LWR1, _WR2, &_LWR2, _WI, &_LWI, _XTOL, &temp,
          &Te, &ND, _DF, _FIJ, _JDIF, &eamb);
  }
  else {
    FORTRAN_NAME(smneutd)(_WR1, &_LWR1, _WR2, &_LWR2, _XTOL, &ND, _DF, _FIJ, _JDIF);
  }

  for (CFint is = 0; is < _NS; ++is) {
    rhoUdiff[is] = _JDIF[is];
  }

//   if (Te > 4000.) {
//     cout << "Te = " << Te << endl;
//     cout << "temp = " << temp << endl;
//     cout << "rhoU = " << rhoUdiff << endl;
//     abort();
//   }
}

//////////////////////////////////////////////////////////////////////////////

void MutationLibrary2OLD::getDij_fick(RealVector& dx,
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
  FORTRAN_NAME(numberd)(_WR1, &_LWR1, &pressure, &temperature, &temperature, _X, &ND);
  FORTRAN_NAME(density)(_WR1,&_LWR1,_X,&ND,&rho);

  // Fill WR2 vector, necessary if temperature is changed
  FORTRAN_NAME(collision)(_WR1, &_LWR1, _WR2, &_LWR2, &temperature, &temperature, &ND, _X);

  // First we fill the matrix of the binary diffusion coefficients
  CFuint ij;
  for (CFint is = 1; is < _NS+1; ++is) {
    for (CFint js = 1; js < _NS+1; ++js) {
      ij = ((is-1)*(2*_NS-is)+2*js)/2;
      Dij(is-1,js-1) = _WR2[_IBINIJ+ij-1] /ND;
    }
  }

  // set the mixture molar mass
  FORTRAN_NAME(molarmass)(_WR1,&_LWR1,&rho,&ND,&MMass);
  
  // Compute the Diffusion coefficient
  // Since we use Fick law with a constant coefficient for all spicies
  // It is necessary to choose a spicies that is present
  CFreal sum = 0.0;
  if (_NS == 5){
    if (_X[0] != 0.0){
      for (CFint js = 1; js < _NS; ++js)
	{
	  sum += (_X[js]/Dij(0,js));
	}
      sum *= _X[0];
      yi = _X[0]*_MOLARMASSP[0]/MMass;
    }
    
    else if (_X[1] != 0.0){
      for (CFint js = 0; js < _NS; ++js)
	{ if (js != 1)
	    sum += (_X[js]/Dij(1,js));
	}
      sum *= _X[1];
      yi = _X[1]*_MOLARMASSP[1]/MMass;
      
    }
    
    else if (_X[2] != 0.0){
      for (CFint js = 0; js < _NS; ++js)
	{	if (js != 2)
	    sum += (_X[js]/Dij(2,js));
	}
      sum *= _X[2];
      yi = _X[2]*_MOLARMASSP[2]/MMass;
      
    }
    
    else if(_X[3] != 0.0){
      for (CFint js = 0; js < _NS; ++js){
	if (js != 2)
	  sum += (_X[js]/Dij(3,js));
      }
      sum *= _X[3];
      yi = _X[3]*_MOLARMASSP[3]/MMass;
      
    }
    
    else if(_X[4] != 0.0){
      for (CFint js = 0; js < _NS; ++js){
	sum += (_X[js]/Dij(4,js));
      }
      sum *= _X[4];
      yi = _X[4]*_MOLARMASSP[4]/MMass;
      
    }
    
    
    else std::cout<<"There is a bug\n";
  }
  
  else if (_NS==11 ){ 
    if (_X[1] != 0.0){
      for (CFint js = 0; js < _NS; ++js)
	{ if (js != 1)
	    sum += (_X[js]/Dij(1,js));
	}
      sum *= _X[1];
      yi = _X[1]*_MOLARMASSP[1]/MMass;
      
    }
    
    else if (_X[2] != 0.0){
      for (CFint js = 0; js < _NS; ++js)
	{	if (js != 2)
	    sum += (_X[js]/Dij(2,js));
	}
      sum *= _X[2];
      yi = _X[2]*_MOLARMASSP[2]/MMass;
      
    }
    
    else if(_X[3] != 0.0){
      for (CFint js = 0; js < _NS; ++js){
	if (js != 3)
	  sum += (_X[js]/Dij(3,js));
      }
      sum *= _X[3];
      yi = _X[3]*_MOLARMASSP[3]/MMass;
      
    }
    
    else if(_X[4] != 0.0){
      for (CFint js = 0; js < _NS; ++js){
	if (js != 4)
	  sum += (_X[js]/Dij(4,js));
      }
      sum *= _X[4];
      yi = _X[4]*_MOLARMASSP[4]/MMass;
      
    }
    else if(_X[5] != 0.0){
      for (CFint js = 0; js < _NS; ++js){
	if (js != 5)
	  sum += (_X[js]/Dij(5,js));
      }
      sum *= _X[5];
      yi = _X[5]*_MOLARMASSP[5]/MMass;
      
    }
    // else if(_X[6] != 0.0){
    // 	for (CFint js = 0; js < _NS; ++js){
    // 	  if (js != 6)
    // 	  sum += (_X[js]/Dij(6,js));
    // 	}
    // 	sum *= _X[6];
    // 	yi = _X[6]*_MOLARMASSP[6]/MMass;
    
    // }
    // else if(_X[7] != 0.0){
    // 	for (CFint js = 0; js < _NS; ++js){
    // 	  if (js != 7)
    // 	  sum += (_X[js]/Dij(7,js));
    // 	}
    // 	sum *= _X[7];
    // 	yi = _X[7]*_MOLARMASSP[7]/MMass;
    
    // }
    //   else if(_X[8] != 0.0){
    // 	for (CFint js = 0; js < _NS; ++js){
    // 	  if (js != 8)
    // 	  sum += (_X[js]/Dij(8,js));
    // 	}
    // 	sum *= _X[8];
    // 	yi = _X[8]*_MOLARMASSP[8]/MMass;
    
    //     }
    //     else if(_X[9] != 0.0){
    // 	for (CFint js = 0; js < _NS; ++js){
    // 	  if (js != 9)
    // 	  sum += (_X[js]/Dij(9,js));
    // 	}
    // 	sum *= _X[9];
    // 	yi = _X[9]*_MOLARMASSP[9]/MMass;
    
    //     }
    // else if(_X[10] != 0.0){
    // 	for (CFint js = 0; js < _NS; ++js){
    // 	  if (js != 10)
    // 	  sum += (_X[js]/Dij(10,js));
    // 	}
    // 	sum *= _X[10];
    // 	yi = _X[10]*_MOLARMASSP[10]/MMass;
    
    //     }
    
    else std::cout<<"There is a bug\n";
    
  }
  
  // Diff_coeff = rho*(1.0 - yi)/sum;
  for (CFint is = 0; is < _NS; ++is){
    rhoUdiff[is] = -dx[is]*rho*(1.0 - yi)/sum;
  }
}

//////////////////////////////////////////////////////////////////////////////

void MutationLibrary2OLD::getGammaN(CFreal& m_GN)
{
  m_GN= _GammaN;
}

//////////////////////////////////////////////////////////////////////////////

void MutationLibrary2OLD::getGammaO(CFreal& m_GO)
{
  m_GO= _GammaO;
}

//////////////////////////////////////////////////////////////////////////////

void MutationLibrary2OLD::setSpeciesMolarFractions(const RealVector& xs){
std::cout<< "setSpeciesMolarFractions is not implemented in mutation2.0.0"<<endl;
}
//////////////////////////////////////////////////////////////////////////////

void MutationLibrary2OLD::getSpeciesTotEnthalpies(CFdouble& temp,
                 RealVector& tVec,
                 CFdouble& pressure,
                 RealVector& hsTot,
                 RealVector* hsVib,
                 RealVector* hsEl)
{

  if (tVec.size()== 0) {
    setDefaultTVarray(temp);
  }

  CFdouble Te = temp;
  if (tVec.size() > 0) {
    setTVarray(tVec);
    const CFuint TeID = tVec.size()-1;
    Te = (_nbTe == 0)? tVec[0] : tVec[TeID];
  }

  FORTRAN_NAME (enthalpy)(_WR1,&_LWR1,_WI,&_LWI,&temp,&Te,&temp,_TVARRAY,&pressure,
        _HTOTAL,_HTRANS,_HELECT,_HROT,_HVIBR,_HFORM);


 /* if (tVec.size() > 0) {
    *hsVib = 0.;
    CFuint icc=0;
    for (CFint i = 0; i < _nbTvib; ++i, ++icc) {
      const CFuint nvSize = _WI[_INVIBSPEI+i-1];
      for (CFuint j = 0; j < nvSize; ++j) {
  const CFuint idx = _WI[_IVIBSPEI+icc-1]-1;
  (*hsVib)[i] = _HVIBR[idx]/_MOLARMASSP[idx];
      }
    }
  }*/

  // WRONG !!!!!!!!!! NOT FLEXIBLE !!
  if (tVec.size() > 0) {
    for (CFuint i = 0; i < _molecIDs.size(); ++i) {
      (*hsVib)[i] = _HVIBR[_molecIDs[i]]/_MOLARMASSP[_molecIDs[i]];
    }
  }

  // returning the total enthalpies per unit mass of species
  for(CFint i = 0; i < _NS; ++i) {
    hsTot[i] = (!_noElectEnergy) ? _HTOTAL[i] / _MOLARMASSP[i] :
      (_HTOTAL[i] -_HELECT[i])/_MOLARMASSP[i];
  }
  
  if (presenceElectron()) {
    if(hsEl != CFNULL) {
      for(CFuint i = 0; i < _boltzmannIDs.size(); ++i) {
	const CFuint elEnergyID = _boltzmannIDs[i];
	(*hsEl)[elEnergyID] = _HELECT[elEnergyID]/_MOLARMASSP[elEnergyID];
      }
      
      (*hsEl)[0] = _HTRANS[0]/_MOLARMASSP[0];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void MutationLibrary2OLD::getSourceTermVT(CFdouble& temperature,
              RealVector& tVec,
              CFdouble& pressure,
              CFdouble& rho,
              RealVector& omegav,
              CFdouble& omegaRad)
{
  assert(omegav.size() >= 1); 

  CFdouble temp = max(_TminFix,temperature);
  CFdouble YTOL = 1e-20;

  setTVec(temp, tVec);

  const CFuint TeID = tVec.size()-1;
  CFreal Te = (_nbTe == 0)? tVec[0] : tVec[TeID];

  setTVarray(tVec);
  for (CFint i = 0; i < _NV; ++i) {
    _TT[i] = temp;
    _TE[i] = Te;
  }

  // N+/O+
  //      _OMEGAI = - MASS(93)*_OMEGA(93)*4.05D8*MASS(93)
  //     &          - MASS(95)*_OMEGA(95)*4.30D8*MASS(95)
  //      _OMEGAI = - MASS(48)*_OMEGA(48)*4.05D8*MASS(48)
  //    &          - MASS(50)*_OMEGA(50)*4.30D8*MASS(50)

  CFdouble ND = 0.0;
  CFdouble omegaet = 0.0;

  FORTRAN_NAME(numberd)(_WR1, &_LWR1, &pressure, &temp,&Te, _X,&ND);

  FORTRAN_NAME (collision)(_WR1, &_LWR1, _WR2, &_LWR2, &temp, &Te, &ND, _X);

  FORTRAN_NAME(enthalpy)(_WR1, &_LWR1, _WI, &_LWI, &temp, &Te, &temp, _TVARRAY, &pressure,
       _HTOTAL, _HTRANS, _HELECT, _HROT, _HVIBR, _HFORM);

  // this must be returned and summed to the global energy equation 
  omegaRad = 0.0;
  if (_escape > 0.) {
     for(CFint i = 0; i < _NS; ++i) {
       _HINT[i] = _HTOTAL[i] - _HTRANS[i];
     }
     FORTRAN_NAME (arrheniusrad)(_WR1, &_LWR1, _WR3, &_LWR3, _WI, &_LWI, _Y, &YTOL, &rho,
                              &Te, _HINT, &_escape, &omegaRad);
  }

  FORTRAN_NAME(enthalpy)(_WR1, &_LWR1, _WI, &_LWI, &temp, &temp, &temp, _TT, &pressure,
       _HT, _HTT, _HTE, _HTR, _HTV, _HTF);

  FORTRAN_NAME(enthalpy)(_WR1, &_LWR1, _WI, &_LWI, &Te, &Te, &Te, _TE, &pressure,
       _HE, _HET, _HEE, _HER, _HEV, _HEF);
  
  CFint iimolF77 = _iimol+1;
  
  FORTRAN_NAME(omegavtransfer) (&_LWR1, _WR1, &_LWR2, _WR2, &_LWR3, _WR3, &_LWR4, _WR4,
        &_LWI, _WI, &iimolF77, &pressure, _TVEC, _Y, &YTOL, &rho, _HTV, _HVIBR, _HEV,
        _STVIB, _OMEGAVE, _OMEGAVV, &omegaet);

  // Vibration-Dissociation coupling models
  switch(_cvModel) {
  case (0):
    cvCandler();
    break;
  case (1):
    cvGnoffo();
    break;
  default:
    cvCandler();
    break;
  }
  
  // omegaI overwritten
  if (_boltzmannIDs.size() == static_cast<CFuint>(_NS)) {
    //  // I have checked the dimensions MARCO
    _OMEGAI = -_OMEGA[8]*4.05e8*_MOLARMASSP[8] - _OMEGA[10]*4.30e8*_MOLARMASSP[10];
  }
/*
  I COMMENT THIS FOR CR TO USE THE SELF CONSISTENT MODEL !!!!
//     _OMEGAI = -_OMEGA[10]*4.05e8*_MOLARMASSP[10] - _OMEGA[12]*4.30e8*_MOLARMASSP[12];
      _OMEGAI = -_OMEGA[15]*4.05e8*_MOLARMASSP[15] - _OMEGA[17]*4.30e8*_MOLARMASSP[17];
*/
//  }

  if (_nbTe == 0) {
    for (CFint i = 0; i < _nbTvib; ++i) {
      omegav[i] = _STVIB[i] + _OMEGAVV[i] + _factorOmega*_OMEGACV[i];
    }

    // Park's model
    if (_includeElectronicEnergy && presenceElectron()) {
      omegav[max(0,_electrEnergyID)] += omegaet + _OMEGAI;
      for(CFuint i = 0; i < _boltzmannIDs.size(); ++i) {
        const CFuint elEnergyID = _boltzmannIDs[i];
	omegav[max(0,_electrEnergyID)] += _OMEGA[elEnergyID]*_HELECT[elEnergyID]/_MOLARMASSP[elEnergyID];
      }
    }
  }
  else {
    CFreal sumOmegaVE = 0.0;
    for (CFint i = 0; i < _nbTvib; ++i) {
      // this is wrong for multi-Tv: only N2 gets the OMEGAVEs to be generalized
//      _OMEGAVE[i] = 0.0;  // to remove
      omegav[i] = _STVIB[i] + _OMEGAVE[i] + _OMEGAVV[i] + _factorOmega*_OMEGACV[i];
      sumOmegaVE += _OMEGAVE[i];      
    }

//      _OMEGAI = 0.0;  // to remove
      omegav[omegav.size()-1] = omegaet - sumOmegaVE + _OMEGAI;

    if (_includeElectronicEnergy) {
      for(CFuint i = 0; i < _boltzmannIDs.size(); ++i) {
          const CFuint elEnergyID = _boltzmannIDs[i];
          omegav[omegav.size()-1] += _OMEGA[elEnergyID]*_HELECT[elEnergyID]/_MOLARMASSP[elEnergyID];
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void MutationLibrary2OLD::getMolarMasses(RealVector& mm)
{
  assert(mm.size() == static_cast<CFuint>(_NS));

  // check the units
  for (CFint i = 0; i < _NS; ++i) {
    mm[i] = _MOLARMASSP[i];
  }
}

//////////////////////////////////////////////////////////////////////////////

void MutationLibrary2OLD::cvCandler()
{
  CFuint icc=0;
  for (CFint i = 0; i < _nbTvib; ++i) {
    _OMEGACV[i] = 0.0;
    const CFuint nvSize = _WI[_INVIBSPEI+i-1];
    for (CFuint j = 0; j < nvSize; ++j, ++icc) {
      const CFuint idx = _WI[_IVIBSPEI+icc-1]-1;
      _OMEGACV[i] += _OMEGA[idx]*_HVIBR[idx]/_MOLARMASSP[idx];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void MutationLibrary2OLD::cvGnoffo()
{
  // in reality it should be  _prefDissFactor.size() == number of molecules!!!
  assert(_prefDissFactor.size() >= static_cast<CFuint>(_nbTvib));

  CFuint icc=0;
  for (CFint i = 0; i < _nbTvib; ++i) {
    _OMEGACV[i] = 0.0;
    const CFuint nvSize = _WI[_INVIBSPEI+i-1];
    for (CFuint j = 0; j < nvSize; ++j, ++icc) {
      const CFuint idx = _WI[_IVIBSPEI+icc-1]-1;

      /// @TODO this is correct only for air!! generalize it !!!
      /// CHECK THE UNITS I AM NOT SURE ABOUT WR3[_IDIS-1+i]
      _OMEGACV[i] += _prefDissFactor[i]*_OMEGA[idx]*_WR3[_IDIS-1+i]/_MOLARMASSP[idx];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void MutationLibrary2OLD::testAndWriteProperties()
{
  CFdouble lambdaTot = 0.0;
  CFdouble lambdaInth = 0.0;
  CFdouble lambdaTh = 0.0;
  CFdouble lambdaTe = 0.0;
  CFdouble ND = 0.0;
  CFdouble p = 1.;
  const CFreal DT = 1.;
//   CFuint counter = 0;

  while (p < 100000.) {
    const std::string fileName = "properties.dat" + Common::StringOps::to_str(static_cast<CFuint>(p));
    ofstream fout(fileName.c_str());

    fout << "TITLE = physical properties\n";
    fout << "VARIABLES = T p ";
    for (CFint i = 0; i < _NS; ++i) {
      fout << "H" << Common::StringOps::to_str(i) << " ";
    }
    fout << "gamma lambda eta" << endl;

    CFdouble T = 100.;

    while (T < 60000.) {
      setDefaultTVarray(T);

      FORTRAN_NAME(numberd)(_WR1,&_LWR1,&p,&T,&T,_X,&ND); //ND needed for collision below

      FORTRAN_NAME(enthalpy)(_WR1, &_LWR1, _WI, &_LWI, &T, &T, &T, _TVARRAY, &p,
           _HTOTAL, _HTRANS, _HELECT, _HROT, _HVIBR, _HFORM);

      CFdouble gamma = 0.0;
      FORTRAN_NAME(frozengammafast)(_WR1, &_LWR1, _WI, &_LWI, &T, &p, _X, &_EPS, &gamma);

      CFdouble tempEps1 = T*(1.0+_EPS);
      CFdouble tempEps = T*_EPS;
      setDefaultTVarray(tempEps1);

      FORTRAN_NAME(enthalpy)(_WR1, &_LWR1, _WI, &_LWI, &tempEps1, &tempEps1, &tempEps1, _TVARRAY,
           &p, _HTOTALP, _HTRANSP, _HELECTP, _HROTP, _HVIBRP, _HFORMP);

      for(int is=0 ; is < _NS ; is++) {
         _CPE[is]   = (_HELECTP[is] - _HELECT[is]) / (tempEps);
         _CPR[is]   = (_HROTP[is] - _HROT[is]) / (tempEps);
         _CPV[is]   = (_HVIBRP[is] - _HVIBR[is]) / (tempEps);
         _CPINT[is] = _CPE[is] + _CPR[is] + _CPV[is];
      }

      FORTRAN_NAME(compotol)(_X, &_TOL, _XTOL);
      //Update of kinetic data
      FORTRAN_NAME(collision)(_WR1, &_LWR1, _WR2, &_LWR2, &T, &T, &ND, _X);

      //A. Heavy particles properties

      //Eucken corrections for internal energy (no inelastic collisions)
      FORTRAN_NAME(lambdaint)(_WR1, &_LWR1, _WR2, &_LWR2, _CPINT, _XTOL, &lambdaInth);

      if (_lambdaAlgo == LAMBDACG) {
         FORTRAN_NAME(lambdachicg)(_WR1, &_LWR1, _WR2, &_LWR2, _XTOL, &lambdaTh, _CHIH);
      }
      else {
         FORTRAN_NAME(lambdachid)(_WR1, &_LWR1, _WR2, &_LWR2, _XTOL, &lambdaTh, _CHIH);
      }

      //B. Electron gas properties
      if(presenceElectron()){
        int i=3;
        FORTRAN_NAME(lambdae)(_WR1, &_LWR1, _WR2, &_LWR2, _X, &T, &lambdaTe, &i);
      }

      lambdaTot = lambdaInth + lambdaTh + lambdaTe;

      FORTRAN_NAME(numberd)(_WR1,&_LWR1,&p,&T,&T,_X,&ND); //ND needed for collision below

      FORTRAN_NAME(compotol)(_X, &_TOL, _XTOL);

      //Update of kinetic data
      FORTRAN_NAME(collision)(_WR1, &_LWR1, _WR2, &_LWR2, &T, &T, &ND, _X);

      //A. Heavy particles properties
      //Eucken corrections for internal energy (no inelastic collisions)
      CFreal eta = 0.0;
      FORTRAN_NAME(etacg)(_WR1, &_LWR1, _WR2, &_LWR2, _XTOL, &eta);

      fout << T << " " << p << " ";
      for (CFint ip = 0; ip < _NS; ++ip) {
  fout << _HTOTAL[ip] << " ";
  // fout << _HTOTAL[ip] << " " << _HELECT[ip] << " " << " " << _HVIBR[ip] << " ";
      }
      fout << gamma << " " <<  lambdaTot << " " << eta << endl;

      T += DT;
    }
    fout.close();
    p = p * 10.;
  }
}

//////////////////////////////////////////////////////////////////////////////

void MutationLibrary2OLD::testAndWriteSourceTerms()
{
  // Limit the mass composition of species !!!
  CFdouble YLIM = 1.0e-15;

  // Fill species mass fractions _Y
  // Fill species molar masses   _X

  RealVector tVec(2);
  CFdouble ND = 0.0;
  CFdouble p = 1.;
  const CFreal DT = 1.;
//   CFuint counter = 0;
  CFdouble rho = 0.;

  while (p < 100000.) {
    const std::string fileName = "sources.dat" + StringOps::to_str(static_cast<CFuint>(p));
    ofstream fout(fileName.c_str());

    fout << "TITLE = source terms\n";
    fout << "VARIABLES = T p ";
    for (CFint i = 0; i < _NS; ++i) {
      fout << "OMEGA" << StringOps::to_str(i) << " ";
    }
    for (CFint i = 0; i < _NS; ++i) {
      fout << "OMEGACV" << StringOps::to_str(i) << " ";
    }
    fout << "OMEGAI\n";

    CFdouble T = 100.;

    while (T < 60000.) {
      tVec = T;
      // set the working _TVEC array
      setTVec(T,tVec);

      FORTRAN_NAME(numberd)(_WR1,&_LWR1,&p,&T,&T,_X,&ND);
      FORTRAN_NAME(density)(_WR1,&_LWR1,_X,&ND,&rho);
 
      CFint flagJ = 0;      
      CFint iimolF77 = _iimol+1;

      FORTRAN_NAME(arrhenius) (_WR1, &_LWR1, _WR2, &_LWR2, _WR3, &_LWR3, _WR4, &_LWR4,
			       _WI, &_LWI, _Y, &YLIM, &p, _TVEC, &rho, &_escape, _OMEGA, 
                               _OMEGACV, &_OMEGAI, _OMEGAJACOB, &iimolF77, &flagJ);       

      fout << T << " " << p << " ";
      for (CFint ip = 0; ip < _NS; ++ip) {
	fout << _factorOmega*_OMEGA[ip] << " ";
      }
      for (CFint ip = 0; ip < _NS; ++ip) {
	fout << _factorOmega*_OMEGACV[ip] << " ";
      }
      fout << _OMEGAI << endl;

      T += DT;
    }
    fout.close();
    p = p * 10.;
  }
}

//////////////////////////////////////////////////////////////////////////////

void MutationLibrary2OLD::transportCoeffNEQ(CFreal& temperature, 
					    CFdouble& pressure,
					    CFreal* tVec, 
					    RealVector& normConcGradients,
					    CFreal& eta,
					    CFreal& lambdaTrRo, 
					    RealVector& lambdaInt,
					    RealVector& rhoUdiff)
{
  //composition must be called before!
  
  // --- at first collision integrals are computed --- //
  CFdouble temp = temperature;
  if (temp < 100.) {
    // cout << " MutationLibrary2OLD::etaCG() => Watch Out !! T <<<<< 100. " << temp << endl;
    temp = 100.;
  }
  
  CFdouble Te = getTe(temp,tVec);
  CFdouble ND = 0.0;
  
  FORTRAN_NAME(numberd)(_WR1,&_LWR1,&pressure,&temp,&Te,_X,&ND); //ND needed for collision below
  FORTRAN_NAME(compotol)(_X, &_TOL, _XTOL);
  FORTRAN_NAME(collision)(_WR1, &_LWR1, _WR2, &_LWR2, &temp, &Te, &ND, _X);
  
  // --- dynamic viscosity --- //
  //A. Heavy particles properties
  //Eucken corrections for internal energy (no inelastic collisions)
  FORTRAN_NAME(etacg)(_WR1, &_LWR1, _WR2, &_LWR2, _XTOL, &eta);
  
  // --- thermal conductivities --- //
  CFdouble lambdaRot = 0.0;
  CFdouble lambdaTh  = 0.0;
  CFdouble lambdaE   = 0.0;
  CFdouble lambdaTe  = 0.0;
  
  CFuint ic = 0;
  for (CFint i = 0; i < _NS; ++i) {
    const CFuint nvibMode = _WI[_IVIBI+i-1];
    for (CFuint j = 0; j < nvibMode; ++j, ++ic) {
      const CFuint ivib = _WI[_IVIBTEMPI+i-1]-2;
      //cf_assert(ivib < tVec->size());
      _TVARRAY[ic] = tVec[ivib];
      _TVIBEPS1[ic] = tVec[ivib] *(1.0+_EPS);
      _TVIBEPS[ic] = tVec[ivib]*_EPS;
    }
  }
  
  CFdouble TeEps = Te*_EPS;
  CFdouble TeEps1 = Te*(1.0+_EPS);
  CFdouble tempEps1 = temp*(1.0+_EPS);
  CFdouble tempEps = temp*_EPS;
  
  //Thermodynamic properties
  // - Species specific heat per unit mole
  FORTRAN_NAME(enthalpy)(_WR1, &_LWR1, _WI, &_LWI, &temp, &Te, &temp, _TVARRAY, &pressure,
			 _HTOTAL, _HTRANS, _HELECT, _HROT, _HVIBR, _HFORM);
  
  FORTRAN_NAME(enthalpy)(_WR1, &_LWR1, _WI, &_LWI, &tempEps1, &TeEps1, &tempEps1, _TVIBEPS1,
			 &pressure, _HTOTALP, _HTRANSP, _HELECTP, _HROTP, _HVIBRP, _HFORMP);
  
  for(CFint is = 0; is < _NS ; is++) {
    // _CPE[is]   = (_HELECTP[is] - _HELECT[is])/TeEps;
    _CPR[is]   = (_HROTP[is] - _HROT[is])/tempEps;
    
    // to generalize this !!!
    if (_WI[_IVIBI+is-1] == 1) {
      _CPVIB[is] = (_HVIBRP[is] - _HVIBR[is])/(tVec[_WI[_IVIBTEMPI+is-1]-2] *_EPS);
    }
    else {
      _CPVIB[is] = 0.0;
    }
  }
  
  // to test ...............
  const CFuint bsize = _boltzmannIDs.size();
  for(CFuint is = 0; is < bsize; is++) {
    const CFuint elEnergyID = _boltzmannIDs[is];
    _CPE[elEnergyID]   = (_HELECTP[elEnergyID] - _HELECT[elEnergyID])/TeEps;
  }
  
  //A. Heavy particles properties
  //Eucken corrections for rotational internal energy (no inelastic collisions)
  FORTRAN_NAME(lambdaint)(_WR1, &_LWR1, _WR2, &_LWR2, _CPR, _XTOL, &lambdaRot);
  
  FORTRAN_NAME(lambdaint)(_WR1, &_LWR1, _WR2, &_LWR2, _CPE, _XTOL, &lambdaE);
  
  FORTRAN_NAME(lambdavib)(_WR1, &_LWR1, _WR2, &_LWR2, _WI, &_LWI, _CPVIB, _XTOL,_LAMBDAVIB);
  
  if (_lambdaAlgo == LAMBDACG) {
    FORTRAN_NAME(lambdachicg)(_WR1, &_LWR1, _WR2, &_LWR2, _XTOL, &lambdaTh, _CHIH);
  }
  else {
    FORTRAN_NAME(lambdachid)(_WR1, &_LWR1, _WR2, &_LWR2, _XTOL, &lambdaTh, _CHIH);
  }
  
  //B. Electron gas properties
  int i=3;
  FORTRAN_NAME(lambdae)(_WR1, &_LWR1, _WR2, &_LWR2, _XTOL, &Te, &lambdaTe, &i);
  lambdaTrRo = lambdaRot + lambdaTh;
  
  for (CFint i = 0; i < _nbTvib; i++ ){
    lambdaInt[i] = _LAMBDAVIB[i];
  }
  
  if(_nbTe == 1) {
    lambdaInt[_nbTvib] = lambdaTe + lambdaE;
  }
  else {
    if (presenceElectron()) {
      lambdaInt[max(0,_electrEnergyID)] +=  lambdaE + lambdaTe;
    }
    else {
      if (_includeElectronicEnergy) {
	lambdaInt[max(0,_electrEnergyID)] +=  lambdaE;
      }
      lambdaTe = 0.0;
    }
  }
  
  
  // --- mass diffusion fluxes --- // 
  // NOTE: X must be set before the calling of this function!
  CFdouble rho = 0.0;
  CFdouble MMass = 0.0;
  
  // Compute FIJ
  FORTRAN_NAME(correction)(_WR1, &_LWR1, _WR2, &_LWR2, _XTOL, &_sonine, _FIJ);
  FORTRAN_NAME(density)(_WR1,&_LWR1,_X,&ND,&rho);  // density of the mixture
  FORTRAN_NAME(molarmass)(_WR1,&_LWR1,&rho,&ND,&MMass); // mixture molar mass
  
  // Set driving forces as gradients of molar fractions
  CFreal normMMassGradient = 0.0;
  for (CFint is = 0; is < _NS; ++is) {
    normMMassGradient += normConcGradients[is] / _MOLARMASSP[is];
  }
  normMMassGradient *= -MMass*MMass;
  
  for (CFint is = 0; is < _NS; ++is) {
    _DF[is] = (MMass*normConcGradients[is] + _Y[is]*normMMassGradient) /
      _MOLARMASSP[is];
  }
  
  CFdouble eamb = 0.;
  if (presenceElectron()) {
    FORTRAN_NAME(smd)(_WR1, &_LWR1, _WR2, &_LWR2, _WI, &_LWI, _XTOL, &temp,
		      &Te, &ND, _DF, _FIJ, _JDIF, &eamb);
    
    //FORTRAN_NAME(smsut)(_WR1, &_LWR1, _WR2, &_LWR2, _WI, &_LWI, _XTOL, &temp,
    //&Te, &ND, _DF, _FIJ, _JDIF, &eamb);
  }
  else {
    FORTRAN_NAME(smneutd)(_WR1, &_LWR1, _WR2, &_LWR2, _XTOL, &ND, _DF, _FIJ, _JDIF);
    
    // Compute the mass diffusion fluxes
    // FORTRAN_NAME(smneutsut)(_WR1, &_LWR1, _WR2, &_LWR2, _XTOL, &ND, _DF, _FIJ, _JDIF);
  }
  
  for (CFint is = 0; is < _NS; ++is) {
    rhoUdiff[is] = _JDIF[is];
  }
}

//////////////////////////////////////////////////////////////////////////////

void MutationLibrary2OLD::getSource(CFdouble& temperature,
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
  cout << "MutationLibrary2OLD::getSource() not implemented" << endl;		
  throw NotImplementedException(FromHere(), "MutationLibrary2OLD::getSource()");
}

//////////////////////////////////////////////////////////////////////////////

void MutationLibrary2OLD::getSourceEE(CFdouble& temperature,
				   RealVector& tVec,
				   CFdouble& pressure,
				   CFdouble& rho,
				   const RealVector& ys,
				   bool flagJac,
				   CFdouble& omegaEE)
  
{
  cout << "MutationLibrary2OLD::getSourceEE() not implemented" << endl;
  throw NotImplementedException(FromHere(), "MutationLibrary2OLD::getSourceEE()");
}

//////////////////////////////////////////////////////////////////////////////

} // namespace Mutation2OLD

} // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

