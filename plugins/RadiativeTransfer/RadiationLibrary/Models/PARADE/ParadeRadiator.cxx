#include <fstream>

#include "RadiativeTransfer/RadiationLibrary/Models/PARADE/ParadeRadiator.hh"
#include "RadiativeTransfer/RadiationLibrary/RadiationPhysicsHandler.hh"

#include "Common/CFLog.hh"
#include "Common/DebugFunctions.hh"
#include "Common/CFPrintContainer.hh"

#include "Environment/ObjectProvider.hh"
#include "Environment/CFEnv.hh"
#include "Environment/FileHandlerOutput.hh"
#include "Environment/FileHandlerInput.hh"
#include "Environment/SingleBehaviorFactory.hh"
#include "Environment/DirPaths.hh"
#include "Common/BadValueException.hh"
#include "Common/Stopwatch.hh"
#include "Common/PEFunctions.hh"
#include "Framework/MeshData.hh"
#include "Framework/PhysicalChemicalLibrary.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/PhysicalConsts.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Environment;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

namespace RadiativeTransfer {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<ParadeRadiator,
			    Radiator,
			    RadiativeTransferModule,
			    1>
paradeRadiatorProvider("ParadeRadiator");

//////////////////////////////////////////////////////////////////////////////

void ParadeRadiator::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< string >("LocalDirName","Name of the local temporary directories where Parade is run.");
  options.addConfigOption< string >                                 
    ("Namespace", "Namespace which must be used to run PARADE in parallel).");
  options.addConfigOption< bool, Config::DynamicOption<> >
    ("ReuseProperties", "Reuse existing radiative data (requires the same number of processors as in the previous run).");
  options.addConfigOption< CFreal >("Tmin","Minimum temperature.");
  options.addConfigOption< CFreal >("NDmin","Minimum number density.");
  options.addConfigOption< CFuint >("nbPoints","Number of Points to discretize the spectra per loop.");
  options.addConfigOption< string >("LibraryPath","Path to Parade data files.");
  options.addConfigOption< bool >                                 
    ("LTE", "True is it is a local thermodynamic equilibrium simulation.");
  options.addConfigOption< vector<bool> >                                 
    ("MolecularSpecies", "Array of flags indicating whether each species is molecular.");
  options.addConfigOption< string >("bandsDistr", "Banding distribution.");
  options.addConfigOption< bool >("Binning","Activation of binning.");
  options.addConfigOption< bool >("Banding","Activation of banding.");
  options.addConfigOption< bool >("WriteRadFileASCII", "Write the radiative coefficients to ASCII file (for debugging).");
  options.addConfigOption< bool >("SaveMemory", "Flag asking to parallelize as much as possible in order to save memory.");
  options.addConfigOption< bool > ("Equilibrium","Activation LTE");
 options.addConfigOption< bool > ("WriteHSNB","Writing HSNB input file");
  options.addConfigOption< bool > ("DataMassFractions","Activating mass fractions reading");
}
  
//////////////////////////////////////////////////////////////////////////////

ParadeRadiator::ParadeRadiator(const std::string& name) :
  Radiator(name),
  m_statesID(),
  m_rank(0),
  m_nbProc(1),
  m_inFileHandle(),
  m_outFileHandle(),
  m_paradeDir(),
  m_gridFile(),
  m_tempFile(),
  m_densFile(),
  m_radFile(),
  m_trTempID(),
  m_elTempID(),
  m_vibTempID(),
  m_isLTE(),
  m_massfraction()
{
  addConfigOptionsTo(this);
  
  m_localDirName = "Parade";
  setParameter("LocalDirName", &m_localDirName);
  
  m_namespace = "Default";
  setParameter("Namespace", &m_namespace);
  
  m_reuseProperties = false;
  setParameter("ReuseProperties", &m_reuseProperties);

  m_TminFix = 300.;
  setParameter("Tmin", &m_TminFix);

  m_ndminFix = 1e+10;
  setParameter("NDmin", &m_ndminFix);

  m_nbPoints= 10000;
  setParameter("nbPoints", &m_nbPoints);

  m_libPath = "";
  setParameter("LibraryPath", &m_libPath);

  m_isLTE = false;
  setParameter("LTE", &m_isLTE);

  m_massfraction = false;
  setParameter("DataMassFractions", &m_massfraction);
  
  m_molecularSpecies = vector<bool>();
  setParameter("MolecularSpecies", &m_molecularSpecies);

  m_bandsDistr = "equally";
  setParameter("bandsDistr", &m_bandsDistr);

  m_binning = false;
  setParameter("Binning",&m_binning);
  
  m_banding = false;
  setParameter("Banding",&m_banding);
    
  m_writeLocalRadCoeffASCII = false;
  setParameter("WriteRadFileASCII",&m_writeLocalRadCoeffASCII);

  m_writeHSNB  = true;
  setParameter("WriteHSNB",&m_writeHSNB);
  
  m_saveMemory = false;
  setParameter("SaveMemory",&m_saveMemory);

  m_Equilibrium = false;
  setParameter("Equilibrium",&m_Equilibrium);
}
  
//////////////////////////////////////////////////////////////////////////////

ParadeRadiator::~ParadeRadiator()
{
}

//////////////////////////////////////////////////////////////////////////////

void ParadeRadiator::configure ( Config::ConfigArgs& args )
{
  Radiator::configure(args);
  ConfigObject::configure(args);
}
      
//////////////////////////////////////////////////////////////////////////////
      
void ParadeRadiator::setup()
{
  CFLog(VERBOSE, "ParadeRadiator::setup() => START\n");
  
  Radiator::setup();
  
  m_rank   = PE::GetPE().GetRank(m_namespace);
  m_nbProc = PE::GetPE().GetProcessorCount(m_namespace);
  
  m_inFileHandle  = Environment::SingleBehaviorFactory<Environment::FileHandlerInput>::getInstance().create();
  m_outFileHandle = Environment::SingleBehaviorFactory<Environment::FileHandlerOutput>::getInstance().create();
    
  // create a m_localDirName-P# directory for the current process
  // cp parade executable and database inside the m_localDirName-P# directory 
  const string paradeDir = m_localDirName + "_" + m_radPhysicsPtr->getTRSname() + 
    "-P" + StringOps::to_str(PE::GetPE().GetRank(m_namespace));
  boost::filesystem::path paradePath(paradeDir);    
 
  m_paradeDir = Environment::DirPaths::getInstance().getWorkingDir() / paradePath; 
  CFLog(VERBOSE, "ParadeRadiator::setup() => m_paradeDir = " << m_paradeDir << "\n");
 
  boost::filesystem::path grid("grid.flo");
  m_gridFile = m_paradeDir / grid;
  
  boost::filesystem::path temp("temp.flo");
  m_tempFile = m_paradeDir / temp;
  
  boost::filesystem::path dens("dens.flo");
  m_densFile = m_paradeDir / dens;
  
  boost::filesystem::path rad("parade.rad");
  m_radFile = m_paradeDir / rad;

  m_library = PhysicalModelStack::getActive()->getImplementor()->
    getPhysicalPropertyLibrary<PhysicalChemicalLibrary>();

  if (m_library.isNotNull()) {
    const CFuint nbSpecies = m_library->getNbSpecies();
    cf_assert(nbSpecies > 0);
    CFLog(VERBOSE, "ParadeRadiator::setup() => nbSpecies = " << nbSpecies << "\n");
    m_mmasses.resize(nbSpecies);
    m_library->getMolarMasses(m_mmasses);
    
    /*
    m_mmasses[10] = 5.4858e-07;
    m_mmasses[0] = 0.028;
    m_mmasses[1] = 0.032;
    m_mmasses[2] = 0.03;
    m_mmasses[3] = 0.014;
    m_mmasses[4] = 0.016;
    m_mmasses[5] = 0.028;
    m_mmasses[6] = 0.032;
    m_mmasses[7] = 0.03;
    m_mmasses[8] = 0.014;
    m_mmasses[9] = 0.016;
    */
    
    for(int i = 0;i<nbSpecies;i++){
      CFLog(INFO,"The molar mass for species " << i << " is " << m_mmasses[i] << "\n");
    }
    
    m_avogadroOvMM.resize(nbSpecies);
    m_avogadroOvMM = PhysicalConsts::Avogadro()/m_mmasses;
   
    if (m_molecularSpecies.size() == 0) {
      vector<CFuint> moleculeIDs;
      m_library->setMoleculesIDs(moleculeIDs);
      cf_assert(moleculeIDs.size() > 0);
      
      m_molecularSpecies.resize(nbSpecies, false);
      for (CFuint i = 0; i < moleculeIDs.size(); ++i) {
	m_molecularSpecies[moleculeIDs[i]] = true;
      }
    }
    else {
      cf_assert(m_molecularSpecies.size() == nbSpecies);
    }
    
    const string msg = "ParadeRadiator::setup() => m_molecularSpecies: ";
    CFLog(VERBOSE, CFPrintContainer<vector<bool> >(msg, &m_molecularSpecies) << "\n");
  }
  
  // if this is a parallel simulation, only ONE process at a time sets the library
  runSerial<void, ParadeRadiator, &ParadeRadiator::setLibrarySequentially>(this, m_namespace);
  
  //get the statesID used for this Radiator
  m_radPhysicsPtr->getCellStateIDs( m_statesID );

  Framework::DataSocketSink < Framework::State* , Framework::GLOBAL > states =
  m_radPhysicsHandlerPtr->getDataSockets()->states;
  DataHandle<State*, GLOBAL> statesHandle = states.getDataHandle();
  
  if (statesHandle.size() != m_statesID.size()) {
    CFLog(ERROR, "ParadeRadiator::setup() => statesHandle.size() != m_statesID.size() => "
	  << statesHandle.size() << " != " <<  m_statesID.size() << "\n");
    cf_assert(statesHandle.size() == m_statesID.size());
  }
  
  // set up the TRS states
  m_pstates.reset
      (new Framework::DofDataHandleIterator<CFreal, State, GLOBAL>(statesHandle, &m_statesID));
  
  CFLog(VERBOSE, "ParadeRadiator::setup() => END\n");
}
  
//////////////////////////////////////////////////////////////////////////////

void ParadeRadiator::setLibrarySequentially()
{
  // the path to the data files must be specified (there is no default)
  if (m_libPath == "") {
    CFLog(ERROR, "ParadeLibrary::setLibrarySequentially() => libpath NOT SET\n");
    exit(1);
  }
  
  if (!m_reuseProperties) {
    // create a m_localDirName-P# directory for the current process
    // cp parade executable and database inside the m_localDirName-P# directory
    std::string command1   = "rm -fr " + m_paradeDir.string() + " ; mkdir " + m_paradeDir.string();
    Common::OSystem::getInstance().executeCommand(command1);
    std::string command2   = "cp " + m_libPath + "/Binary/parade " + m_paradeDir.string();
    Common::OSystem::getInstance().executeCommand(command2);
    std::string command3   = "cp -R " + m_libPath + "/Data " + m_paradeDir.string();
    Common::OSystem::getInstance().executeCommand(command3);
  }
  
  // all the following can fail if the format of the file parade.con changes
  // initialize some private data in this class (all processors, one by one, must run this)
}
  
//////////////////////////////////////////////////////////////////////////////
            
void ParadeRadiator::unsetup()
{
  Radiator::unsetup();
}
  
/////////////////////////////////////////////////////////////////////////////

void ParadeRadiator::setupSpectra(CFreal wavMin, CFreal wavMax)
{ 
  CFLog(INFO,"ParadeRadiator::computeProperties() => START\n");
  
  m_wavMin = wavMin;
  m_wavMax = wavMax;
  m_dWav = (m_wavMax-m_wavMin)/((static_cast<CFreal>(m_nbPoints)-1));
  
  Stopwatch<WallTime> stp;
  
  CFLog(VERBOSE, "ParadeRadiator::setupSpectra(wavmin: "
	<<m_wavMin<<", wavMax: "<<m_wavMax<<"), dWav: "<<m_dWav<<
        " nbPoints: "<<m_nbPoints<<"\n");
  
  if (!m_reuseProperties) {
    stp.start();
    // update the wavelength range inside parade.con
    // copy the modified files into the local Parade directories
    updateWavRange(wavMin, wavMax);
    
    // write the local grid, temperature and densities fields 
    writeLocalData();
    
    // run concurrently Parade in each local directory, one per process 
    runLibraryInParallel();
  }
  
  PE::GetPE().setBarrier(m_namespace);
  
  if (!m_reuseProperties) {
    CFLog(INFO,"ParadeRadiator::runLibraryInParallel() took " << stp.read() << "s\n");
  }
  
  // read in the radiative properties from parade.rad
  readLocalRadCoeff();
  
  if (m_binning && !m_banding) {
    if (m_nbBins > m_nbPoints) { CFLog(WARN, "ParadeRadiator:: WARNING - The number of Bins is higher than the number of spectrum points\n"); }
    computeBinning();
  }
  else if (!m_binning && m_banding) {
    if (m_nbBands > m_nbPoints) { CFLog(WARN, "ParadeRadiator:: WARNING - The number of Bands is higher than the number of spectrum points\n"); }
    computeBanding();
  }
  else if (m_binning && m_banding) {
    if (m_nbBins*m_nbBands > m_nbPoints) { CFLog(WARN, "ParadeRadiator:: WARNING - The total number of Bins is higher than the number of spectrum points\n"); }
    computeBinningBanding();
  }
  
  PE::GetPE().setBarrier(m_namespace);
  
  CFLog(INFO,"ParadeRadiator::computeProperties() => END\n");
}
  
//////////////////////////////////////////////////////////////////////////////
      
void ParadeRadiator::updateWavRange(CFreal wavMin, CFreal wavMax)
{
  CFLog(VERBOSE,"ParadeRadiator::updateWavRange() => START\n");

  boost::filesystem::path confFile = Environment::DirPaths::getInstance().getWorkingDir() / "parade.con";

  if (PE::GetPE().GetRank(m_namespace) == 0) {
    // all the following can fail if the format of the file parade.con changes

    boost::filesystem::path confBkp  = Environment::DirPaths::getInstance().getWorkingDir() / "parade.con.bkp";

    CFLog(INFO, "ParadeRadiator: computing wavelength range [" << wavMin << ", " << wavMax << "]\n");
    
    // back up the last parade.con file
    std::string command = "cp " + confFile.string() + " " + confBkp.string();
    Common::OSystem::getInstance().executeCommand(command); 
   
    SelfRegistPtr<FileHandlerInput> fhandleIn   = SingleBehaviorFactory<FileHandlerInput>::getInstance().create();
    SelfRegistPtr<FileHandlerOutput> fhandleOut = SingleBehaviorFactory<FileHandlerOutput>::getInstance().create();

    ifstream& fin  = fhandleIn->open(confBkp);  
    ofstream& fout = fhandleOut->open(confFile);     
  
    // read from the original parade.con file 
   // ifstream fin("parade.con.bkp");
    //ofstream fout("parade.con");
    
    string line;
    while (getline(fin,line)) {
      // detect and store the position in the input file where cell data will be inserted
      bool isModified = false;
      
      // update the min wavelength
      size_t pos1 = line.find("wavlo");
      if (pos1 != string::npos) {
      fout << wavMin << " wavlo [A]" << endl;
	isModified = true;
      }
      
      // update the max wavelength
      size_t pos2 = line.find("wavhi");
      if (pos2 != string::npos) {
        fout << wavMax << " wavhi [A] " << endl;
	isModified = true;
      }
      
      // update the number of spectral points (= number of wavelengths considered)
      size_t pos3 = line.find("npoints");
      if (pos3 != string::npos) {
        fout << m_nbPoints << " npoints" << endl;
	isModified = true;
      }
      
      if (!isModified) {
	fout << line << endl;
      }
    }
    fin.close();
    fout.close();
  }
  
  PE::GetPE().setBarrier(m_namespace);
  
  // each processor copies the newly updated parade.con to its own directory
  std::string command   = "cp " + confFile.string() + " " + m_paradeDir.string();
  Common::OSystem::getInstance().executeCommand(command);

  CFLog(VERBOSE,"ParadeRadiator::updateWavRange() => END\n");
}
      
//////////////////////////////////////////////////////////////////////////////
      
void ParadeRadiator::writeLocalData()
{ 
  CFLog(VERBOSE, "ParadeRadiator::writeLocalData() => START\n");
  
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  const CFuint nbTemps = m_radPhysicsHandlerPtr->getNbTemps();
  const CFuint tempID = m_radPhysicsHandlerPtr->getTempID();
  
  const CFuint totalNbPoints = m_pstates->getSize();
  CFuint nbPoints = totalNbPoints;
  CFuint countNode = 0;
  // if the full mesh is in this process than you only write a part of its data
  if (fullGridInProcess()) {
    const CFuint nbPointsPerProc = nbPoints/m_nbProc;
    nbPoints = (m_rank < m_nbProc-1) ? nbPointsPerProc : nbPointsPerProc + nbPoints%m_nbProc;
    for (CFuint i = 0; i < m_rank; ++i) {
      countNode += nbPointsPerProc;
    }
  }
  
  const CFuint startNode = countNode;
  
  CFLog(INFO, "ParadeRadiator::writeLocalData() => nbPoints = " 
	<< nbPoints << ", nbTemps = " << nbTemps 
	<< ", tempID = " << tempID << "\n");
  
  // write the mesh file
  ofstream& foutG = m_outFileHandle->open(m_gridFile);
  foutG << "TINA" << endl;
  foutG << 1 << " " << nbPoints << endl;
  for (CFuint i =0; i < nbPoints; ++i, ++countNode) {
    foutG.precision(14);
    foutG.setf(ios::scientific,ios::floatfield);
    cf_assert(countNode < totalNbPoints);
    CFreal *const node = m_pstates->getNode(countNode);
    
    if (dim == DIM_1D) {
      foutG << node[XX] << " " << 0.0 << " " << 0.0  << endl;
    }
    if (dim == DIM_2D) {
      foutG << node[XX] << " " << node[YY] << " " << 0.0  << endl;
    }
    if (dim == DIM_3D) {
      foutG << node[XX] << " " << node[YY] << " " << node[ZZ] << endl;
    }
  } 
  foutG.close();
  
  CFLog(INFO, "ParadeRadiator::writeLocalData() => written coordinates for cells [" 
	<< startNode << ", " << countNode << "]\n");
  
  // write the temperatures
  ofstream& foutT = m_outFileHandle->open(m_tempFile);
  //const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  foutT << 1 << " " << nbPoints << " " <<  nbTemps << endl;
  // here it is assumed that the temperatures are the LAST variables 
  countNode = startNode;
  for (CFuint i = 0; i < nbPoints; ++i, ++countNode) {
    cf_assert(countNode < totalNbPoints);
    CFreal *const currState = m_pstates->getState(countNode);
    for (CFuint t = 0; t < nbTemps; ++t) {
      foutT.precision(14);
      foutT.setf(ios::scientific,ios::floatfield);
      const CFreal temp = std::max(currState[tempID + t], m_TminFix);
      foutT << temp << " ";
    }
    foutT << endl;
  }
  foutT.close();
  
  CFLog(INFO, "ParadeRadiator::writeLocalData() => written temperature for cells [" 
	<< startNode << ", " << countNode << "]\n");
  
  // write the number densities
  ofstream& foutD = m_outFileHandle->open(m_densFile);
  const CFuint nbSpecies = m_library->getNbSpecies();
  foutD << 1 << " " << nbPoints << " " << nbSpecies << endl;  

  //write species partial density
 foutD.precision(14);
 foutD.setf(ios::scientific,ios::floatfield);
 countNode = startNode;

 if(!m_massfraction){
 
  if(m_isLTE) {	  
    // LTE: composition is computed though the chemical library
    RealVector x(nbSpecies);
    RealVector y(nbSpecies);	 
    CFreal rhoi = 0.;
    for (CFuint i = 0; i < nbPoints; ++i, ++countNode) {
      cf_assert(countNode < totalNbPoints);
      CFreal *const currState = m_pstates->getState(countNode);
      CFreal temp  = currState[tempID];
      CFreal press = currState[0];
      m_library->setComposition(temp, press, &x);
      m_library->getSpeciesMassFractions(x,y);
      const CFreal rho = m_library->density(temp, press, CFNULL);

      for (CFuint t = 0; t < nbSpecies; ++t) {
        // number Density = partial density/ molar mass * Avogadro number
        rhoi = rho*y[t];
        const CFreal nb = std::max(rhoi*m_avogadroOvMM[t],m_ndminFix);
        foutD << nb << " ";
      }
      foutD << endl;
    }
  } 
  else{ 
    // NEQ: species partial density is a system state
    // here it is assumed that the species densities are the FIRST variables 
    for (CFuint i =0; i < nbPoints; ++i, ++countNode) {
      cf_assert(countNode < totalNbPoints);
      CFreal *const currState = m_pstates->getState(countNode);
      for(CFuint t=0; t<nbSpecies; ++t){
	const CFreal nb = std::max(currState[t]*m_avogadroOvMM[t],m_ndminFix);
	foutD << nb << " ";
      }
      foutD << endl;
    }
  }
 }
 else{
  if(m_isLTE) {	  
    // LTE: composition is computed though the chemical library
    RealVector x(nbSpecies);
    RealVector y(nbSpecies);	 
    CFreal rhoi = 0.;
    for (CFuint i = 0; i < nbPoints; ++i, ++countNode) {
      cf_assert(countNode < totalNbPoints);
      CFreal *const currState = m_pstates->getState(countNode);
      CFreal temp  = currState[tempID];
      CFreal press = currState[0];
      m_library->setComposition(temp, press, &x);
      m_library->getSpeciesMassFractions(x,y);
      const CFreal rho = m_library->density(temp, press, CFNULL);

      for (CFuint t = 0; t < nbSpecies; ++t) {
        // number Density = partial density/ molar mass * Avogadro number
        rhoi = rho*y[t];
        const CFreal nb = std::max(rhoi*m_avogadroOvMM[t],m_ndminFix);
        foutD << nb << " ";
      }
      foutD << endl;
    }
  } 
  else{ 
    // NEQ: species partial density is a system state
    // here it is assumed that the species densities are the FIRST variables 
    for (CFuint i =0; i < nbPoints; ++i, ++countNode) {
      cf_assert(countNode < totalNbPoints);
      CFreal *const currState = m_pstates->getState(countNode);
      CFreal temp  = currState[tempID];
      CFuint pressID = currState[0];
      CFreal press = currState[pressID];
CFLog(INFO,"The pressure in this state is = " << press << "\n");
      const CFreal rho = m_library->density(temp, press, CFNULL);
      for(CFuint t=0; t<nbSpecies; ++t){
	const CFreal nb = std::max(currState[t]*rho*m_avogadroOvMM[t],m_ndminFix);
	foutD << nb << " ";
      }
      foutD << endl;
    }
  }
 }

  

if(m_writeHSNB){

ofstream fout("HSNB_input.out");

 fout << "Number of points" << endl;
 fout << nbPoints << endl;

 fout << "Geometry type" << endl;
 fout << "CARTESIAN" << endl;

 // boost::filesystem::path confFile = Environment::DirPaths::getInstance().getWorkingDir() / "parade.con";
 // SelfRegistPtr<FileHandlerInput> fhandleIn   = SingleBehaviorFactory<FileHandlerInput>::getInstance().create();
 // ifstream& fin  = fhandleIn->open(confFile);

 fout << " #Loc Tr Tv P N2 O2 NO N O N2+ O2+ NO+ N+ O+ e-" << endl;
CFLog(INFO,"HSNB_radiator::writeLocalData()=> Writing file" << "\n");

//fout.precission(14);
//fout.setf(ios::scientific,ios::floatfield);

const CFuint nbSpecies = m_library->getNbSpecies();


for (CFuint i=0; i < nbPoints; ++i) {

if(m_isLTE) {
//In LTE, the composition is computed through the chemical library
RealVector x(nbSpecies);
RealVector y(nbSpecies);
CFreal rhoi = 0.;

      CFreal *const currState = m_pstates->getState(i);
      CFreal temp  = currState[tempID];
      CFreal press = currState[0];
      m_library->setComposition(temp, press, &x);
      m_library->getSpeciesMassFractions(x,y);
      const CFreal rho = m_library->density(temp, press, CFNULL);

      fout << i << " ";
      fout << temp << " ";
      fout << temp << " ";
      fout << press << " ";

// Temperature is written two times because we make equal the rotational and vibrational temperatures for the moment

      for (CFuint t = 0; t < nbSpecies; ++t) {
        // number Density = partial density/ molar mass * Avogadro number
        rhoi = rho*y[t];
        //const CFreal nb = std::max(rhoi*m_avogadroOvMM[t],m_ndminFix);
        fout << rhoi << " ";
      }
 }


else{

CFreal *const currState = m_pstates->getState(i);
CFreal temp  = currState[tempID];
CFreal press = currState[0];

     fout << i << " ";
     fout << temp << " ";
     fout << temp << " ";
     fout << press << " ";

  for(CFuint t=0; t<nbSpecies; ++t){
     const CFreal nb = std::max(currState[t]*m_avogadroOvMM[t],m_ndminFix);
     fout << currState[t] << " ";
  }

 }

fout << endl;
 }
 fout.close();
 }

  CFLog(INFO, "ParadeRadiator::writeLocalData() => written densities for cells [" 
	<< startNode << ", " << countNode << "]\n");
  
  CFLog(VERBOSE, "ParadeRadiator::writeLocalData() => START\n");
}   
      
//////////////////////////////////////////////////////////////////////////////
  
void ParadeRadiator::readLocalRadCoeff()
{
  CFLog(VERBOSE, "ParadeRadiator::readLocalRadCoeff() => START\n");
  
  /// array storing absorption and emission coefficients
  /// wavelength       = m_radCoeff(local state ID, spectral point idx*3)
  /// emission coeff   = m_radCoeff(local state ID, spectral point idx*3+1)
  /// absorption coeff = m_radCoeff(local state ID, spectral point idx*3+2)

  fstream& fin = m_inFileHandle->openBinary(m_radFile);
  //string nam = "file-" + StringOps::to_str(PE::GetPE().GetRank());
  
  int one = 0;
  fin.read((char*)&one, sizeof(int));
  //fout << "one = "<< one;
  cf_assert(one == 1);
  
  int nbCells = 0;
  fin.read((char*)&nbCells, sizeof(int));
  
  CFLog(VERBOSE,"ParadeRadiator::readLocalRadCoeff() => nbCells = " << nbCells << "\n");
  if (!fullGridInProcess()) {
    cf_assert(nbCells == (int)m_pstates->getSize()  );
  }
  else {
    cf_assert(nbCells <= (int)m_pstates->getSize()  );
  }
  
  vector<int> wavptsmx(3);
  fin.read((char*)&wavptsmx[0], 3*sizeof(int));
  //fout << "wav3 =  "<< wavptsmx[0] << " " << wavptsmx[1] << " " << wavptsmx[2]<< endl;
  cf_assert(wavptsmx[0] == (int) m_nbPoints);
  const CFuint totalNbCells = m_pstates->getSize();
  const CFuint sizeLocalCells = (!m_saveMemory) ? totalNbCells : nbCells;
  m_data.resize(sizeLocalCells*m_nbPoints*3);
  
  // here if the process stores the full mesh (as in the FV DOM algorithm), each processor 
  // reads only a portion of the data
  LocalArray<CFreal>::TYPE partialData;
  LocalArray<CFreal>::TYPE* currData = &m_data;
  if (fullGridInProcess() && !m_saveMemory) {
    partialData.resize(nbCells*m_nbPoints*3);
    currData = &partialData;
  }
  
  double etot = 0.;
  int wavpts = 0;
  const CFuint sizeCoeff = m_nbPoints*3;
  for (int iCell = 0 ; iCell < nbCells; ++iCell) {
    fin.read((char*)&etot, sizeof(double));
    //fout << "etot = " << etot << endl;
    fin.read((char*)&wavpts, sizeof(int));
    cf_assert(wavpts == (int)m_nbPoints);
    //fout << "wavpt = " <<  wavpts << endl;
    // this reads [wavelength, emission, absorption] for each cell
    fin.read((char*)&((*currData)[iCell*sizeCoeff]), sizeCoeff*sizeof(double));
  }
  fin.close();
  
  if (fullGridInProcess() && !m_saveMemory) {
    // in case the full mesh is stored in each process, since we have read in only a part 
    // of spectral data, we need now to gather all data so that each process has the full 
    // spectra before performing binning 
    // AL: this is just an intermediate step towards a truly parallel binning calculation
    CFuint minSizeToSend = 0;
    CFuint maxSizeToSend = 0;
    vector<int> recvCounts(m_nbProc, 0);
    vector<int> displs(m_nbProc, 0);
    computeRecvCountsDispls(totalNbCells, sizeCoeff, minSizeToSend, maxSizeToSend, recvCounts, displs);
    cf_assert(currData->size() <= maxSizeToSend);
    cf_assert(currData->size() >= minSizeToSend);
    
    MPIError::getInstance().check
      ("MPI_Allgatherv", "ParadeRadiator::readLocalRadCoeff()",
       MPI_Allgatherv(&(*currData)[0], currData->size(), MPIStructDef::getMPIType(&(*currData)[0]),
		      &m_data[0], &recvCounts[0], &displs[0],  MPIStructDef::getMPIType(&m_data[0]),
		      PE::GetPE().GetCommunicator(m_namespace)));
    
    CFLog(VERBOSE, "ParadeRadiator::readLocalRadCoeff() => after MPI_Allgatherv()\n");
  }
  
  CFLog(INFO, "ParadeRadiator::readLocalRadCoeff() => read data for all cells\n");
  
  if (m_writeLocalRadCoeffASCII) {
    writeLocalRadCoeffASCII(nbCells);
  }
  
  CFLog(VERBOSE, "ParadeRadiator::readLocalRadCoeff() => END\n");
}

//////////////////////////////////////////////////////////////////////////////
  
void ParadeRadiator::writeLocalRadCoeffASCII(const CFuint nbCells)
{
  CFLog(VERBOSE,"ParadeRadiator::writeLocalRadCoeffASCII() => START\n");
  
  boost::filesystem::path rad("parade.rad.ASCII");
  boost::filesystem::path outfile = m_paradeDir / rad;
  ofstream& fout = m_outFileHandle->open(outfile);
  
  // Writing m_data into a ASCII file.
  const CFuint nbCols = m_nbPoints*3;
  RealMatrix dataMat(m_data.size()/nbCols, nbCols, &m_data[0]);
  for(CFuint p=0;p<m_nbPoints;++p) {
    const CFuint p3 = p*3;
    for(CFuint m=0;m<nbCells;++m) {
      fout << dataMat(m,p3) << " ";
      fout << dataMat(m,p3+2) << " ";
      fout << dataMat(m,p3+1) << " ";
    }
  }
  fout.close();
  
  // if (PE::GetPE().GetRank("Default") == 0) {
  //   ofstream fout1("inwav.plt"); 
  //   fout1 << "TITLE = Original spectrum of radiative properties\n";
  //   fout1 << "VARIABLES = Wavl EmCoef AbCoef \n";
  //   for (CFuint i = 0; i < data.nbCols()/3; ++i) {
  //     fout1 << m_data(0,i*3) << " " << m_data(0,i*3+1) << " " << m_data(0,i*3+2) << endl;
  //   }
  //   fout1.close();
  // }
  
  //   // Convert spectrum to frequency space
  //   const CFreal c = PhysicalConsts::LightSpeed();
  //   const CFuint stride = data.nbCols()/3;
  //   for (CFuint iCell = 0 ; iCell < nbCells; ++iCell) {
  //     for (CFuint iw = 0; iw < stride; ++iw) {
  //       const CFreal lambda = data(iCell, iw*3)*1.e-10; // conversion from [A] to [m]
  //       // convert emission from wavelength spectrum to frequency spectrum
  //       data(iCell, iw*3+1) *= c/(lambda*lambda);
  //       // override wavelength with frequency (nu = c/lambda)
  //       data(iCell, iw*3) = c/lambda; 
  //     }
  //   }
  
  //   if (PE::GetPE().GetRank("Default") == 0) {
  //     ofstream fout1("infreq.txt"); 
  //     fout1 << "TITLE = Original spectrum of radiative properties\n";
  //     fout1 << "VARIABLES = Freq EmCoef AbCoef \n";
  //     for (CFuint i = 0; i < data.nbCols()/3; ++i) {
  //       fout1 << data(0,i*3) << " " << data(0,i*3+1) << " " << data(0,i*3+2) << endl;
  //     }
  //     fout1.close();
//   }
  
  CFLog(VERBOSE, "ParadeRadiator::writeLocalRadCoeffASCII() => END\n");
}

//////////////////////////////////////////////////////////////////////////////
  
void ParadeRadiator::getSpectralIdxs(CFreal lambda, CFuint& idx1, CFuint& idx2)
{
  //assumes constant wavelength discretization
  cf_assert(lambda <= m_wavMax);
  const CFreal idx = (m_nbPoints-1) * ( lambda - m_wavMin )/( m_wavMax - m_wavMin );
  
  idx1 = (CFuint)std::floor(idx);
  idx1 = std::min(m_nbPoints -1, idx1);
  idx2 = (CFuint)std::ceil(idx);
  idx2 = std::min(m_nbPoints -1, idx2);
  
  cf_assert(idx1 < m_nbPoints);
  cf_assert(idx2 < m_nbPoints);
}
  
//////////////////////////////////////////////////////////////////////////////
  
CFreal ParadeRadiator::getEmission(CFreal lambda, RealVector &s_o)
{
  CFuint spectralIdx1 = 0;
  CFuint spectralIdx2 = 0;
  const CFuint stateIdx = m_radPhysicsHandlerPtr->getCurrentCellTrsIdx();
  getSpectralIdxs(lambda, spectralIdx1, spectralIdx2);
  
  const CFuint nbCols = m_nbPoints*3;
  const CFuint nbRows = m_data.size()/nbCols;
 
  cf_assert(stateIdx < nbRows);
  cf_assert(spectralIdx1*3   < nbCols);
  cf_assert(spectralIdx1*3+1 < nbCols);
  cf_assert(spectralIdx2*3   < nbCols);
  cf_assert(spectralIdx2*3+1 < nbCols);

  const CFreal x0 = m_data[stateIdx*nbCols + spectralIdx1*3];
  const CFreal y0 = m_data[stateIdx*nbCols + spectralIdx1*3+1];
  const CFreal x1 = m_data[stateIdx*nbCols + spectralIdx2*3];
  const CFreal y1 = m_data[stateIdx*nbCols + spectralIdx2*3+1];

  //linear interpolation
  return y0 + (y1-y0) * (lambda - x0) / (x1-x0);
}

//////////////////////////////////////////////////////////////////////////////

CFreal ParadeRadiator::getAbsorption(CFreal lambda, RealVector &s_o)
{
  CFuint spectralIdx1 = 0; 
  CFuint spectralIdx2 = 0;
  const CFuint stateIdx = m_radPhysicsHandlerPtr->getCurrentCellTrsIdx();
  getSpectralIdxs(lambda, spectralIdx1, spectralIdx2);
  
  const CFuint nbCols = m_nbPoints*3;
  const CFuint nbRows = m_data.size()/nbCols;
  cf_assert(stateIdx < nbRows);
  cf_assert(spectralIdx1*3   < nbCols); 
  cf_assert(spectralIdx1*3+2 < nbCols);
  cf_assert(spectralIdx2*3   < nbCols);
  cf_assert(spectralIdx2*3+2 < nbCols);
  
  const CFreal x0 = m_data[stateIdx*nbCols + spectralIdx1*3];
  const CFreal y0 = m_data[stateIdx*nbCols + spectralIdx1*3+2];
  const CFreal x1 = m_data[stateIdx*nbCols + spectralIdx2*3];
  const CFreal y1 = m_data[stateIdx*nbCols + spectralIdx2*3+2];
  
  //linear interpolation
  return y0 + (y1-y0) * (lambda - x0) / (x1-x0);
}

//////////////////////////////////////////////////////////////////////////////

void ParadeRadiator::computeEmissionCPD()
{
  CFLog(VERBOSE, "ParadeRadiator::computeEmissionCPD() => start\n");
  const CFuint totalNbCells = m_pstates->getSize();
  const CFuint nbCpdPoints = m_nbPoints;
  m_spectralLoopPowers.resize( totalNbCells );
  m_cpdEms.resize( totalNbCells * nbCpdPoints );

  const CFuint nbCols = m_nbPoints*3;
  RealMatrix dataMat(m_data.size()/nbCols, nbCols, &m_data[0]);
  
  CFuint nbCells = totalNbCells;
  CFuint offsetStateID = 0;
  CFreal* spectralLoopPowersCurr = &m_spectralLoopPowers[0];
  CFreal* cpdEmsCurr = &m_cpdEms[0];
  vector<CFreal> spectralLoopPowersLocal;
  vector<CFreal> cpdEmsLocal;
  if (m_saveMemory) {
    nbCells = dataMat.nbRows(); // only local data
    const CFuint nbCellsPerProc = totalNbCells/m_nbProc;
    for (CFuint rank = 0; rank < m_rank; ++rank) {
      offsetStateID += nbCellsPerProc;
    }
    spectralLoopPowersLocal.resize(nbCells);
    spectralLoopPowersCurr = &spectralLoopPowersLocal[0];
    cpdEmsLocal.resize(nbCells * nbCpdPoints);
    cpdEmsCurr = &cpdEmsLocal[0];
  }
  
  CFreal emInt = 0.;
  for (CFuint s = 0; s < nbCells; ++s) {
    //first pass to get the total emission coefficient
    //use this information to get the cell's total Radiative Power
    emInt=0.;
    cpdEmsCurr[s*nbCpdPoints] = 0;
    for (CFuint j=1; j < nbCpdPoints; j++ ) {
      //cout<<"spectral point "<<j<<" of "<<nbCpdPoints<<endl;
      emInt += ( dataMat(s,(j-1)*3+1) + dataMat(s,j*3+1) )/2. * m_dWav*1e-10;
      cpdEmsCurr[s*nbCpdPoints + j] = emInt;
    }
    
    const CFuint stateID = s+offsetStateID;
    spectralLoopPowersCurr[s] = emInt*getCellVolume( m_pstates->getStateLocalID(stateID)) * 
      m_angstrom*4.0* MathConsts::CFrealPi();
    cf_assert(getCellVolume( m_pstates->getStateLocalID(stateID)) > 0.);
    
    //second pass to get the [0->1] cumulative probably distribution
    cf_assert(emInt > 0.);
    for (CFuint j=1; j < nbCpdPoints; j++ ){
      cpdEmsCurr[ s*nbCpdPoints + j ] /= emInt;
    }
  }
  
  if (m_saveMemory) {
    CFuint minSizeToSend = 0;
    CFuint maxSizeToSend = 0;
    vector<int> recvCounts(m_nbProc, 0);
    vector<int> displs(m_nbProc, 0);
    computeRecvCountsDispls(totalNbCells, 1, minSizeToSend, maxSizeToSend, recvCounts, displs);
    CFuint sendSize = nbCells;
    cf_assert(sendSize <= maxSizeToSend);
    cf_assert(sendSize >= minSizeToSend);
    cf_assert(sendSize <= maxSizeToSend);
    cf_assert(sendSize >= minSizeToSend);
    
    MPIError::getInstance().check
      ("MPI_Allgatherv", "ParadeRadiator::computeEmissionCPD() => m_spectralLoopPowers",
       MPI_Allgatherv(&spectralLoopPowersCurr[0], sendSize, MPIStructDef::getMPIType(&spectralLoopPowersCurr[0]),
		      &m_spectralLoopPowers[0], &recvCounts[0], &displs[0],  MPIStructDef::getMPIType(&m_spectralLoopPowers[0]),
		      PE::GetPE().GetCommunicator(m_namespace)));
    
    computeRecvCountsDispls(totalNbCells, nbCpdPoints, minSizeToSend, maxSizeToSend, recvCounts, displs);
    sendSize = nbCells*nbCpdPoints;
    cf_assert(sendSize <= maxSizeToSend);
    cf_assert(sendSize >= minSizeToSend);
    cf_assert(sendSize <= maxSizeToSend);
    cf_assert(sendSize >= minSizeToSend);
    
    MPIError::getInstance().check
      ("MPI_Allgatherv", "ParadeRadiator::computeEmissionCPD() => m_cpdEms",
       MPI_Allgatherv(&cpdEmsCurr[0], sendSize, MPIStructDef::getMPIType(&cpdEmsCurr[0]),
		      &m_cpdEms[0], &recvCounts[0], &displs[0],  MPIStructDef::getMPIType(&m_cpdEms[0]),
		      PE::GetPE().GetCommunicator(m_namespace))); 
  } 
  
/*  if (  PE::GetPE().GetRank()  == 0) {
  //test for the first cell
  cout<<endl<<" wav = [";
  for (CFuint j=0; j < m_nbPoints; j++ ){
    cout<<dataMat(0,j*3+0)<<' ';
  }
  cout<<" ];" <<endl<<" em = [";
  for (CFuint j=0; j < m_nbPoints; j++ ){
    cout<<dataMat(0,j*3+1)<<' ';
  }
  cout<<" ];"<<endl<<" am = [";
  for (CFuint j=0; j < m_nbPoints; j++ ){
    cout<<dataMat(0,j*3+2)<<' ';
  }

  cout<<" ];"<<endl<<" cmd = [ ";
  for (CFuint j=0; j < nbCpdPoints; j++ ){
    cout<<m_cpdEms[ j ] <<' ';
  }

  cout<<" ];"<<endl<<" integral= "<<m_spectralLoopPowers[0]/getCellVolume( m_pstates->getStateLocalID(0))<<' '<<endl;

  cout<<"state idx= "<<m_radPhysicsHandlerPtr->getCurrentCellTrsIdx()<<endl;
  cout<<"PDF = [";
  RealVector so(3);
  CFreal wav;
  for(CFuint i=0;i<1000;++i){
    for(CFuint j=0;j<10;++j){
      getRandomEmission(wav,so);
      cout<< wav <<' ';
    }
    cout<<" ..."<<endl;
  }
  cout<<" ];"<<endl;
}
*/
  CFLog(VERBOSE, "ParadeRadiator::computeEmissionCPD() => end\n");
}

//////////////////////////////////////////////////////////////////////////////

CFreal ParadeRadiator::getSpectraLoopPower()
{
  return m_spectralLoopPowers[ m_radPhysicsHandlerPtr->getCurrentCellTrsIdx() ];
}

//////////////////////////////////////////////////////////////////////////////

void ParadeRadiator::getRandomEmission(CFreal &lambda, RealVector &s_o)
{
  //cout<<"get emission"<<endl;
  static CFuint dim = Framework::PhysicalModelStack::getActive()->getDim();
  static CFuint dim2 = m_radPhysicsHandlerPtr->isAxi() ? 3 : dim;
  
  CFuint nbCpdPoints = m_nbPoints;
  CFuint stateIdx = m_radPhysicsHandlerPtr->getCurrentCellTrsIdx();
  
  //cout<<"cpd, stateIdx "<<stateIdx<<" :";
  //vector<CFreal>::iterator it;
  //for(it = it_start; it<=it_end; ++it ){
    //cout<<*it<<' ';
  //}
  //cout<<endl;
  
  CFreal rand = m_rand.uniformRand();
  CFreal* it_start = &m_cpdEms[ stateIdx*nbCpdPoints ];
  CFreal* it_end = &m_cpdEms[ (stateIdx+1)*nbCpdPoints-1 ];
  CFreal* it_upp = std::upper_bound(it_start, it_end, rand);
  CFuint spectralIdx1 = it_upp - it_start;
  CFuint spectralIdx2 = spectralIdx1 - 1;
  
  //cout<<stateIdx<<dx<<' '<<spectralIdx1<<' '<<spectralIdx2<<endl;

  const CFuint nbCols = m_nbPoints*3;
  const CFreal x0 = m_cpdEms[ stateIdx*nbCpdPoints + spectralIdx1 ];
  const CFreal y0 = m_data[stateIdx*nbCols + spectralIdx1*3];
  const CFreal x1 = m_cpdEms[ stateIdx*nbCpdPoints + spectralIdx2 ];
  const CFreal y1 = m_data[stateIdx*nbCols + spectralIdx2*3];

  lambda = y0 + (y1-y0) * (rand - x0) / (x1-x0);

  //cout<<x0<<' '<<x1<<' '<<y0<<' '<<y1<<' '<<rand<<' '<<lambda<<endl;

  m_rand.sphereDirections(dim2, s_o);
}

//////////////////////////////////////////////////////////////////////////////

void ParadeRadiator::computeBinning()
{    
  CFLog(VERBOSE, "ParadeRadiator::computeBinning() => START\n");

  Stopwatch<WallTime> stp;
  stp.start();
  
  const CFuint nbCols = m_nbPoints*3;
  RealMatrix dataMat(m_data.size()/nbCols, nbCols, &m_data[0]);
  
  const CFuint totalNbCells = m_pstates->getSize();
  CFuint nbCells = dataMat.nbRows();
  cf_assert((!m_saveMemory && nbCells == totalNbCells) || 
	    ( m_saveMemory && nbCells <= totalNbCells));

  if (m_saveMemory) {
    CFuint totalNbCells = 0;
    MPIError::getInstance().check
      ("MPI_Allreduce", "ParadeRadiator::computeBinning()",
       MPI_Allreduce(&nbCells, &totalNbCells, 1, MPIStructDef::getMPIType(&nbCells), 
		     MPI_SUM, PE::GetPE().GetCommunicator(m_namespace)));
    cf_assert(totalNbCells == m_pstates->getSize());
  }

  /*for(CFuint i=0;i<m_nbPoints;++i) {
    for(CFuint j=0;j<nbCells;++j) {
      CFLog(VERBOSE, "m_data( " << j << "," << i << ") = " << dataMat(j,i*3+2) << "\n");
      CFLog(VERBOSE, "m_data( " << j << "," << i << ") = " << dataMat(j,i*3+1) << "\n");
    }
  }*/

  //Function to calculate the binning algorithm with absorption & emission coefficients
  const CFuint nbBins = m_nbBins;
  CFreal num_alphatot = 0.0;
  CFreal den_alphatot = 0.0;
  m_alphaav.resize(m_nbPoints, 0.);
  
  vector<CFreal>* num_alphatotVec = CFNULL;
  vector<CFreal>* den_alphatotVec = CFNULL;
  if (m_saveMemory) {
    num_alphatotVec = new vector<CFreal>(m_nbPoints, 0.);
    den_alphatotVec = new vector<CFreal>(m_nbPoints, 0.);
  }
  
  // offset for the ID from which the cells need to be counted 
  CFuint offsetStateID = 0;
  if  (m_saveMemory) {
    const CFuint nbCellsPerProc = totalNbCells/m_nbProc;
    for (CFuint rank = 0; rank < m_rank; ++rank) {
      offsetStateID += nbCellsPerProc;
    }
  }
  
  const CFreal h = 6.626070040e-34;     //SI units Js
  const CFuint c = 3e08; //SI units m/s
  const CFreal k_b = 1.3806485279e-23;  //SI units J/K
  const CFuint tempID = m_radPhysicsHandlerPtr->getTempID();
  DeviceFunc df;
  
  for(CFuint i=0;i<m_nbPoints;++i) {
    num_alphatot = 0.;
    den_alphatot = 0.;
    CFreal Bs = 0.;
    for(CFuint s=0;s<nbCells;++s) {
      const CFreal alpha = dataMat(s,i*3+2);
      const CFreal epsilon = dataMat(s,i*3+1);
      CFreal *const currState = m_pstates->getState(s);
      CFreal temp  = currState[tempID];

      if(!m_Equilibrium){
	if(std::abs(alpha) > 0.){
          Bs = epsilon/alpha;
    	}
	else{
	  Bs = df.computePlanckFun(h,c,k_b,dataMat(s,i*3),temp);
	  //(2*h*pow(c,2)/pow((dataMat(s,i*3)*1e-10),5))*
	  //(1/(exp(h*c/((dataMat(s,i*3)*1e-10)*k_b*temp))-1));
	}
      }
      else{
	// Planck function if we are in equilibrium
	//
	Bs = df.computePlanckFun(h,c,k_b,dataMat(s,i*3),temp);
	CFLog(DEBUG_MIN,"The source term computed with Planck function for temperature " << temp << "and wavelength " 
	      << dataMat(s,i*3) << " is = " << Bs << "\n");
      }
      
      const CFuint stateID = s+offsetStateID;
      const CFreal volume = getCellVolume(m_pstates->getStateLocalID(stateID));
      cf_assert(volume > 0.); 
      const CFreal BsVolume = Bs*volume;
      const CFreal num_alpha_vol = alpha*BsVolume;
      const CFreal den_alpha_vol = BsVolume;
      num_alphatot += num_alpha_vol;
      den_alphatot += den_alpha_vol;
    }
    
    if (!m_saveMemory) {
      m_alphaav[i] = num_alphatot / den_alphatot;
    }
    else {
      (*num_alphatotVec)[i] = num_alphatot;
      (*den_alphatotVec)[i] = den_alphatot;
    }
    CFLog(DEBUG_MIN,"ParadeLibrary::computeBinning () => m_alphaav = " << m_alphaav[i] <<"\n");
  }
  
  if (m_saveMemory) {
    // compute the total numerator and denominator per spectral point across all ranks
    vector<CFreal> num_alphatotGlobal(m_nbPoints, 0.);
    vector<CFreal> den_alphatotGlobal(m_nbPoints, 0.);
    MPIError::getInstance().check
      ("MPI_Allreduce", "ParadeRadiator::computeBinning()",
       MPI_Allreduce(&(*num_alphatotVec)[0], &num_alphatotGlobal[0], m_nbPoints, 
		     MPIStructDef::getMPIType(&(*num_alphatotVec)[0]), 
		     MPI_SUM, PE::GetPE().GetCommunicator(m_namespace)));
    
    MPIError::getInstance().check
      ("MPI_Allreduce", "ParadeRadiator::computeBinning()",
       MPI_Allreduce(&(*den_alphatotVec)[0], &den_alphatotGlobal[0], m_nbPoints, 
		     MPIStructDef::getMPIType(&(*den_alphatotVec)[0]), 
		     MPI_SUM, PE::GetPE().GetCommunicator(m_namespace)));
    
    cf_assert(m_alphaav.size() == m_nbPoints);
    for (CFuint i = 0; i < m_alphaav.size(); ++i) {
      cf_assert(std::abs(den_alphatotGlobal[i]) > 0.);
      m_alphaav[i] = num_alphatotGlobal[i] / den_alphatotGlobal[i];
      // CFLog(INFO, "m_alphaav[" << i << "] = " <<  m_alphaav[i] << "\n");
    }
  } 
  
  // free memory
  deletePtr(num_alphatotVec);
  deletePtr(den_alphatotVec);
    
  CFreal alphamin = m_alphaav[0];
  CFreal alphamax = m_alphaav[0];
  for(CFuint r=0;r<m_nbPoints;++r) {
    if(m_alphaav[r]<alphamin) {
      alphamin = m_alphaav[r];
    }
    if(m_alphaav[r]>alphamax) {
      alphamax = m_alphaav[r];
    }
  }
  
  if(alphamin == 0){alphamin = 1.e-30;}
  if(alphamax == 0){alphamax = 1.e-30;}  

  CFLog(INFO,"ParadeLibrary::computeBinning () => [alphamin, alphamax] = [" 
	<< alphamin <<", " << alphamax << "]\n");
  
  /* string alphaFile = "alpha.txt_" + StringOps::to_str(m_rank); 
  ofstream fout(alphaFile.c_str());
  for (CFuint i = 0; i < m_alphaav.size(); ++i) {
  fout << m_alphaav[i] << endl;
  }*/
    
  // Logarithmic spacing
  //
  const CFreal alpha_minlog = std::log(alphamin);
  const CFreal alpha_maxlog = std::log(alphamax);
  
  CFLog(INFO,"ParadeLibrary::computeBinning () => [alpha_minlog, alpha_maxlog] = [" 
	<< alpha_minlog << ", " << alpha_maxlog << "]\n");
  
  const CFreal dy = (alpha_maxlog-alpha_minlog) / (nbBins-1);
  CFLog(VERBOSE,"ParadeLibrary::computeBinning () => dy = " << dy <<"\n");

#ifdef CF_HAVE_CUDA
  LocalArray<CFreal>::TYPE vctBins(0., nbBins);
#else
  LocalArray<CFreal>::TYPE vctBins(nbBins, 0.);
#endif
  cf_assert(vctBins.size() == nbBins);  
  for(CFuint i = 0; i<nbBins; ++i) {
    vctBins[i] = std::exp(alpha_minlog + (dy * i));
    CFLog(VERBOSE,"ParadeLibrary::computeBinning () => vctBins(" << i << ") = " << vctBins[i] <<"\n");
  }
  
  // To search for minimum and maximum alpha
  CFreal alphamin_tot = m_alphaav[0];
  CFreal alphamax_tot = m_alphaav[0];
  for(CFuint i=0;i<m_nbPoints;++i) {
    for(CFuint s=0;s<nbCells;s++) {
      const CFreal alphaA = dataMat(s,i*3+2);
      if(alphaA < alphamin_tot) {alphamin_tot = alphaA;}
      if(alphaA > alphamax_tot) {alphamax_tot = alphaA;}
    }
  }

  if (m_saveMemory) {
    CFreal alphaMinTotLocal = alphamin_tot;
    CFreal alphaMaxTotLocal = alphamax_tot;
    CFreal alphaMinTotGlobal = 0.;
    CFreal alphaMaxTotGlobal = 0.;
    MPIError::getInstance().check
      ("MPI_Allreduce", "ParadeRadiator::computeBinning()",
       MPI_Allreduce(&alphaMinTotLocal, &alphaMinTotGlobal, 1, 
		     MPIStructDef::getMPIType(&alphaMinTotLocal), 
		     MPI_MIN, PE::GetPE().GetCommunicator(m_namespace)));
    MPIError::getInstance().check
      ("MPI_Allreduce", "ParadeRadiator::computeBinning()",
       MPI_Allreduce(&alphaMaxTotLocal, &alphaMaxTotGlobal, 1, 
		     MPIStructDef::getMPIType(&alphaMaxTotLocal), 
		     MPI_MAX, PE::GetPE().GetCommunicator(m_namespace)));
    
    alphamin_tot = alphaMinTotGlobal;
    alphamax_tot = alphaMaxTotGlobal;
  }

  if(alphamin_tot == 0){alphamin_tot == 1e-30;}
  if(alphamax_tot == 0){alphamax_tot == 1e-30;}
  
  CFLog(INFO,"ParadeLibrary::computeBinning () => [alphamin_tot, alphamax_tot] = [" 
	<< alphamin_tot << ", " << alphamax_tot << "]\n");
  
  vctBins[0] = alphamin_tot;
  vctBins[nbBins-1] = alphamax_tot;

  // AL: up to here it is correct
  /*string vctFile = "vct.txt_" + StringOps::to_str(m_rank); 
    ofstream fout(vctFile.c_str());
    for (CFuint i = 0; i < vctBins.size(); ++i) {
    fout << vctBins[i] << endl;
  }*/
  
  computeAveragedBins(nbBins,2, vctBins);
    
  /*
    const CFuint tempID2 = m_radPhysicsHandlerPtr->getTempID();
    CFreal *const currState2 = m_pstates->getState(24);
    CFreal temp2  = currState2[tempID2];
    CFLog(INFO, "The temperature for cell 25 is = " << temp2 << "\n");
    CFreal sigma_check = 5.67e-08;
    CFreal Planckintegral = sigma_check*pow(temp2,4);
    CFLog(INFO, "The Plank integral value for cell 25 is = " << Planckintegral << "\n");
    CFreal Totalsource = 0.0;
    for(CFuint w =1; w<nbBinsre;w++){
    CFLog(VERBOSE, "The source term for cell 25 and band " << w << " is = " <<  B_bin[w+nbBinsre*24] << "\n");
    Totalsource += B_bin[w+nbBinsre*24];
    }
    CFLog(INFO, "The Planck integral value banded for cell 28 is = " << Totalsource << "\n");
  */
  
  /*
    ofstream fout1("inwav.txt"); 
    fout1 << "TITLE = Original spectrum of radiative properties\n";
     fout1 << "VARIABLES = Wavl EmCoef AbCoef \n";
     for (CFuint i = 0; i < m_data.nbCols()/3; ++i) {
     fout1 << dataMat(0,i*3) << " " << dataMat(0,i*3+1) << " " << dataMat(0,i*3+2) << endl;
     }
     fout1.close();
     
     ofstream fout3("notbinnedabs.txt"); 
     fout3 << "TITLE = Original spectrum of radiative properties\n";
     fout3 << "VARIABLES = AbCoef \n";
     for (CFuint i = 0; i < m_data.nbCols()/3; ++i) {
     fout3 << dataMat(0,i*3+2) << " " << endl;
     }
     fout3.close();

     ofstream fout2("binnedabs.txt"); 
     fout2 << "TITLE = Binned spectrum of radiative properties\n";
     fout2 << "VARIABLES = Wavl EmCoef AbCoef \n";
     for (CFuint k = 1; k < nbBinsre; ++k) {
       fout2 << alpha_avbin[k] << " " << endl;
     }
     fout2.close();
   
*/
  CFLog(INFO, "ParadeRadiator::computeBinning() took " << stp.read() << "s\n");
  CFLog(VERBOSE, "ParadeRadiator::computeBinning() => END\n");
}

//////////////////////////////////////////////////////////////////////////////
  
void ParadeRadiator::computeBanding() 
{  
  CFLog(VERBOSE, "ParadeRadiator::computeBanding() => START\n");
  
  const CFuint Bandcoeff = (m_nbBands + 1);
  
  Stopwatch<WallTime> stp;
  stp.start();
  
  const CFuint nbCols = m_nbPoints*3;
  RealMatrix dataMat(m_data.size()/nbCols, nbCols, &m_data[0]);
  
  const CFuint totalNbCells = m_pstates->getSize();
  CFuint nbCells = dataMat.nbRows();
  cf_assert((!m_saveMemory && nbCells == totalNbCells) || 
	    ( m_saveMemory && nbCells <= totalNbCells));
  
  if (m_saveMemory) {
    CFuint totalNbCells = 0;
    MPIError::getInstance().check
      ("MPI_Allreduce", "ParadeRadiator::computeBinning()",
       MPI_Allreduce(&nbCells, &totalNbCells, 1, MPIStructDef::getMPIType(&nbCells), 
		     MPI_SUM, PE::GetPE().GetCommunicator(m_namespace)));
    cf_assert(totalNbCells == m_pstates->getSize());
  }

  // Banding
  CFreal alphamin_tot = dataMat(0,0);
  CFreal alphamax_tot = dataMat(0,0);
  
  /// The following lines are to check everything
  CFLog(INFO,"Number of spectral points is: " << m_nbPoints << "\n");
  CFLog(INFO,"Number of cells is: " << nbCells << "\n");
  
  for(CFuint i=0;i<m_nbPoints;++i) {
    for(CFuint s=0;s<nbCells;s++) {
      const CFreal wavelength = dataMat(s,i*3);
      if(wavelength < alphamin_tot) {
	alphamin_tot = wavelength;
      }
      if(wavelength > alphamax_tot) {
	alphamax_tot = wavelength;
      }
    }
  }

  CFreal alphamax_totGlobal = 0.;
  CFreal alphamin_totGlobal = 1e10;

 MPIError::getInstance().check
    ("MPI_Allreduce", "ParadeRadiator::computeBanding()",
     MPI_Allreduce(&alphamin_tot, &alphamin_totGlobal, 1, 
		   MPIStructDef::getMPIType(&alphamin_tot), 
		   MPI_MIN, PE::GetPE().GetCommunicator(m_namespace)));
  
  MPIError::getInstance().check
    ("MPI_Allreduce", "ParadeRadiator::computeBanding()",
     MPI_Allreduce(&alphamax_tot, &alphamax_totGlobal, 1, 
		   MPIStructDef::getMPIType(&alphamax_tot), 
		   MPI_MAX, PE::GetPE().GetCommunicator(m_namespace)));
  
  const CFreal alpha_minlog = std::log(alphamin_totGlobal);
  const CFreal alpha_maxlog = std::log(alphamax_totGlobal);
  
  CFLog(VERBOSE,"ParadeLibrary::computeBanding () => alpha_minlog = " << alpha_minlog <<"\n");
  CFLog(VERBOSE,"ParadeLibrary::computeBanding () => alpha_maxlog = " << alpha_maxlog <<"\n");

  const CFreal dy = (alpha_maxlog-alpha_minlog) / (m_nbBands);
  
  CFLog(INFO,"ParadeLibrary::computeBanding() => dy = " << dy <<"\n");

#ifdef CF_HAVE_CUDA
  LocalArray<CFreal>::TYPE vctBins(0., Bandcoeff);
#else
  LocalArray<CFreal>::TYPE vctBins(Bandcoeff, 0.);
#endif
  for(int i = 0; i<Bandcoeff; ++i) {
    vctBins[i] = std::exp(alpha_minlog + (dy * i));
    CFLog(VERBOSE,"ParadeLibrary::computeBanding() => vctBins(" << i << ") = " << vctBins[i] <<"\n");
  }
  
   computeAveragedBins(m_nbBands, 0, vctBins);

  // To be commented out after verification
  if(PE::GetPE().GetRank(m_namespace) == 0) {
    DataHandle<CFreal> alpha_avbin = m_radPhysicsHandlerPtr->getDataSockets()->alpha_avbin;
    DataHandle<CFreal> B_bin = m_radPhysicsHandlerPtr->getDataSockets()->B_bin;  

    ofstream fout1("alpha.txt");
    for(CFuint j=0;j<totalNbCells;++j) {
      for(CFuint k=1;k<m_nbBands;++k) {
	fout1 << "alpha (" << k << "," << j << ") = " << 
	  alpha_avbin[k + m_nbBands*j] << "\n";
      }
    }
    fout1.close();
    
    ofstream fout2("beta.txt");
    for(CFuint j=0;j<totalNbCells;++j) {
      for(CFuint k=1;k<m_nbBands;++k) {
	fout2 << "beta (" << k << "," << j << ") = " << 
	  B_bin[k + m_nbBands*j] << "\n";
      }
    }
    fout2.close();
  }
  
  //INTEGRAL OF THE PLANCK FUNCTION
  /*
  const CFuint tempID2 = m_radPhysicsHandlerPtr->getTempID();
      CFreal *const currState2 = m_pstates->getState(29);
      CFreal temp2  = currState2[tempID2];

      CFLog(INFO, "The temperature for cell 29 is = " << temp2 << "\n");

      CFreal sigma_check = 5.67e-08;
      CFreal Planckintegral = sigma_check*pow(temp2,4);
      CFLog(INFO, "The Plank integral value for cell 29 is = " << Planckintegral << "\n");

      CFreal Totalsource = 0.0;
       for(CFuint w =1; w<m_nbBands;w++){
       CFLog(VERBOSE, "The source term for cell 24 and band " << w << " is = " <<  B_bin[w+m_nbBands*29] << "\n");
      	Totalsource += B_bin[w+m_nbBands*29];
       }

      CFLog(INFO, "The Planck integral value banded for cell 29 is = " << Totalsource << "\n");


	CFLog(INFO, "The source term for cell 29 and band 99  is = " <<  B_bin[99+m_nbBands*29] << "\n");
  */




  CFLog(INFO, "ParadeRadiator::computeBanding() took " << stp.read() << "s\n");
  CFLog(VERBOSE, "ParadeRadiator::computeBanding() => END\n");
}

/////////////////////////////////////////////////////////////////////////////

void ParadeRadiator::computeBinningBanding() 
{ 
  CFLog(VERBOSE, "ParadeRadiator::computeBinningBanding() => START\n");
 
  Stopwatch<WallTime> stp;
  stp.start();
  const CFuint nbCells = m_pstates->getSize();
  
  CFuint nbBins = m_nbBins; //for each band, to be set in the *.CFcase
  RealVector vctBands(m_nbBands + 1); //vector containing the limits of bands
  
  if (nbBins == 1) {
    CFLog(INFO, "ParadeLibrary::computeproprieties() => nbBins = 1 is not allowed, nbBins = 2 is set\n"); //otherwise dy = inf
    nbBins = 2;
  }
  
  //variables for binning
  const CFuint nbBinsre = nbBins;
  const CFuint nbBandsre = m_nbBands;
  m_alpha_bin.resize(nbBinsre*nbCells);
  m_emission_bin.resize(nbBinsre*nbCells);
  m_B_bin.resize(nbBinsre*nbCells*nbBandsre);       //bands included as bins
  m_alpha_avbin.resize(nbBinsre*nbCells*nbBandsre); //bands included as bins
  RealVector vctBins(nbBins);
  CFreal alpha_bincoeff = 0.;
  CFreal emission_bincoeff = 0.;
  CFreal B_bincoeff = 0.;
  
  const CFuint nbCols = m_nbPoints*3;
  RealMatrix dataMat(m_data.size()/nbCols, nbCols, &m_data[0]);
  
  DataHandle<CFreal> alpha_avbin = m_radPhysicsHandlerPtr->getDataSockets()->alpha_avbin;
  DataHandle<CFreal> B_bin = m_radPhysicsHandlerPtr->getDataSockets()->B_bin;
  
  //logarithmic spacing for banding
  if (m_bandsDistr == "logarithmic") {
    CFLog(VERBOSE,"ParadeLibrary::computeproperties () => Logarithmic banding\n");
    const CFreal wavMin_log = std::log(m_wavMin);
    const CFreal wavMax_log = std::log(m_wavMax);
    const CFreal dx = (wavMax_log-wavMin_log)/m_nbBands;
    vctBands[0] = m_wavMin;
    for(CFuint i = 1; i<m_nbBands; ++i) {
      vctBands[i] = std::exp(wavMin_log + (dx * i));
    }
    vctBands[m_nbBands] = m_wavMax;
  }
  
  //equally-spaced distribution for banding
  else if(m_bandsDistr == "equally") {
    CFLog(VERBOSE,"ParadeLibrary::computeproperties () => Equally-spaced banding\n");
    CFreal dx = (m_wavMax-m_wavMin)/m_nbBands;
    vctBands[0] = m_wavMin;
    for(CFuint i = 1; i<m_nbBands; ++i) {
      vctBands[i] = m_wavMin + (dx * i);
    }
    vctBands[m_nbBands] = m_wavMax;
  }
  
  for(CFuint i = 0; i<(m_nbBands + 1); ++i) {
    CFLog(INFO,"ParadeLibrary::computeproperties () => vctBands(" << i << ") = " << vctBands[i] <<"\n");
  }
  
  //Number of points in each band for each cell
  //initialized to 1 to take into accout the starting point of each band
  std::vector< int > nb_bandPoints(m_nbBands, 0); 
  
  for(CFuint ii=0; ii<m_nbPoints; ++ii) {
    for(CFuint iii=0; iii<m_nbBands; ++iii) {
      if(dataMat(1, ii*3)>vctBands[iii] && dataMat(1, ii*3)<vctBands[iii+1]) {
	//counting points within the band
	//on the 1st cell 'cause the overall points distribution does not change over cells
	nb_bandPoints[iii]++;
      }
    }
  }
  nb_bandPoints[0]++;
  nb_bandPoints[m_nbBands-1]++; //To include wavMax in the last band
  
  for(CFuint iii=0; iii<m_nbBands; ++iii) {
    CFLog(INFO,"ParadeLibrary::computeproperties () => nb_bandPoints(" << iii  << ") = " << nb_bandPoints[iii] << "\n");
  }
  
  CFuint cc = 0; //counter to pick the proper proprieties from m_data
  for(CFuint b=0; b<m_nbBands; ++b) { //the loop over bands starts
    //Function to calculate the binning algorithm with the absorption and emission coefficients
    CFreal num_alphatot = 0.0;
    CFreal den_alphatot = 0.0;
    m_alphaav.clear();
    
    for(CFuint i=0;i<nb_bandPoints[b]; ++i) {
      for(CFuint s=0;s<nbCells; ++s) {
	const CFreal alpha = dataMat(s,(i+cc)*3+2);
	const CFreal epsilon = dataMat(s,(i+cc)*3+1);
	const CFreal Bs = epsilon/alpha;
	
	// Planck function if we are not in equilibrium
	//
	// CFreal h = 6.626070040e-34;     //SI units Js
	// CFuint c = 3e08; //SI units m/s
	// CFreal k_b = 1.3806485279e-23;  //SI units J/K
	//
	//
	// T needs to be defined from the values extracted
	// const  B = (2*h*pow(c,2)/pow(dataMat(j,i*3),5))*(1/(exp(h*c/(dataMat(j,i*3)*k_b*T))-1));	
	
	const CFreal volume = getCellVolume( m_pstates->getStateLocalID(s));
	cf_assert(volume > 0.); 
	const CFreal BsVolume = Bs*volume;
	const CFreal num_alpha_vol = alpha*BsVolume;
	const CFreal den_alpha_vol = BsVolume;
	num_alphatot += num_alpha_vol;
	den_alphatot += den_alpha_vol;
      }
      
      m_alphaav.push_back(num_alphatot/den_alphatot);
      CFLog(VERBOSE,"ParadeLibrary::computeproperties () => m_alphaav = " << m_alphaav[i] <<"\n");
    }
    
    CFreal alphamin = m_alphaav[0];
    CFreal alphamax = m_alphaav[0];
    
    for(CFuint r=0; r<nb_bandPoints[b]; ++r) {
      if(m_alphaav[r]<alphamin) {
	alphamin = m_alphaav[r];
      }
      if(m_alphaav[r]>alphamax) {
	alphamax = m_alphaav[r];
      }
    }
    
    CFLog(INFO,"ParadeLibrary::computeBinningBanding () => [alphamin, alphamax] = [" 
	  << alphamin <<", " << alphamax << "]\n");
    
    // Logarithmic spacing for binning
    //
    
    CFreal alpha_minlog = std::log(alphamin);
    CFreal alpha_maxlog = std::log(alphamax);
    
    CFLog(VERBOSE,"ParadeLibrary::computeproperties () => alpha_minlog = " << alpha_minlog <<"\n");
    CFLog(VERBOSE,"ParadeLibrary::computeproperties () => alpha_maxlog = " << alpha_maxlog <<"\n");
    
    const CFreal dy = (alpha_maxlog-alpha_minlog) / (nbBins-1);
    CFLog(INFO,"ParadeLibrary::computeproperties () => dy = " << dy <<"\n");
    
    vctBins = 0.;
    for(int i = 0; i<nbBins; ++i) {
      vctBins[i] = std::exp(alpha_minlog + (dy * i));
      CFLog(VERBOSE,"ParadeLibrary::computeproperties () => vctBins(" << i << ") = " << vctBins[i] <<"\n");
    }
    
    //To search for minimum alpha and maximum
    
    CFreal alphamin_tot = m_alphaav[0];
    CFreal alphamax_tot = m_alphaav[0];
    
    for(CFuint i=0;i<nb_bandPoints[b]; ++i) {
      for(CFuint s=0;s<nbCells; ++s) {
	if((dataMat(s,(i+cc)*3+2))<alphamin_tot) {
	  alphamin_tot = dataMat(s,(i+cc)*3+2);
	}
	if((dataMat(s,(i+cc)*3+2))>alphamax_tot) {
	  alphamax_tot = dataMat(s,(i+cc)*3+2);
	}
      }
    }
    
    vctBins[0] = alphamin_tot;
    vctBins[nbBins-1] = alphamax_tot;
    
    for( CFuint i = 0; i<nbBins; ++i) {
      CFLog(INFO, " vctBins = " << vctBins[i] << "\n");
    }
    
    // Binning part of the function
    for(CFuint i=0;i<nb_bandPoints[b]; ++i) {
      for(CFuint j=0;j<nbCells; ++j){
	CFreal Bs = dataMat(j,(i+cc)*3+1)/dataMat(j,(i+cc)*3+2);
	
	m_alpha_bin[0+nbBinsre*j] = 0;
	m_emission_bin[0+nbBinsre*j] = 0;
	B_bin[0+nbBinsre*j] = 0;
	
	for(CFuint k=1;k<nbBinsre; ++k) {
	  if((dataMat(j,(i+cc)*3+2)>=vctBins[k-1]) && (dataMat(j,(i+cc)*3+2)<vctBins[k])) {
	    alpha_bincoeff = dataMat(j,(i+cc)*3+2)*Bs*m_dWav*1.e-10;
	    emission_bincoeff = dataMat(j,(i+cc)*3+1)*m_dWav*1e-10; // AL: here there was no 1e10
	    B_bincoeff = Bs*m_dWav*1.e-10;
	    /*
	      CFLog(VERBOSE,"ParadeLibrary::computeproperties () => m_emission_bin(" 
	      << k << "," << j << ") = " << m_emission_bin[nbBins*j+k] <<"\n");
	      CFLog(VERBOSE,"ParadeLibrary::computeproperties () => m_alpha_bin(" 
	      << k << "," << j << ") = " <<  m_alpha_bin[nbBins*j+k] <<"\n");
	      CFLog(VERBOSE,"ParadeLibrary::computeproperties () => m_B_bin(" 
	      << k << "," << j << ") = " << m_B_bin[nbBins*j+k] <<"\n");
	    */
	  }

	  else{
	    alpha_bincoeff = 0.;
	    emission_bincoeff = 0.;
	    B_bincoeff = 0.;
	  }
	  
	  m_alpha_bin[k + nbBinsre*j] += alpha_bincoeff;
	  m_emission_bin[k + nbBinsre*j] += emission_bincoeff;
	  m_B_bin[k + nbBinsre*j] += B_bincoeff;
	}
      }
    }
    
    for(CFuint k=1;k<nbBinsre; ++k) {
      for(CFuint j=0;j<nbCells; ++j) {
	CFLog(VERBOSE,"ParadeLibrary::computeproperties () => m_alpha_bin(" << k << "," << j << ") = " <<  m_alpha_bin[nbBinsre*j+k] <<"\n");
	CFLog(VERBOSE,"ParadeLibrary::computeproperties () => m_B_bin(" << k << "," << j << ") = " << m_B_bin[nbBinsre*j+k] <<"\n");
      }
    }
    
    for(CFuint j=0;j<nbCells; ++j) {
      for(CFuint k=1;k<nbBinsre;++k) {
	m_alpha_avbin[0+nbBinsre*j] = 0.;
	
	if(m_B_bin[k + nbBinsre*j] != 0.) {
	  CFreal m_alpha_avbin_a = (m_alpha_bin[k + nbBinsre*j]) / (m_B_bin[k + nbBinsre*j]);
	  m_alpha_avbin[k + nbBinsre*j] = m_alpha_avbin_a;
	  CFLog(VERBOSE,"ParadeRadiator::computeProperties()=>alpha_avbin(" << k << "," << j << ") = "<<m_alpha_avbin[k + nbBinsre*j]<<"\n");
	}
	else {
	  alpha_avbin[k + nbBinsre*j] = 0.;
	  CFLog(VERBOSE,"ParadeRadiator::computeProperties()=>alpha_avbin(" << k << ","<< j << ") = "<<alpha_avbin[k + nbBinsre*j]<<"\n");
	}
      }
    }
    
    //alpha_avbin, emission_bin are the average values for absorptivity and emissivity for each bin.
    
    for(CFuint j=0;j<nbCells; ++j) {
      for(CFuint k=1;k<nbBinsre; ++k) {
	//the index "b" added for banding
	alpha_avbin[(k+b) + nbBinsre*j] = m_alpha_avbin[k + nbBinsre*j];
	B_bin[(k+b) + nbBinsre*j] = m_B_bin[k + nbBinsre*j];
	CFLog(VERBOSE,"alpha (" << b << "," << k << "," << j << ") = " << alpha_avbin[(k+b) + nbBinsre*j] << "\n");
      }
    }
    cc = cc + nb_bandPoints[b];
  }
  
  CFLog(INFO, "ParadeRadiator::computeBinningBanding() took "<< stp.read() << "s\n");  CFLog(VERBOSE, "ParadeRadiator::computeBinningBanding() => END\n");
}
  
//////////////////////////////////////////////////////////////////////////////

void ParadeRadiator::computeRecvCountsDispls(const CFuint totalNbCells, const CFuint sizeCoeff, 
					     CFuint& minSizeToSend, CFuint& maxSizeToSend,
					     vector<int>& recvCounts, vector<int>& displs)
{
  const CFuint nbCellsPerProc = totalNbCells/m_nbProc;
  minSizeToSend = nbCellsPerProc*sizeCoeff;
  maxSizeToSend = minSizeToSend + (totalNbCells%m_nbProc)*sizeCoeff;
  CFuint count = 0;
  for (CFuint rank = 0; rank < m_nbProc; ++rank) {
    recvCounts[rank] = (rank < m_nbProc-1) ? minSizeToSend : maxSizeToSend;
    // displacements start at 0
    displs[rank] = count;
    count += minSizeToSend;
  }
}

//////////////////////////////////////////////////////////////////////////////

void ParadeRadiator::computeAveragedBins(const CFuint nbBinsre, 
					 const CFuint testID,
					 LocalArray<CFreal>::TYPE& vctBins)
{
  const CFuint totalNbCells = m_pstates->getSize();
  const CFuint nbCols = m_nbPoints*3;
  CFuint nbCells = m_data.size()/nbCols;
  
  //alpha_avbin is average value for absorptivity for each bin
  SafePtr<SocketBundle> sockets = m_radPhysicsHandlerPtr->getDataSockets();
  DataHandle<CFreal> alpha_avbin = sockets->alpha_avbin; // array w GLOBAL cell size
  DataHandle<CFreal> B_bin = sockets->B_bin;             // array w GLOBAL cell size
  
  if (alpha_avbin.size() != nbBinsre*totalNbCells) {
    CFLog(ERROR, "alpha_avbin.size() != nbBinsre*totalNbCells => " << 
	  alpha_avbin.size() << " != " <<  nbBinsre*totalNbCells << "\n");
    cf_assert(alpha_avbin.size() == nbBinsre*totalNbCells);
  }
  if (B_bin.size() != nbBinsre*totalNbCells) {
    CFLog(ERROR, "B_bin.size() != nbBinsre*totalNbCells => " << B_bin.size() << 
	  " != " <<  nbBinsre*totalNbCells << "\n");
    cf_assert(B_bin.size() == nbBinsre*totalNbCells);
  }
  
  // from now only local arrays (=global if m_saveMemory==false) are used
  m_alpha_bin.resize(nbBinsre*nbCells);     // array w LOCAL cell size
  m_emission_bin.resize(nbBinsre*nbCells);  // array w LOCAL cell size
  
  vector<CFreal> alpha_avbinLocal;
  vector<CFreal> B_binLocal;
  // by default (not optimized mode) you point to the global storage of bin arrays
  CFreal* alpha_avbinCurr = &alpha_avbin[0];
  CFreal* B_binCurr = &B_bin[0];
  
  if (m_saveMemory) {
    // if you run in optimized parallel mode, point to a local copy of the bin arrays
    // of size nbCells < totalNbCells
    cf_assert((nbCells < totalNbCells && m_nbProc > 1) || (nbCells == totalNbCells && m_nbProc == 1));
    
    alpha_avbinLocal.resize(nbBinsre*nbCells, 0.);
    alpha_avbinCurr = &alpha_avbinLocal[0];
    
    B_binLocal.resize(nbBinsre*nbCells, 0.);
    B_binCurr = &B_binLocal[0];
  }  

  const CFuint tempID = m_radPhysicsHandlerPtr->getTempID();
  DeviceFunc df;
  const CFuint isEquil = (const CFreal)m_Equilibrium;
  
  // fill the alpha and emission arrays for each local cell in parallel simulations
  for(CFuint i=0;i<m_nbPoints;++i) {
    // m_alpha_bin, m_emission_bin, B_bin are logically 2D arrays filled 
    // rows: each row gives all the bins for a cell
    for(CFuint j=0;j<nbCells;++j) {
      CFreal *const currState = m_pstates->getState(j);
      CFreal temp  = currState[tempID];
      df.computeCellBins<CPU>(isEquil, m_nbPoints, i,j, nbBinsre, testID, temp, m_dWav,
			      &vctBins[0], &m_data[0], &m_alpha_bin[0], &m_emission_bin[0],
			      B_binCurr);
    }
  }
  
  /*for(CFuint k=1;k<nbBinsre;++k) {
    for(CFuint j=0;j<nbCells;++j) {
    CFLog(DEBUG_MAX,"ParadeLibrary::computeproperties () => m_alpha_bin(" << k << "," << j << ") = " << m_alpha_bin[nbBinsre*j+k] <<"\n");
    CFLog(DEBUG_MAX,"ParadeLibrary::computeproperties () => B_binCurr(" << k << "," << j << ") = " << B_binCurr[nbBinsre*j+k] <<"\n");
    }
    }*/
  
  for(CFuint j=0;j<nbCells;++j) {
    for(CFuint k=1;k<nbBinsre;++k) {
      const CFuint idx0 = nbBinsre*j;
      alpha_avbinCurr[idx0] = 0.;
      const CFuint idx = k + idx0;
      // AL: is this fix needed to mask an error or it is supposed to be like this?
      if(B_binCurr[idx] != 0.) {
	alpha_avbinCurr[idx] = m_alpha_bin[idx] / B_binCurr[idx];
	CFLog(DEBUG_MED,"ParadeRadiator::computeAverageBins() => alpha_avbinCurr(" << k << "," << j << ") = "<< alpha_avbinCurr[idx] <<"\n");
      }
      else {
	alpha_avbinCurr[idx] = 0.;
	CFLog(DEBUG_MED,"ParadeRadiator::computeAverageBins() => alpha_avbinCurr(" << k << ","<< j << ") = "<< alpha_avbinCurr[idx] <<"\n");
      }
    }
  }
  
  if (m_saveMemory) {
    // here we need to gather all the entries for alpha_avbin and B_bin from all processes, so that
    // every process keeps a global storage of them
    CFuint minSizeToSend = 0;
    CFuint maxSizeToSend = 0;
    vector<int> recvCounts(m_nbProc, 0);
    vector<int> displs(m_nbProc, 0);
    computeRecvCountsDispls(totalNbCells, nbBinsre, minSizeToSend, maxSizeToSend, recvCounts, displs);
    const CFuint sendSize = nbCells*nbBinsre;
    cf_assert(sendSize <= maxSizeToSend);
    cf_assert(sendSize >= minSizeToSend);
    cf_assert(sendSize <= maxSizeToSend);
    cf_assert(sendSize >= minSizeToSend);
    
    MPIError::getInstance().check
      ("MPI_Allgatherv", "ParadeRadiator::computeAverageBins() => alpha_avbin",
       MPI_Allgatherv(&alpha_avbinCurr[0], sendSize, MPIStructDef::getMPIType(&alpha_avbinCurr[0]),
		      &alpha_avbin[0], &recvCounts[0], &displs[0],  MPIStructDef::getMPIType(&alpha_avbin[0]),
		      PE::GetPE().GetCommunicator(m_namespace)));
    
    MPIError::getInstance().check
      ("MPI_Allgatherv", "ParadeRadiator::computeAverageBins() => B_bin",
       MPI_Allgatherv(&B_binCurr[0], sendSize, MPIStructDef::getMPIType(&B_binCurr[0]),
		      &B_bin[0], &recvCounts[0], &displs[0],  MPIStructDef::getMPIType(&B_bin[0]),
		      PE::GetPE().GetCommunicator(m_namespace)));
  }
  
  // To be commented out after verification
  if(PE::GetPE().GetRank(m_namespace) == 0) {
   
    /*for(CFuint j=0;j<totalNbCells;++j) {
	for(CFuint k=1;k<nbBinsre;++k) {
	CFLog(DEBUG_MAX,"alpha (" << k << "," << j << ") = " << alpha_avbin[k + nbBinsre*j] << "\n");
	}
	}*/
    ofstream fout1("alpha.txt");
    for(CFuint j=0;j<totalNbCells;++j) {
      for(CFuint k=1;k<nbBinsre;++k) {
	fout1 << "alpha (" << k << "," << j << ") = " << 
	  alpha_avbin[k + nbBinsre*j] << "\n";
      }
    }
    fout1.close();
    
    ofstream fout2("beta.txt");
    for(CFuint j=0;j<totalNbCells;++j) {
      for(CFuint k=1;k<nbBinsre;++k) {
	fout2 << "beta (" << k << "," << j << ") = " << 
	  B_bin[k + nbBinsre*j] << "\n";
      }
    }
    fout2.close();
  }
  
  PE::GetPE().setBarrier(m_namespace);
}

//////////////////////////////////////////////////////////////////////////////

} // namespace RadiativeTransfer

} // namespace COOLFluiD

