#include <fstream>

#include "RadiativeTransfer/RadiationLibrary/Models/PARADE/ParadeRadiator.hh"

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
  options.addConfigOption< string >("LibraryPath","Path to Parade data files");
  options.addConfigOption< CFreal >("NDmin","Minimum number density.");
  options.addConfigOption< CFuint >("nbPoints","Number of Points to discretize the spectra per loop");
  options.addConfigOption< bool >                                 
    ("LTE", "True is it is a local thermodynamic equilibrium simulation");
  options.addConfigOption< vector<bool> >                                 
    ("MolecularSpecies", "Array of flags indicating whether each species is molecular.");
  options.addConfigOption< CFuint >("nbBands", "Number of Bands");
  options.addConfigOption< string >("bandsDistr", "Banding distribution");
  options.addConfigOption< CFuint >("nbBinsPARADE", "Number of Bins using PARADE");

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
  m_isLTE()
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
  
  m_molecularSpecies = vector<bool>();
  setParameter("MolecularSpecies", &m_molecularSpecies);

  m_nbBands = 1;
  setParameter("nbBands", &m_nbBands);

  m_nbBinsPARADE = 100;
  setParameter("nbBinsPARADE", &m_nbBinsPARADE);

  m_bandsDistr = "equally";
  setParameter("bandsDistr", &m_bandsDistr);
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
    
    m_avogadroOvMM.resize(nbSpecies);
    m_avogadroOvMM = PhysicalConsts::Avogadro()/m_mmasses;
   /* 
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
    */
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

CFuint m_nbBins = m_nbBinsPARADE; //for each band, to be set in the *.CFcase
CFuint nbCells = m_pstates->getSize();
RealVector vctBands(m_nbBands + 1); //vector containing the limits of bands

 if(m_nbBins == 1){
  CFLog(INFO, "ParadeLibrary::computeproprieties() => nbBins = 1 is not allowed, nbBins = 2 is set\n"); //otherwise dy = inf
  m_nbBins = 2;
 }

//variables for binning
CFuint nbBinsre = m_nbBins;
CFuint nbBandsre = m_nbBands;
m_alpha_bin.resize(nbBinsre*nbCells);
m_emission_bin.resize(nbBinsre*nbCells);
m_B_bin.resize(nbBinsre*nbCells*nbBandsre);       //bands included as bins
m_alpha_avbin.resize(nbBinsre*nbCells*nbBandsre); //bands included as bins
alpha.resize(nbCells);
epsilon.resize(nbCells);
B.resize(nbCells);
den_alpha_vol.resize(nbCells);
num_alpha_vol.resize(nbCells); 
RealVector vctBins(m_nbBins);
CFreal m_alpha_bincoeff;
CFreal m_emission_bincoeff;
CFreal m_B_bincoeff;

DataHandle<CFreal> alpha_avbin = m_radPhysicsHandlerPtr->getDataSockets()->alpha_avbin;
DataHandle<CFreal> B_bin = m_radPhysicsHandlerPtr->getDataSockets()->B_bin;

//logarithmic spacing for banding
 if(m_bandsDistr == "logarithmic") {
  CFLog(VERBOSE,"ParadeLibrary::computeproperties () => Logarithmic banding\n");
  CFreal wavMin_log = std::log(wavMin);
  CFreal wavMax_log = std::log(wavMax);
  CFreal dx = (wavMax_log-wavMin_log)/m_nbBands;
  vctBands[0] = wavMin;
  for(CFuint i = 1; i<m_nbBands; i++)
   {
    vctBands[i] = std::exp(wavMin_log + (dx * i));
   }
  vctBands[m_nbBands] = wavMax;
 }

 //equally-spaced distribution for banding
 else if(m_bandsDistr == "equally") {
  CFLog(VERBOSE,"ParadeLibrary::computeproperties () => Equally-spaced banding\n");
  CFreal dx = (wavMax-wavMin)/m_nbBands;
  vctBands[0] = wavMin;
  for(CFuint i = 1; i<m_nbBands; i++)
   {
    vctBands[i] = wavMin + (dx * i);
   }
  vctBands[m_nbBands] = wavMax;
 }

 for(CFuint i = 0; i<(m_nbBands + 1); i++) {
  CFLog(INFO,"ParadeLibrary::computeproperties () => vctBands(" << i << ") = " << vctBands[i] <<"\n");
 }

 //Number of points in each band for each cell
 std::vector< int > nb_bandPoints(m_nbBands, 0); //initialized to 1 to take into accout the starting point of each band
 
 for(CFuint ii=0; ii<m_nbPoints; ii++) {
  for(CFuint iii=0; iii<m_nbBands; iii++) {
    if(m_data(1, ii*3)>vctBands[iii] && m_data(1, ii*3)<vctBands[iii+1]) {
      //counting points within the band
      //on the 1st cell 'cause the overall points distribution does not change over cells
      nb_bandPoints[iii]++;
    }
   }
  }
 nb_bandPoints[0]++;
 nb_bandPoints[m_nbBands-1]++; //To include wavMax in the last band

 for(CFuint iii=0; iii<m_nbBands; iii++) {
   CFLog(INFO,"ParadeLibrary::computeproperties () => nb_bandPoints(" << iii  << ") = " << nb_bandPoints[iii] << "\n");
 }

 CFuint cc = 0; //counter to pick the proper proprieties from m_data
 for(CFuint b=0; b<m_nbBands; b++) { //the loop over bands starts

//Function to calculate the binning algorithm with the absorption and emission coefficients

 CFreal num_alphatot = 0.0;
 CFreal den_alphatot = 0.0;
 
 alphaav.clear();

 for(CFuint i=0;i<nb_bandPoints[b];i++)
 {

  for(CFuint s=0;s<nbCells;s++)
  {
   alpha[s] = m_data(s,(i+cc)*3+2);
   epsilon[s] = m_data(s,(i+cc)*3+1);
   B[s] = epsilon[s] / alpha[s];
     

// Planck function if we are not in equilibrium
//
// CFreal h = 6.626070040e-34;     //SI units Js
// CFuint c = 3e08; //SI units m/s
// CFreal k_b = 1.3806485279e-23;  //SI units J/K
//
//
// T needs to be defined from the values extracted
// const  B = (2*h*pow(c,2)/pow(m_data(j,i*3),5))*(1/(exp(h*c/(m_data(j,i*3)*k_b*T))-1));



 if(getCellVolume( m_pstates->getStateLocalID(s)) > 0)
 {
 num_alpha_vol[s] = alpha[s]*B[s]*getCellVolume( m_pstates->getStateLocalID(s));
 den_alpha_vol[s] = B[s]*getCellVolume( m_pstates->getStateLocalID(s));
 }
 else
 {
 CFLog(VERBOSE,"In binning function, negative or 0 value found" << "\n");
 CFLog(VERBOSE,"State Local ID = " << m_pstates->getStateLocalID(s) << "\n");
 num_alpha_vol[s] = alpha[s]*B[s]*-getCellVolume( m_pstates->getStateLocalID(s));
 den_alpha_vol[s] = B[s]*-getCellVolume( m_pstates->getStateLocalID(s));
 }

 }

 for(CFuint s = 0; s < nbCells; s++)
 {
   CFLog(VERBOSE,"Alpha = " << alpha[s] << "\n");
 }

 for(CFuint it = 0;it < nbCells; it++)
 {
    num_alphatot += num_alpha_vol[it];
    den_alphatot += den_alpha_vol[it];

 //CFLog(INFO,"ParadeLibrary::computeproperties () => num_alphatot" << num_alphatot <<"\n");
 
 }
 

 alphaav.push_back(num_alphatot/den_alphatot);


CFLog(VERBOSE,"ParadeLibrary::computeproperties () => alphaav = " << alphaav[i] <<"\n");

}


CFreal alphamin = alphaav[0];
CFreal alphamax = alphaav[0];

for(CFuint r=0;r<nb_bandPoints[b];r++)
{

 if(alphaav[r]<alphamin)
   {
 alphamin = alphaav[r];
   }
 if(alphaav[r]>alphamax)
   {
 alphamax = alphaav[r];
   }

}


CFLog(INFO,"ParadeLibrary::computeproperties () => alphamax = " << alphamax <<"\n");
CFLog(INFO,"ParadeLibrary::computeproperties () => alphamin = " << alphamin <<"\n");



// Logarithmic spacing for binning
//

 CFreal alpha_minlog = std::log(alphamin);
 CFreal alpha_maxlog = std::log(alphamax);

CFLog(VERBOSE,"ParadeLibrary::computeproperties () => alpha_minlog = " << alpha_minlog <<"\n");
CFLog(VERBOSE,"ParadeLibrary::computeproperties () => alpha_maxlog = " << alpha_maxlog <<"\n");

CFreal dy = (alpha_maxlog-alpha_minlog) / (m_nbBins-1);

CFLog(INFO,"ParadeLibrary::computeproperties () => dy = " << dy <<"\n");

 vctBins = 0.;
  for(int i = 0; i<m_nbBins; i++)
  {
  vctBins[i] = std::exp(alpha_minlog + (dy * i));
  CFLog(VERBOSE,"ParadeLibrary::computeproperties () => vctBins(" << i << ") = " << vctBins[i] <<"\n");
  }

  //To search for minimum alpha and maximum

CFreal alphamin_tot = alphaav[0];
CFreal alphamax_tot = alphaav[0];

for(CFuint i=0;i<nb_bandPoints[b];i++)
{
 for(CFuint s=0;s<nbCells;s++)
 {
   if((m_data(s,(i+cc)*3+2))<alphamin_tot)
   {
     alphamin_tot = m_data(s,(i+cc)*3+2);
   }
   if((m_data(s,(i+cc)*3+2))>alphamax_tot)
   {
     alphamax_tot = m_data(s,(i+cc)*3+2);
   }
 }
}

  vctBins[0] = alphamin_tot;
  vctBins[m_nbBins-1] = alphamax_tot;
    
  for( CFuint i = 0; i<m_nbBins;i++)
    {
  CFLog(INFO, " vctBins = " << vctBins[i] << "\n");
    }

// Binning part of the function
 
for(CFuint i=0;i<nb_bandPoints[b];i++)
{
 for(CFuint j=0;j<nbCells;j++){
    
   CFreal Bs = m_data(j,(i+cc)*3+1)/m_data(j,(i+cc)*3+2);

 m_alpha_bin[0+nbBinsre*j] = 0;
 m_emission_bin[0+nbBinsre*j] = 0;
 B_bin[0+nbBinsre*j] = 0;

   for(CFuint k=1;k<nbBinsre;k++){
        
     if((m_data(j,(i+cc)*3+2)>=vctBins[k-1]) && (m_data(j,(i+cc)*3+2)<vctBins[k]))
    {

      m_alpha_bincoeff = m_data(j,(i+cc)*3+2)*Bs*m_dWav*1.e-10;
      m_emission_bincoeff = m_data(j,(i+cc)*3+1)*m_dWav;
      m_B_bincoeff = Bs*m_dWav*1.e-10;
    /*
    CFLog(VERBOSE,"ParadeLibrary::computeproperties () => m_emission_bin(" << k << "," << j << ") = " << m_emission_bin[m_nbBins*j+k] <<"\n");
    CFLog(VERBOSE,"ParadeLibrary::computeproperties () => m_alpha_bin(" << k << "," << j << ") = " <<  m_alpha_bin[m_nbBins*j+k] <<"\n");
    CFLog(VERBOSE,"ParadeLibrary::computeproperties () => m_B_bin(" << k << "," << j << ") = " << m_B_bin[m_nbBins*j+k] <<"\n");
    */
    }

   
    m_alpha_bin[k + nbBinsre*j] += m_alpha_bincoeff;
    m_emission_bin[k + nbBinsre*j] += m_emission_bincoeff;
    m_B_bin[k + nbBinsre*j] += m_B_bincoeff;

  }
 }
}

for(CFuint k=1;k<nbBinsre;k++)
{
 for(CFuint j=0;j<nbCells;j++)
 {
    CFLog(VERBOSE,"ParadeLibrary::computeproperties () => m_alpha_bin(" << k << "," << j << ") = " <<  m_alpha_bin[nbBinsre*j+k] <<"\n");
    CFLog(VERBOSE,"ParadeLibrary::computeproperties () => m_B_bin(" << k << "," << j << ") = " << m_B_bin[nbBinsre*j+k] <<"\n");
  }
}

 for(int j=0;j<nbCells;j++)
   {
  for(int k=1;k<nbBinsre;k++)
  {
   m_alpha_avbin[0+nbBinsre*j] = 0.;
   
   if(m_B_bin[k + nbBinsre*j]!= 0.)
     {
     CFreal m_alpha_avbin_a = (m_alpha_bin[k + nbBinsre*j]) / (m_B_bin[k + nbBinsre*j]);
     m_alpha_avbin[k + nbBinsre*j] = m_alpha_avbin_a;
     CFLog(VERBOSE,"ParadeRadiator::computeProperties()=>alpha_avbin(" << k << "," << j << ") = "<<m_alpha_avbin[k + nbBinsre*j]<<"\n");
     }
   else
     {
     alpha_avbin[k + nbBinsre*j] = 0.;
     CFLog(VERBOSE,"ParadeRadiator::computeProperties()=>alpha_avbin(" << k << ","<< j << ") = "<<alpha_avbin[k + nbBinsre*j]<<"\n");
     }
 }
}

//alpha_avbin, emission_bin are the average values for absorptivity and emissivity for each bin.

for(int j=0;j<nbCells;j++)
 {
   for(int k=1;k<nbBinsre;k++)
  {
   //the index "b" added for banding
    alpha_avbin[(k+b) + nbBinsre*j] = m_alpha_avbin[k + nbBinsre*j];
    B_bin[(k+b) + nbBinsre*j] = m_B_bin[k + nbBinsre*j];
    CFLog(VERBOSE,"alpha (" << b << "," << k << "," << j << ") = " << alpha_avbin[(k+b) + nbBinsre*j] << "\n");
  }
 }
  cc = cc + nb_bandPoints[b];

  PE::GetPE().setBarrier(m_namespace);

 }
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
  
  CFLog(INFO, "ParadeRadiator::writeLocalData() => written densities for cells [" 
	<< startNode << ", " << countNode << "]\n");
  
  CFLog(VERBOSE, "ParadeRadiator::writeLocalData() => START\n");
}   
      
//////////////////////////////////////////////////////////////////////////////
  
void ParadeRadiator::readLocalRadCoeff()
{
  CFLog(VERBOSE, "ParadeRadiator::readLocalRadCoeff() => START\n");
  
  fstream& fin = m_inFileHandle->openBinary(m_radFile);
  //string nam = "file-" + StringOps::to_str(PE::GetPE().GetRank());
  
  int one = 0;
  fin.read((char*)&one, sizeof(int));
  //fout << "one = "<< one;
  cf_assert(one == 1);
  
  int nbCells = 0;
  fin.read((char*)&nbCells, sizeof(int));
  CFLog(INFO,"nbCells="<<nbCells<<"\n");
    //fout << " nbCells = " << nbCells << endl;
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
  m_data.resize(totalNbCells, m_nbPoints*3, 0.);
  
  // here if the process stores the full mesh, each processor reads 
  // only a portion of the data
  RealMatrix partialData;
  RealMatrix* currData = &m_data;
  if (fullGridInProcess()) {
    partialData.resize(nbCells, m_nbPoints*3, 0.);
    currData = &partialData;
  }
  
  double etot = 0.;
  int wavpts = 0;
  const CFuint sizeCoeff = m_nbPoints*3;
  for (int iPoint = 0 ; iPoint < nbCells; ++iPoint) {
    fin.read((char*)&etot, sizeof(double));
    //fout << "etot = " << etot << endl;
    fin.read((char*)&wavpts, sizeof(int));
    cf_assert(wavpts == (int)m_nbPoints);
    //fout << "wavpt = " <<  wavpts << endl;
    // this reads [wavelength, emission, absorption] for each cell
    fin.read((char*)&((*currData)[iPoint*sizeCoeff]), sizeCoeff*sizeof(double));
  }
  fin.close();
  
  if (fullGridInProcess()) {
    // in case the full mesh is stored in each process, since we have read in only a part 
    // of spectral data, we need now to gather all data so that each process has the full 
    // spectra before performing binning 
    // AL: this is just an intermediate step towards a truly parallel binning calculation
    vector<int> recvCounts(m_nbProc, 0);
    vector<int> displs(m_nbProc, 0);
    const CFuint nbCellsPerProc = totalNbCells/m_nbProc;
    const CFuint minSizeToSend = nbCellsPerProc*sizeCoeff;
    const CFuint maxSizeToSend = minSizeToSend + (totalNbCells%m_nbProc)*sizeCoeff;
    CFuint count = 0;
    for (CFuint rank = 0; rank < m_nbProc; ++rank) {
      recvCounts[rank] = (rank < m_nbProc-1) ? minSizeToSend : maxSizeToSend;
      // displacements start at 0
      displs[rank] = count;
      count += minSizeToSend;
    }
    
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

  /*
if(m_writePARADEToFile)
{
CFLog(INFO,"Writing absorption and emission coefficients to a file" << "\n");

boost::filesystem::path file = m_dirName / boost::filesystem::path(m_outTabName);
file = PathAppender::getInstance().appendParallel(file);

 SelfRegistPtr<Environment::FileHandlerOutput>fhandle = 

    Environment::SingleBehaviorFactory<Environment::FileHandlerOutput>::getInstance().create();


ofstream& fout = fhandle->open(file);

//Writing m_data into a file.

for(CFuint p=0;p<m_nbPoints;++p)
{
 for(CFuint m=0;m<nbCells;++m)
 {
 fout << m_data(m,p*3+2) << " ";
 fout << m_data(m,p*3+1) << " ";
 }
}

fhandle->close();

}
*/
    /*    
  boost::filesystem::path file = m_dirName / boost::filesystem::path("Alfa_eps_cell1500.plt");
  file = PathAppender::getInstance().appendParallel(file);

  SelfRegistPtr<Environment::FileHandlerOutput>fhandle = 

    Environment::SingleBehaviorFactory<Environment::FileHandlerOutput>::getInstance().create();

  ofstream& outputFile = fhandle->open(file);

  
  for(CFuint i=0;i<m_nbPoints;++i)
    {
  outputFile << m_data(1500,i*3) << " " << m_data(1500,i*3+1) << " " << m_data(1500,i*3+2) << endl;
    }
  
  fhandle->close();

    */

  //fout.close();
  
  // if (PE::GetPE().GetRank() == 0) {
  //     ofstream fout1("inwav.txt"); 
//     fout1 << "TITLE = Original spectrum of radiative properties\n";
//     fout1 << "VARIABLES = Wavl EmCoef AbCoef \n";
//     for (CFuint i = 0; i < data.nbCols()/3; ++i) {
//       fout1 << data(0,i*3) << " " << data(0,i*3+1) << " " << data(0,i*3+2) << endl;
//     }
//     fout1.close();
//   }
  
//   // Convert spectrum to frequency space
//   const CFreal c = PhysicalConsts::LightSpeed();
//   const CFuint stride = data.nbCols()/3;
//   for (CFuint iPoint = 0 ; iPoint < nbCells; ++iPoint) {
//     for (CFuint iw = 0; iw < stride; ++iw) {
//       const CFreal lambda = data(iPoint, iw*3)*1.e-10; // conversion from [A] to [m]
//       // convert emission from wavelength spectrum to frequency spectrum
//       data(iPoint, iw*3+1) *= c/(lambda*lambda);
//       // override wavelength with frequency (nu = c/lambda)
//       data(iPoint, iw*3) = c/lambda; 
//     }
//   }
  
//   if (PE::GetPE().GetRank() == 0) {
//     ofstream fout1("infreq.txt"); 
//     fout1 << "TITLE = Original spectrum of radiative properties\n";
//     fout1 << "VARIABLES = Freq EmCoef AbCoef \n";
//     for (CFuint i = 0; i < data.nbCols()/3; ++i) {
//       fout1 << data(0,i*3) << " " << data(0,i*3+1) << " " << data(0,i*3+2) << endl;
//     }
//     fout1.close();
//   }
  
  CFLog(VERBOSE, "ParadeRadiator::readLocalRadCoeff() => END\n");
}

//////////////////////////////////////////////////////////////////////////////
  
inline void ParadeRadiator::getSpectralIdxs(CFreal lambda, CFuint& idx1, CFuint& idx2)
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
  
/// array storing absorption and emission coefficients
/// wavelength       = m_radCoeff(local state ID, spectral point idx*3)
/// emission coeff   = m_radCoeff(local state ID, spectral point idx*3+1)
/// absorption coeff = m_radCoeff(local state ID, spectral point idx*3+2)
///

//////////////////////////////////////////////////////////////////////////////
  
CFreal ParadeRadiator::getEmission(CFreal lambda, RealVector &s_o)
{
  CFuint spectralIdx1 = 0;
  CFuint spectralIdx2 = 0;
  const CFuint stateIdx = m_radPhysicsHandlerPtr->getCurrentCellTrsIdx();
  getSpectralIdxs(lambda, spectralIdx1, spectralIdx2);
 
  cf_assert(stateIdx < m_data.nbRows());
  cf_assert(spectralIdx1*3   < m_data.nbCols());
  cf_assert(spectralIdx1*3+1 < m_data.nbCols());
  cf_assert(spectralIdx2*3   < m_data.nbCols());
  cf_assert(spectralIdx2*3+1 < m_data.nbCols());

  const CFreal x0 = m_data( stateIdx, spectralIdx1*3   );
  const CFreal y0 = m_data( stateIdx, spectralIdx1*3+1 );

  const CFreal x1 = m_data( stateIdx, spectralIdx2*3 );
  const CFreal y1 = m_data( stateIdx, spectralIdx2*3+1 );

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
  
  cf_assert(stateIdx < m_data.nbRows());
  cf_assert(spectralIdx1*3   < m_data.nbCols()); 
  cf_assert(spectralIdx1*3+2 < m_data.nbCols());
  cf_assert(spectralIdx2*3   < m_data.nbCols());
  cf_assert(spectralIdx2*3+2 < m_data.nbCols());
  
  const CFreal x0 = m_data( stateIdx, spectralIdx1*3   );
  const CFreal y0 = m_data( stateIdx, spectralIdx1*3+2 );

  const CFreal x1 = m_data( stateIdx, spectralIdx2*3 );
  const CFreal y1 = m_data( stateIdx, spectralIdx2*3+2 );

  //linear interpolation
  return y0 + (y1-y0) * (lambda - x0) / (x1-x0);
}

//////////////////////////////////////////////////////////////////////////////

void ParadeRadiator::computeEmissionCPD()
{
  CFuint nbCells = m_pstates->getSize();
  CFuint nbCpdPoints = m_nbPoints;

  m_spectralLoopPowers.resize( nbCells );
  m_cpdEms.resize( nbCells * nbCpdPoints );

  CFreal emInt;
  
  for (CFuint s = 0; s < nbCells; ++s) {
    //first pass to get the total emission coefficient
    //use this information to get the cell's total Radiative Power

    //cout<<"firrst pass, cell "<<s<<' '<<m_dWav<<endl;
    emInt=0.;
    m_cpdEms[ s*nbCpdPoints + 0 ] = 0;
    for (CFuint j=1; j < nbCpdPoints; j++ ){
      //cout<<"spectral point "<<j<<" of "<<nbCpdPoints<<endl;
      emInt += ( m_data(s,(j-1)*3+1) + m_data(s,(j)*3+1) )/2. * m_dWav;
      m_cpdEms[ s*nbCpdPoints + j ] = emInt;
    }
    //cout<<"done, cell "<<s<<' '<<emInt <<endl;
    m_spectralLoopPowers[s]=emInt*getCellVolume( m_pstates->getStateLocalID(s)) * m_angstrom * 4.0 * 3.1415926535897;
    if(getCellVolume( m_pstates->getStateLocalID(s)) <= 0)
    {
    CFLog(VERBOSE,"In compute Emission CPD, negative or 0 value found" << "\n");
    CFLog(VERBOSE,"State Local ID = " << m_pstates->getStateLocalID(s) << "\n");
    }


    //cout<<"second pass, cell "<<s<<endl;
    //second pass to get the [0->1] comulative probably distribution
    for (CFuint j=1; j < nbCpdPoints; j++ ){
      m_cpdEms[ s*nbCpdPoints + j ] /= emInt;
    }
  }

/*  if (  PE::GetPE().GetRank()  == 0) {
  //test for the first cell
  cout<<endl<<" wav = [";
  for (CFuint j=0; j < m_nbPoints; j++ ){
    cout<<m_data(0,j*3+0)<<' ';
  }
  cout<<" ];" <<endl<<" em = [";
  for (CFuint j=0; j < m_nbPoints; j++ ){
    cout<<m_data(0,j*3+1)<<' ';
  }
  cout<<" ];"<<endl<<" am = [";
  for (CFuint j=0; j < m_nbPoints; j++ ){
    cout<<m_data(0,j*3+2)<<' ';
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

  const CFreal x0 = m_cpdEms[ stateIdx*nbCpdPoints + spectralIdx1 ];
  const CFreal y0 = m_data( stateIdx, spectralIdx1*3 );

  const CFreal x1 = m_cpdEms[ stateIdx*nbCpdPoints + spectralIdx2 ];
  const CFreal y1 = m_data( stateIdx, spectralIdx2*3 );

  lambda = y0 + (y1-y0) * (rand - x0) / (x1-x0);

  //cout<<x0<<' '<<x1<<' '<<y0<<' '<<y1<<' '<<rand<<' '<<lambda<<endl;

  m_rand.sphereDirections(dim2, s_o);
}

//////////////////////////////////////////////////////////////////////////////

void ParadeRadiator::getData()
{
}

//////////////////////////////////////////////////////////////////////////////

} // namespace RadiativeTransfer

} // namespace COOLFluiD

