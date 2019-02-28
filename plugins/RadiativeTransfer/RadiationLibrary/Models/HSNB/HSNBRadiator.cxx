#include "RadiativeTransfer/RadiationLibrary/Models/HSNB/HSNBRadiator.hh"

#include <fstream>

#include "RadiativeTransfer/RadiationLibrary/RadiationPhysicsHandler.hh"
#include "RadiativeTransfer/RadiationLibrary/RadiationPhysics.hh"

#include "Common/CFLog.hh"
#include "Common/DebugFunctions.hh"
#include "Common/CFPrintContainer.hh"
#include "Common/OSystem.hh"

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

using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Environment;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

namespace RadiativeTransfer {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<HSNBRadiator,
                Radiator,
                RadiativeTransferModule,
                1>
HSNBRadiatorProvider("HSNBRadiator");

//////////////////////////////////////////////////////////////////////////////

void HSNBRadiator::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< string >("LibraryPath","Path to HSNB data files.");
  options.addConfigOption< std::string >("LocalDirName","Name of the local temporary directories where HSNB is run.");
  options.addConfigOption< CFuint >("TvIndex", "Index of Tv in the state vector");
  options.addConfigOption< CFuint >("TrIndex", "Index of Tr in the state vector");
  options.addConfigOption< CFuint >("pIndex", "Index of the total pressure in the state vector");
  options.addConfigOption< bool, Config::DynamicOption<> >
    ("ReuseProperties", "Reuse existing radiative data (requires the same number of processors as in the previous run).");
  options.addConfigOption< std::string >("NonUniformPathTreatment","Specifies how to treat radiative properties along a non-uniform path");
  options.addConfigOption< std::string >("Namespace","Namespace that the HSNBRadiator is run in in parallel");
  options.addConfigOption< std::string >("compositionType", "Describe how species composition is described (moleFractions/partialDensity)");
  options.addConfigOption< CFreal >("WavenumberMin", "Specify a global lower bound for the spectra which significantly influences memory reqs.");
  options.addConfigOption< CFreal >("WavenumberMax", "Specify a global lower bound for the spectra which significantly influences memory reqs.");
  options.addConfigOption< bool, Config::DynamicOption<> >
    ("UsePrecomputedProperties", "Precomputes parameters for emission and absorption for all cells and bands in all species EXCEPT continua. Potentially takes up a lot of memory.");
  options.addConfigOption< bool, Config::DynamicOption<> >
    ("UsePrecomputedContinuumProperties", "Precomputes parameters for emission and absorption for all cells and bands for continua systems. Potentially takes up a lot of memory but saves computation time since evaluation of continua is generally most expensive");


}

void HSNBRadiator::configure ( Config::ConfigArgs& args )
{
  Radiator::configure(args);
  ConfigObject::configure(args);
}

void HSNBRadiator::unsetup()
{
  Radiator::unsetup();
}

void HSNBRadiator::setupSpectra(CFreal wavMin, CFreal wavMax)
{

    //NOTE Maybe care for thread safety, also move outside public function

//    determineBandRange();

//    // Compute radiative properties in each cell for each system
//    std::cout << "Computing spectral features... ";
//    for (int i = 0; i < m_diatomics.size(); ++i)
//        m_diatomics[i].setupLocalParameters(m_thermoData);
//    for (int i = 0; i < m_continua.size(); ++i)
//        m_continua[i].setupLocalParameters(m_thermoData);
//    for (int i = 0; i < m_atoms.size(); ++i)
//        m_atoms[i]->setupLocalParameters(m_thermoData);
//    std::cout << "Done!" << std::endl;

//    // Prepare the line emission array
//    int nlines = 0;
//    for (int i = 0; i < m_atoms.size(); ++i)
//        nlines += m_atoms[i]->nLines();
//    if (nlines > 0) {
//        m_line_emis = new RealMatrix(m_atoms.size(),nlines);
//    }

//    //Load black body cumulative energy
//    load_BBCE_data();

}

CFreal HSNBRadiator::getEmission(CFreal lambda, RealVector &s_o)
{

}



void HSNBRadiator::computeEmissionCPD()
{

}

void HSNBRadiator::getRandomEmission(CFreal &lambda, RealVector &s_o)
{
    getRandomDirection(s_o);
}

void HSNBRadiator::getRandomDirection(RealVector &s_o)
{
    static CFuint dim = Framework::PhysicalModelStack::getActive()->getDim();
    static CFuint dim2 = m_radPhysicsHandlerPtr->isAxi() ? 3 : dim;
    m_rand.sphereDirections(dim2, s_o);
}

void HSNBRadiator::getRandomBBSigma(CFreal &sig, CFreal temperature)
{
    double f=Random::uniform(0.00032, 0.9989);
    for (int i=0; i<m_BBCE.cumulative.size(); i++){
        if (f<m_BBCE.cumulative[i]) {
            double deltaf=m_BBCE.cumulative[i]-f;
            double deltafTot=m_BBCE.cumulative[i]-m_BBCE.cumulative[i-1];
            double deltaLamdaTTot=m_BBCE.TTimesLamda[i]-m_BBCE.TTimesLamda[i-1];
            double lamdaT=m_BBCE.TTimesLamda[i]-(deltaf/deltafTot)*deltaLamdaTTot;
            //check if is within the range
            cf_assert(lamdaT<=50000);
            cf_assert(lamdaT>=1000);

            sig=temperature*10000/lamdaT;
            return;
            }
    }
    std::cout << "ERROR: OUT    ";
    sig=0;
}

void HSNBRadiator::addNewCellAbsorption(CFreal &ku, CFreal &param2, CFreal &param3, CFreal& wavenumber, CFreal& distance, CFuint& mechanismID)
{
//    const CFuint cellID = m_radPhysicsHandlerPtr->getCurrentCellStateID();


//    if (mechanismIsDiatomic(mechanismID)) {
//        m_diatomics[mechanismID].addNewCellAbsorption(CFreal &ku, CFreal &param2, CFreal &param3, CFreal& wavenumber, CFreal& distance, CFuint& mechanismID);
//    }
//    else if (mechanismIsContinous(mechanismID)) {
//        CFuint localMechanism=mechanismID-m_nbDiatomics;
//        m_continua[localMechanism].addNewCellAbsorption(CFreal &ku, CFreal &param2, CFreal &param3, CFreal& wavenumber, CFreal& distance, CFuint& mechanismID);
//    }
//    else {
//        m_atoms
//    }

//    int b = int(wavenumber / 1000.0);

//    // Convert band to local indexing
//    if (b < lowBand() || b > highBand())
//        return 0.0;
//    b -= lowBand();

//    switch (m_nup) {

//    case THIN:

//        ku = mp_locparams[cellID*2*m_nbands + b*2]*distance;

//    case CURTIS_GODSON:

//        ku = mp_locparams[cellID*4*m_nbands+b*4+0]*distance;
//        param2=mp_locparams[cellID*4*m_nbands+b*4+2];
//        param3=mp_locparams[cellID*4*m_nbands+b*4+3];

//    case LINDQUIST_SIMMONS:

//        ku = mp_locparams[icell*4*m_nbands+b*4+0]*distance;
//        param2=mp_locparams[cellID*4*m_nbands+b*4+2];
//        param3=mp_locparams[cellID*4*m_nbands+b*4+3];

//    }




}


CFreal HSNBRadiator::getTotalEmissivePower() const
{
    return m_totalPower;
}

CFreal HSNBRadiator::getCurrentCellEmissivePower() const
{
    return m_totalEmissionCurrentCell;
}

CFuint HSNBRadiator::getNbNonThickDiatomics() const
{
    return m_nbNonThickDiatomics;
}

CFuint HSNBRadiator::getNbThickDiatomics() const
{
    return m_nbThickDiatomics;
}

CFuint HSNBRadiator::getNbAtoms() const
{
    return m_atoms.size();
}

CFuint HSNBRadiator::getNbContinua() const
{
    return m_continua.size();
}

bool HSNBRadiator::readEmissionData(LagrangianSolver::ParallelVector<CFreal> &stateRadPower, std::vector<CFreal> &gStateRadPower, CFreal &totalStatePower, CFreal &totalGhostStateRadPower)
{
    if (m_reuseProperties) {
        boost::filesystem::path loc("precomputedEmissionData");
        boost::filesystem::path infile = m_hsnbDir / loc;
        bool readSuccessful=false;

        m_inFileHandle  = COOLFluiD::Environment::SingleBehaviorFactory<COOLFluiD::Environment::FileHandlerInput>::getInstance().create();

        // Open the file
        std::ifstream& fin = m_inFileHandle->open(infile);



        if (fin.is_open()) {

            std::string line;

            //Load total state power
            std::getline(fin, line);
            totalStatePower=atof(line.c_str());

            //Load total ghost-state power
            std::getline(fin, line);
            totalGhostStateRadPower=atof(line.c_str());

            for (int i=0; i<stateRadPower.size(); i++) {
                std::getline(fin, line);
                stateRadPower[i]=atof(line.c_str());
            }

            for (int i=0; i<gStateRadPower.size(); i++) {
                std::getline(fin, line);
                gStateRadPower[i]=atof(line.c_str());
            }

            readSuccessful=true;
        }

        else {
            CFLog(ERROR, "HSNBRadiator::readEmissionData => File '"<< infile.c_str() << " couldn't be opened. " << "\n");
        }

        m_inFileHandle->close();
        return readSuccessful;
    }
    else {
        return false;
    }
}

void HSNBRadiator::exportEmissionData(LagrangianSolver::ParallelVector<CFreal> &stateRadPower, std::vector<CFreal> &gStateRadPower, CFreal totalStatePower, CFreal totalGhostStateRadPower)
{
    if (m_reuseProperties==false) {
        CFLog(INFO, "HSNBRadiator::exportEmissionData =>  = Exporting emission data... \n");
        COOLFluiD::Common::SelfRegistPtr<COOLFluiD::Environment::FileHandlerOutput> outFileHandle  = COOLFluiD::Environment::SingleBehaviorFactory<COOLFluiD::Environment::FileHandlerOutput>::getInstance().create();

        boost::filesystem::path loc("precomputedEmissionData");
        boost::filesystem::path outfile = m_hsnbDir / loc;
        ofstream& fout = outFileHandle->open(outfile);

        fout << totalStatePower << std::endl;
        fout << totalGhostStateRadPower << std::endl;

        for (int i=0; i<stateRadPower.size(); i++) {
            fout << stateRadPower[i] << std::endl;
        }

        for (int i=0; i<gStateRadPower.size(); i++) {
            fout << gStateRadPower[i] << std::endl;
        }

        fout.close();
    }
}

void HSNBRadiator::setUpRadiationField()
{
//    Framework::DataSocketSink < Framework::State* , Framework::GLOBAL > states = m_radPhysicsHandlerPtr->getDataSockets()->states;
//    ///Handle to states associated with this radiator
//    Framework::DataHandle<Framework::State*, Framework::GLOBAL>& statesHandle=states.getDataHandle();
//    RealVector* stateData;

    Framework::DataSocketSink < Framework::State* , Framework::GLOBAL > states =
    m_radPhysicsHandlerPtr->getDataSockets()->states;
    DataHandle<State*, GLOBAL> statesHandle = states.getDataHandle();

    RealVector* stateData;
    CFLog(INFO, "HSNBRadiator::setUpRadiationField => START\n");
    this->m_thermoData.setup(m_pI,m_trI,m_tvI,&m_avogadroOvMM, m_convertPartialDensity);

    CFLog(INFO, "HSNBRadiator::setUpRadiationField => Setup " << m_statesID.size() << " states. \n");

    std::cout << "HSNBRadiator::setUpRadiationField => P"<< m_rank <<" Setup " << m_statesID.size() << " states. \n";
    CFLog(INFO, "HSNBRadiator::setup() => P " << m_rank << " nbSpecies="<<  m_thermoData.nSpecies() << "\n");


    for (int i=0; i<m_statesID.size(); i++) {
        //Get current state
//        std::cout << "HSNBRadiator::setUpRadiationField => i=" << i << ", m_statesID[" << i<< "]="
//                  << m_statesID[i] << std::endl;
        stateData= (statesHandle[m_statesID[i]]->getData());

        m_thermoData.addState(stateData, m_statesID[i]);
    }
//    std::cout << "HSNBRadiator::setUpRadiationField => P"<< m_rank<<" HOLDING " << m_thermoData.nCells() << " CELLS \n";


    CFLog(INFO, "HSNBRadiator::setUpRadiationField => HOLDING " << m_thermoData.nCells() << " CELLS \n");
    CFLog(INFO, "HSNBRadiator::setUpRadiationField => P0 FINISHED SETUP\n");


}

void HSNBRadiator::setUpMechanisms()
{
    double diatomicPreomputedSizeMB=0.0;

    determineBandRange();


//    cf_assert(m_usePrecomputedParameters==false);

    if (m_usePrecomputedParameters) {
        m_thermoData.precomputeParameters(true);
    }
    else {
        m_thermoData.precomputeParameters(false);
//        m_reuseProperties=false;
    }


    if (m_usePrecomputedContinuumParameters) {
        m_thermoData.precomputeContinuumParameters(true);
    }
    else {
        m_thermoData.precomputeContinuumParameters(false);
        m_reuseProperties=false;
    }


    Stopwatch<WallTime> stp;
    Stopwatch<WallTime> stpMechanism;
    stp.start();
    stpMechanism.start();
    // Compute radiative properties in each cell for each system   
    CFLog(INFO, "HSNBRadiator::setUpMechanisms => Computing spectral features... \n");



    if (m_co2Exists) {

        std::vector<CFreal> volumeVec;
        for (int i=0; i<m_pstates->getSize();i++) {
            volumeVec.push_back(getCellVolume( m_pstates->getStateLocalID(i)));
        }

        m_co2->setupProfileType(volumeVec, m_thermoData);
        m_co2->setupLocalParameters(m_thermoData);
        CFLog(INFO, "HSNBRadiator::setUpMechanisms => CO2 System loaded in " << stpMechanism.read() << " \n");
    }


    for (int i = 0; i < m_diatomics.size(); ++i) {
        m_diatomics[i].setupLocalParameters(m_thermoData);
        diatomicPreomputedSizeMB+=m_diatomics[i].getPrecomputedBufferByteSize()*(1e-6);

        CFLog(DEBUG_MAX, "HSNBRadiator::setUpMechanisms => Precomputed parameters in diatomic system " << m_diatomics[i].systemName() << " use " << m_diatomics[i].getPrecomputedBufferByteSize()*(1e-6) << "MB of memory \n");
    }
    stpMechanism.stop();
    CFLog(INFO, "HSNBRadiator::setUpMechanisms => Diatomic features loaded in " << stpMechanism.read() << " \n");


    std::cout << "HSNBRadiator::setUpMechanisms => P"<< m_rank << " Total precomputed size of diatomics is " << diatomicPreomputedSizeMB << " MB." << std::endl;

 //    CFLog(INFO, "HSNBRadiator::setUpMechanisms => m_continua.size()=" << m_continua.size()<< " \n");
    stpMechanism.restart();
    if (m_reuseProperties==false) {
        for (int i = 0; i < m_continua.size(); ++i) {
            m_continua[i].setupLocalParameters(m_thermoData);
            m_continua[i].exportLocalParameters(m_hsnbDir, m_thermoData);
            CFLog(INFO, "HSNBRadiator::setUpMechanisms => Load" << m_continua[i].speciesName() << " / " << m_continua[i].systemName() <<"\n");
            CFLog(INFO, "HSNBRadiator::setUpMechanisms => Setting up local rad data... " << std::setprecision(3) << (i/CFreal(m_continua.size()))*100.0 << "% \n");
        }
    }
    else {
        for (int i = 0; i < m_continua.size(); ++i) {
            m_continua[i].readLocalParameters(m_hsnbDir);
            CFLog(INFO, "HSNBRadiator::setUpMechanisms => Load" << m_continua[i].speciesName() << " / " << m_continua[i].systemName() <<"\n");
            CFLog(INFO, "HSNBRadiator::setUpMechanisms => Loading local rad data... " << std::setprecision(3) << (i/CFreal(m_continua.size()))*100.0 << "% \n");
        }
    }
    stpMechanism.stop();
    CFLog(INFO, "HSNBRadiator::setUpMechanisms => Continua features loaded in " << stpMechanism.read() <<  " \n");
    stpMechanism.restart();
    for (int i = 0; i < m_atoms.size(); ++i) {
        m_atoms[i].setupLocalParameters(m_thermoData);
    }
    stpMechanism.stop();
    CFLog(INFO, "HSNBRadiator::setUpMechanisms => Atomic features loaded in " << stpMechanism.read() <<  " \n");

    // Prepare the line emission array
    int nlines = 0;
    for (int i = 0; i < m_atoms.size(); ++i)
        nlines += m_atoms[i].nLines();
    if (nlines > 0) {
        m_line_emis=std::vector<CFreal>(nlines,0.0);
    }

    std::cout << "m_line_emis.size=" << m_line_emis.size() << std::endl;

    stp.stop();

    CFLog(INFO, "HSNBRadiator::setUpMechanisms => Spectral features loaded in "<<stp.read()<< "\n");

}

void HSNBRadiator::load_BBCE_data()
{

    // Open the file
    std::ifstream& BBCE_file = m_inFileHandle->open(m_blackBodyFile);

    if (BBCE_file.is_open()) {

    std::string line;
    //skip first two line
    std::getline(BBCE_file, line);
    std::getline(BBCE_file, line);
    //load all the data
    std::getline(BBCE_file, line);

    while (!line.empty()) {
        double temp_data;
        std::stringstream ss(line);
        ss >> temp_data;
        m_BBCE.TTimesLamda.push_back(temp_data);
        ss >> temp_data;
        m_BBCE.cumulative.push_back(temp_data);
        std::getline(BBCE_file, line);
    }

    }

    else {std::cout << "file BBCE not open";}

    m_inFileHandle->close();
}

void HSNBRadiator::determineBandRange()
{
    //Arbitrary initializer
    m_bmin=1;
    m_bmax=1;

    if (m_co2Exists) {
        m_bmin = m_co2->lowBand();
        m_bmax = m_co2->highBand();
    }
    else if (m_diatomics.size() != 0) {
        m_bmin = m_diatomics[0].lowBand();
        m_bmax = m_diatomics[0].highBand();
    }
    else if (m_continua.size() != 0)
    {
        m_bmin = m_continua[0].lowBand();
        m_bmax = m_continua[0].highBand();
    }

    for (int i = 0; i < m_diatomics.size(); ++i) {
        m_bmin = std::min(m_bmin, m_diatomics[i].lowBand());
        m_bmax = std::max(m_bmax, m_diatomics[i].highBand());
    }

    for (int i = 0; i < m_continua.size(); ++i) {
        m_bmin = std::min(m_bmin, m_continua[i].lowBand());
        m_bmax = std::max(m_bmax, m_continua[i].highBand());
    }

    m_nbands = m_bmax-m_bmin+1;
}

HSNBRadiator::HSNBRadiator(const std::string& name) :
    Radiator(name),
    m_statesID(),
    m_rank(0),
    m_nbProc(1),
    m_inFileHandle()
  {
    addConfigOptionsTo(this);

    m_libPath = "";
    setParameter("LibraryPath", &m_libPath);
    
    m_localDirName = "HSNB";
    setParameter("LocalDirName", &m_localDirName);

    m_tvI=0;
    setParameter("TvIndex", &m_tvI);

    m_trI=0;
    setParameter("TrIndex", &m_trI);
    
    m_pI=0;
    setParameter("pIndex", &m_pI);
    
    m_namespace = "Default";
    setParameter("Namespace", &m_namespace);

    m_reuseProperties = false;
    setParameter("ReuseProperties", &m_reuseProperties);

    m_usePrecomputedParameters =true;
    setParameter("UsePrecomputedProperties", &m_usePrecomputedParameters);

    m_usePrecomputedContinuumParameters = true;
    setParameter("UsePrecomputedContinuumProperties", &m_usePrecomputedContinuumParameters);


    m_nuPathTreatment = "UNIFORM";
    setParameter("NonUniformPathTreatment", &m_nuPathTreatment);

    m_compositionType= "partialDensity";
    setParameter("compositionType", &m_compositionType);

    m_sigMin = -1.0; // default value is on purpose <0
    setParameter("WavenumberMin", &m_sigMin);

    m_sigMax = -1.0; // default value is on purpose <0
    setParameter("WavenumberMax", &m_sigMax);

    if (m_compositionType=="partialDensity") {
        m_convertPartialDensity=true;
    }
    else {
        m_convertPartialDensity=false;
    }

}

HSNBRadiator::~HSNBRadiator()
{

}


void HSNBRadiator::setup(){
    CFLog(INFO, "HSNBRadiator::setup() => START\n");

    Radiator::setup();

    // Be sure to seed the random number generator


    // create a m_localDirName-P# directory for the current process
    const string hsnbDir = m_localDirName + "_" + m_radPhysicsPtr->getTRSname() +
      "-P" + StringOps::to_str(PE::GetPE().GetRank(m_namespace));
    boost::filesystem::path hsnbPath(hsnbDir);

    m_hsnbDir = Environment::DirPaths::getInstance().getWorkingDir() / hsnbPath;
    CFLog(VERBOSE, "HSNBRadiator::setup() => m_hsnbDir = " << m_hsnbDir << "\n");


    if (m_reuseProperties==false) {
        //Create inout directory (and delete old data)
        std::string command1   = "rm -fr " + m_hsnbDir.string() + " ; mkdir " + m_hsnbDir.string();
        Common::OSystem::getInstance().executeCommand(command1);
    }


    m_rank   = PE::GetPE().GetRank(m_namespace);
    m_nbProc = PE::GetPE().GetProcessorCount(m_namespace);

    //Seed the random number generator
    Random::seed(time(NULL)*(m_rank+1));

    m_inFileHandle  = Environment::SingleBehaviorFactory<Environment::FileHandlerInput>::getInstance().create();

    if (m_libPath == "") {
        CFLog(ERROR, "ParadeLibrary::setLibrarySequentially() => libpath NOT SET\n");
        boost::filesystem::path HSNBPath("/plugins/RadiativeTransfer/RadiationLibrary/Models/HSNB");
        m_HSNBPath = Environment::DirPaths::getInstance().getBaseDir() / HSNBPath;
    }
    else {
        m_HSNBPath = m_libPath;
    }


    CFLog(VERBOSE, "HSNBRadiator::setup() => m_HSNBPath = " << m_HSNBPath << "\n");


    MPI_Request sendRequest, recvRequest;

    int active=0;
    if (m_rank == 0) {
        active = 1;
        //Copy datafiles into proc. directory
        std::string command2   = "cp -r " + m_HSNBPath.string() + "/data" + " " + m_hsnbDir.string();
        Common::OSystem::getInstance().executeCommand(command2);
        command2 = "cp -r " + Environment::DirPaths::getInstance().getWorkingDir().string() + "/hsnb.con" + " " + m_hsnbDir.string();
        Common::OSystem::getInstance().executeCommand(command2);

        m_HSNBPath = m_hsnbDir;

        CFLog(VERBOSE, "HSNBRadiator::setup() => Copy datafiles into = " << m_hsnbDir.string() << "\n");

        MPI_Isend(&active, 1, MPIStructDef::getMPIType(&active), m_rank+1, 0, PE::GetPE().GetCommunicator(m_namespace), &sendRequest);
    }
    else {

        MPI_Recv(&active, 1, MPIStructDef::getMPIType(&active), m_rank-1, 0, PE::GetPE().GetCommunicator(m_namespace), MPI_STATUS_IGNORE);
        std::cout << "HSNBRadiator::setup() => P " << m_rank << " is set active. \n";
        //Copy datafiles into proc. directory
        std::string command2   = "cp -r " + m_HSNBPath.string() + "/data" + " " + m_hsnbDir.string();
        Common::OSystem::getInstance().executeCommand(command2);
        command2 = "cp -r " + Environment::DirPaths::getInstance().getWorkingDir().string() + "/hsnb.con" + " " + m_hsnbDir.string();
        Common::OSystem::getInstance().executeCommand(command2);

        m_HSNBPath = m_hsnbDir;

        if (m_rank!=m_nbProc) {
            MPI_Isend(&active, 1, MPIStructDef::getMPIType(&active), m_rank+1, 0, PE::GetPE().GetCommunicator(m_namespace),&sendRequest);
        }
        std::cout << "HSNBRadiator::setup() => P " << m_rank << " starts setting up files... \n";

    }

    boost::filesystem::path processFile("/hsnb.con");
    m_processFile = m_HSNBPath / processFile;

    boost::filesystem::path blackBodyFile("/data/Black_body_cumulated_energy.dat");
    m_blackBodyFile = m_HSNBPath / blackBodyFile;

    loadProcessData();




    //get the statesID used for this Radiator
    m_radPhysicsPtr->getCellStateIDs( m_statesID );

    Framework::DataSocketSink < Framework::State* , Framework::GLOBAL > states =
    m_radPhysicsHandlerPtr->getDataSockets()->states;
    DataHandle<State*, GLOBAL> statesHandle = states.getDataHandle();

    // set up the TRS states
    m_pstates.reset
        (new Framework::DofDataHandleIterator<CFreal, State, GLOBAL>(statesHandle, &m_statesID));




    //TODO Figure initialisation of matrices by ID

    if (statesHandle.size() != m_statesID.size()) {
      CFLog(ERROR, "HSNBRadiator::setup() => statesHandle.size() != m_statesID.size() => "
        << statesHandle.size() << " != " <<  m_statesID.size() << "\n");
      cf_assert(statesHandle.size() == m_statesID.size());
    }

    m_library = PhysicalModelStack::getActive()->getImplementor()->
      getPhysicalPropertyLibrary<PhysicalChemicalLibrary>();

    if (m_library.isNotNull()) {
      const CFuint nbSpecies = m_library->getNbSpecies();
      cf_assert(nbSpecies > 0);
      CFLog(INFO, "HSNBRadiator::setup() => nbSpecies = " << nbSpecies << "\n");
      m_mmasses.resize(nbSpecies);
      m_library->getMolarMasses(m_mmasses);


//      for (int i=0; i<nbSpecies; i++) {
//          CFLog(VERBOSE, "Molar mass of species " << i << "=" <<  std::setprecision(16) << m_mmasses[i] << "\n");
//      }



      // Needed to convert partial densities into [number_density]
      // X_i: Number density, rho_i partial density
      // X_i=rho_i*(avogadroOvMM)=rho_i*(Avogadro/M_i)=rho_i/m_molec_i
      m_avogadroOvMM.resize(nbSpecies);
      m_avogadroOvMM = PhysicalConsts::Avogadro()/m_mmasses;



    }



//    // set up the TRS states
//    m_pstates.reset(new Framework::DofDataHandleIterator<CFreal, Framework::State, Framework::GLOBAL>(statesHandle, &m_statesID));

//    CFreal *const currState = m_pstates->getState(0);
//    for (CFuint t = 0; t < 13; ++t) {
//      std::cout << currState[t] << std::endl;
//    }
//    std::cout << m_pI << " " << m_trI << std::endl;



    // Setup ThermoData class to administrate state / composition information



    setUpRadiationField();



    CFLog(INFO, "HSNBRadiator::setup => Loading BBCE data... \n");
    //Load black body cumulative energy
    load_BBCE_data();
    CFLog(INFO, "HSNBRadiator::setup => BBCE Data loaded. \n");



    setUpMechanisms();
    
    std::cout << "HSNBRadiator::setup => Rank " << m_rank << " finished setup of mechanisms. \n";

    m_curMechanismID=0;
    m_curMechanismRayCount=0;
    m_totalPower=0.0;
    m_atomicSelfAbsorption=0;
    m_avgStepDistance=0.0;
    m_nbSteps=0;

    //NOTE: ADD FLAG IN CONFIG FILE LATER
    m_tolerance = 0.001;




//    CFLog(VERBOSE, "HSNBRadiator::setup => specific emissive power " << emissivePowerDebug(0) << " \n");

//    while (true) {

//    }
    
}



void HSNBRadiator::moleFractionsFromPartialDensity(std::vector<CFreal> &convertVector)
{
 for (int i=0; i<convertVector.size(); i++) {
     convertVector[i]=convertVector[i]*this->m_avogadroOvMM[i];
 }
}

CFreal HSNBRadiator::getSpectraLoopPower()
{

    //NOTE: NOT SURE IF MAYBE THE CURRENT CELL HAS TO BE SET ELSEHOW
    const CFuint cellID =
      m_radPhysicsHandlerPtr->getCurrentCellStateID();

    //std::cout << "Radiator State " << cellID << std::endl;

    m_nbDiatomics=m_diatomics.size();
    m_nbContinua=m_continua.size();
    m_nbAtoms= m_atoms.size();

    CFLog(DEBUG_MAX, "HSNBRadiator::getSpectraLoopPower => m_nbDiatomics=" <<m_nbDiatomics
          << ", m_nbContinua=" <<m_nbContinua<< ", m_nbAtoms=" << m_nbAtoms << "\n");

    m_thermoData.setState(cellID);

    //m_thermoData.printState(cellID);

    m_tempLocalID=m_thermoData.getCurrentLocalCellID();

    cf_assert(m_radPhysicsHandlerPtr->getCurrentCellTrsIdx()==m_tempLocalID);

    m_emissionByMechanisms.clear();
    m_emissionByMechanisms.resize(m_nbDiatomics+m_nbContinua+(m_nbAtoms>0)+m_co2Exists, 0.0);

    CFLog(DEBUG_MAX, "HSNBRadiator::getSpectraLoopPower => m_emissionByMechanisms.size()=" << m_emissionByMechanisms.size() << "\n");


    m_totalEmissionCurrentCell=0.0;

    m_nbMechanisms=m_nbDiatomics+m_nbContinua+(m_nbAtoms>0)+m_co2Exists;


    for (int k = 0; k < m_nbDiatomics; ++k){
        m_emissionByMechanisms[k]=m_diatomics[k].emittedPower(m_tempLocalID,m_thermoData);
        m_totalEmissionCurrentCell+=m_emissionByMechanisms[k];
    }

    for (int k = m_nbDiatomics; k < m_nbContinua+m_nbDiatomics; ++k){
        m_emissionByMechanisms[k]=m_continua[k-m_nbDiatomics].emittedPower(m_thermoData, m_tempLocalID);
        m_totalEmissionCurrentCell+=m_emissionByMechanisms[k];
    }

    if (m_atoms.size() > 0) {
        m_emissionByMechanisms[m_nbDiatomics+m_nbContinua]=updateLineEmission();
        m_totalEmissionCurrentCell+=m_emissionByMechanisms[m_nbDiatomics+m_nbContinua];
    }

    if (m_co2Exists) {
        m_emissionByMechanisms[m_nbMechanisms-1]=m_co2->emittedPower(m_tempLocalID,m_thermoData);
        m_totalEmissionCurrentCell+=m_emissionByMechanisms[m_nbMechanisms-1];
    }



    //Volume in m3 convert to cm3
    CFreal currentCellPower=m_totalEmissionCurrentCell*getCurrentCellVolume() * 4.0* MathConsts::CFrealPi()*1e6;
    m_totalPower+=currentCellPower;
//    CFLog(INFO, "HSNBRadiator::getSpectraLoopPower => m_totalEmissionCurrentCell="<< m_totalEmissionCurrentCell <<", currentCellVolume= "<< getCurrentCellVolume()<< ", m_totalPower="<<m_totalPower <<" \n");
    //std::cout << "HSNBRadiator::getSpectraLoopPower => m_totalEmissionCurrentCell=" << m_totalEmissionCurrentCell << ", currentCellVolume= "<< getCurrentCellVolume()<< ", m_totalPower="<< m_totalPower <<" \n";
    // P=4*pi*V*sum_mechanisms(emission)
    return currentCellPower;

}

void HSNBRadiator::generateCellPhotonData(CFuint stateID, CFreal& energyFraction, CFuint &mechanismIndex,CFint& mechType, CFreal &lambda)
{
    // Keeps track of the number of rays assigned to each system in each cell
        CFreal waveNumber;
        //const CFuint cellID = m_radPhysicsHandlerPtr->getCurrentCellStateID();

        //Correct state has been set in setState
        m_thermoData.setState(stateID);
        m_tempLocalID=m_thermoData.getCurrentLocalCellID();


        ProbDistFunc line_pdf;

        CFLog(DEBUG_MAX, "HSNBRadiator::generateCellPhotonData => START, current cellID= "<< m_tempLocalID << "\n");

        // Update the line emission PDF if need be
        if (m_nbAtoms > 0 && m_nBRaysPerSys[m_nbContinua+m_nbDiatomics] > 0) {
            updateLineEmission();
            line_pdf = ProbDistFunc(m_line_emis);
        }


        //Skip over non-emitting processes
        while (m_curMechanismRayCount==m_nBRaysPerSys[m_curMechanismID]) {
          m_curMechanismID++;
          m_curMechanismRayCount=0;
        }

        mechanismIndex=m_curMechanismID;

//        std::cout << "HSNBRadiator::generateCellPhotonData => Mechanism Index " << mechanismIndex << std::endl;
//        std::cout << "m_curMechanismRayCount " << m_curMechanismRayCount << std::endl;
//        std::cout << "m_nBRaysPerSys[m_curMechanismID] " << m_nBRaysPerSys[m_curMechanismID] << std::endl;
//        std::cout << "m_emissionByMechanisms[m_curMechanismID] " << m_emissionByMechanisms[m_curMechanismID] << std::endl;

//        for (int i=0; i<m_nBRaysPerSys.size(); i++) {
//            std::cout << "m_nBRaysPerSys["<< i<<"] " << m_nBRaysPerSys[i] << std::endl;
//            std::cout << "m_emissionByMechanisms["<< i<<"] " << m_emissionByMechanisms[i] << std::endl;

//        }


        CFLog(DEBUG_MAX, "HSNBRadiator::generateCellPhotonData => Mechanism Index " << mechanismIndex << " \n");

        if (m_curMechanismID < m_nbMechanisms) {
            // Energy associated with each ray fired for this mechanism

            m_tempEnergyFraction=m_emissionByMechanisms[m_curMechanismID]/m_nBRaysPerSys[m_curMechanismID];
            CFLog(DEBUG_MAX, "HSNBRadiator::generateCellPhotonData => m_tempEnergyFraction=" << m_tempEnergyFraction << ", m_emissionByMechanisms[" << m_curMechanismID
                  << "]=" << m_emissionByMechanisms[m_curMechanismID] << " / m_nBRaysPerSys[" << m_curMechanismID <<"]=" << m_nBRaysPerSys[m_curMechanismID] << " \n");
            energyFraction =( m_tempEnergyFraction > 0. ) ? m_tempEnergyFraction : 0 ;


            //std::cout << "Energy fraction "<<energyFraction<<"\n";

            // Non-atomic emission first
            if (m_curMechanismID < m_nbDiatomics+m_nbContinua) {
                // Form the PDF for the band emission
                std::vector<double> band_emission;

                if (mechanismIsDiatomic(m_curMechanismID)) {

                    if (m_diatomics[m_curMechanismID].getMechanismType()==THIN) {
                        //Thin diatomic system: Nonthick
                        mechType=THINDIATOMIC;
                    }
                    else {
                        //CURTIS_GODSON or LINDQUIST_SIMMONS: Thick
                        mechType=THICKDIATOMIC;
                    }

//                    std::cout << "EMISSION BY " << m_diatomics[m_curMechanismID].speciesName() << ", energyFraction: " << energyFraction << std::endl;
                    band_emission.resize(m_diatomics[m_curMechanismID].spectralGridSize());
                    m_diatomics[m_curMechanismID].bandEmission(m_tempLocalID, m_thermoData, &band_emission[0]);
                } else {
                    //Continuum system: Nonthick
                    mechType=CONTINUUM;

//                    std::cout << "CONT EMISSION BY " << m_continua[m_curMechanismID-m_diatomics.size()].speciesName() << ", energyFraction: " << energyFraction << std::endl;
                    band_emission.resize(m_continua[m_curMechanismID-m_diatomics.size()].spectralGridSize());
                    m_continua[m_curMechanismID-m_diatomics.size()].bandEmission(m_thermoData,m_tempLocalID, &band_emission[0]);
                }

//                for (int i=0; i<band_emission.size(); i++) {
//                    CFLog(INFO, "HSNBRadiator::generateCellPhotonData => band_emission[" << i << "]=" << band_emission[i] << "\n");
//                }




                ProbDistFunc pdf(band_emission);

                // Fire rays in cell i for system k
                int band_offset = (m_curMechanismID < m_diatomics.size() ?
                    m_diatomics[m_curMechanismID].lowBand() :
                    m_continua[m_curMechanismID-m_diatomics.size()].lowBand());

//                CFLog(INFO, "HSNBRadiator::generateCellPhotonData => band_offset[" << m_curMechanismID << "]=" << band_offset << "\n");

//                    std::size_t pdeindex=pdf.index();
////                    double uniform = Random::uniform(0.000001, 999.999999);
//                    double uniform = Random::uniform(0.000001, 999.999999)-500.0;

//                    // Uniform random location in the band
                    waveNumber = (pdf.index()+band_offset)*1000.0 +
                        Random::uniform(0.000001, 999.999999)-500.0;

//                    waveNumber = (pdeindex+band_offset)*1000.0 + uniform;

                    //Force positive waveNumber (shouldn't occur anyway!)
//                    while (waveNumber<0.0) {
//                        waveNumber = (pdf.index()+band_offset)*1000.0 +
//                           Random::uniform(0.000001, 999.999999);
//                    }

//                    cout << "Non-Atomic wavenumber=" << waveNumber << ", pdeindex=" << pdeindex << ", uniform=" << uniform << ", band_offset=" << band_offset << endl;

                    while (wavenumberIsValid(waveNumber)==false) {
                           CFLog(DEBUG_MAX, "HSNBRadiator::generateCellPhotonData => Diatomic waveNumber " << waveNumber << " is invalid." "\n");


                        waveNumber = (pdf.index()+band_offset)*1000.0 +
                           Random::uniform(0.000001, 999.999999)-500.0;
                    }


            } else if (m_curMechanismID==(m_nbMechanisms-1-m_co2Exists)) {
                    //Atomic line mechanism
                    mechType=ATOMICLINE;

                    // Get line data for the selected line
                    int line_index = line_pdf.index();
                    int a = 0, end = m_atoms[0].nLines();

                    while (line_index >= end)
                        end += m_atoms[++a].nLines();
                    end -= m_atoms[a].nLines();  // point to beginning of atom

                    // Assuming the line data corresponds to cell i
                    const LineData& line =
                        m_atoms[a].lineData()[line_index - end];

                    waveNumber = Random::voigt(line.gaml, line.gamd, line.sigc);

//                    cout << "Atomic wavenumber=" << waveNumber << ", line.sigc=" << line.sigc << endl;


//                    while (waveNumber<0) {
//                        waveNumber = Random::voigt(line.gaml, line.gamd, line.sigc);
//                    }

                    while (wavenumberIsValid(waveNumber)==false) {
//                        cout << "NONVALID wavenumber=" << waveNumber << ", line.sigc=" << line.sigc << endl;

                        waveNumber = Random::voigt(line.gaml, line.gamd, line.sigc);
                    }


                }
                else { //CO2

                mechType=CO2;

                std::vector<double> band_emission;
                band_emission.resize(m_co2->spectralGridSize());
                m_co2->bandEmission(m_tempLocalID, m_thermoData, &band_emission[0]);

                ProbDistFunc pdf(band_emission);

                // Fire rays in cell i for system k
                int band_offset = m_co2->lowBand();

                // Uniform random location in the band
                waveNumber = (pdf.index()+band_offset)*25.0 +
                        Random::uniform(0.000001, 24.999999)-12.5;

                while (wavenumberIsValid(waveNumber)==false) {
                    CFLog(DEBUG_MAX, "HSNBRadiator::generateCellPhotonData => CO2 waveNumber " << waveNumber << " is invalid." "\n");

                    waveNumber = (pdf.index()+band_offset)*25.0 +
                            Random::uniform(0.000001, 24.999999)-12.5;
                }
            }

            }

         // cout << "MECHANISM ID " << m_curMechanismID << ", m_nBRaysPerSys[m_curMechanismID]: " << m_nBRaysPerSys[m_curMechanismID] << ", m_curMechanismRayCount " << m_curMechanismRayCount << " \n";

          m_curMechanismRayCount++;

          CFLog(DEBUG_MAX, "HSNBRadiator::generateCellPhotonData => waveNumber=" << waveNumber << "\n");


          lambda=1/waveNumber;

}

void HSNBRadiator::setState(CFuint nbRays)
{

    m_totalNbRaysCurrentCell=nbRays;

    //COOLFluiD is computing volumes in [m^3] but HSNB emission requires volumes in [cm^3]
    //thus we multiply with a factor of a million
    CFreal powerFactor= 4.0* MathConsts::CFrealPi()*getCurrentCellVolume()*1e6;;

    m_curMechanismID=0;
    m_curMechanismRayCount=0;
    m_nBRaysPerSys.clear();

    m_nbRaysPerCell=nbRays;

    const CFuint cellID =
      m_radPhysicsHandlerPtr->getCurrentCellStateID();

    m_nbDiatomics=m_diatomics.size();
    m_nbContinua=m_continua.size();
    m_nbAtoms= m_atoms.size();


    m_thermoData.setState(cellID);
    m_tempLocalID=m_thermoData.getCurrentLocalCellID();

    m_emissionByMechanisms.clear();
    m_nbMechanisms=m_nbDiatomics+m_nbContinua+(m_nbAtoms>0)+m_co2Exists;

    m_emissionByMechanisms.resize(m_nbMechanisms, 0.0);

   // CFreal emissionByAtomicLines;
    m_totalEmissionCurrentCell=0.0;


    for (int k = 0; k < m_nbDiatomics; ++k){
        m_emissionByMechanisms[k]=m_diatomics[k].emittedPower(m_tempLocalID,m_thermoData)*powerFactor;
        m_totalEmissionCurrentCell+=m_emissionByMechanisms[k];
    }

    //std::cout << "DIATOMICS FINISHED \n";

    for (int k = m_nbDiatomics; k < m_nbContinua+m_nbDiatomics; ++k){
        //std::cout << "k - m_nbDiatomics " << (k-m_nbDiatomics) << " m_nbContinua " << m_nbContinua << " \n";
        m_emissionByMechanisms[k]=m_continua[k-m_nbDiatomics].emittedPower(m_thermoData,m_tempLocalID)*powerFactor;
        m_totalEmissionCurrentCell+=m_emissionByMechanisms[k];
    }

    if (m_atoms.size() > 0) {
        m_emissionByMechanisms[m_nbDiatomics+m_nbContinua]=updateLineEmission()*powerFactor;
        m_totalEmissionCurrentCell+=m_emissionByMechanisms[m_nbDiatomics+m_nbContinua];
    }

    if (m_co2Exists) {
        m_emissionByMechanisms[m_nbMechanisms-1]=m_co2->emittedPower(m_tempLocalID,m_thermoData)*powerFactor;
        m_totalEmissionCurrentCell+=m_emissionByMechanisms[m_nbMechanisms-1];
    }

    CFLog(DEBUG_MAX, "HSNBRadiator::setState => m_totalEmissionCurrentCell=" << m_totalEmissionCurrentCell << "\n" );
    // Compute the number of rays associated with each system

    CFLog(DEBUG_MAX, "HSNBRadiator::setState => m_emissionByMechanisms.size()=" << m_emissionByMechanisms.size() << "\n" );

    m_nBRaysPerSys = ProbDistFunc(m_emissionByMechanisms).sample(nbRays);


//    for (int i=0; i<m_nBRaysPerSys.size(); i++) {
//        CFLog(VERBOSE, "HSNBRadiator::setState => m_nBRaysPerSys["<< i << "]=" << m_nBRaysPerSys[i] << "\n" );
//    }



//    //Preserve total energy
//    for (int i=0; i<m_nBRaysPerSys.size(); i++) {
//        //CFLog(INFO, "HSNBRadiator::setState => m_nBRaysPerSys[" << i << "]=" << m_nBRaysPerSys[i] << " / m_emissionByMechanisms["<< i << "]=" <<m_emissionByMechanisms[i] << "\n");
//        m_emissionByMechanisms[i]=m_totalEmissionCurrentCell*(CFreal(m_nBRaysPerSys[i])/CFreal(m_totalNbRaysCurrentCell));
//    }

//    cf_assert(sum==nbRays);


}

CFreal HSNBRadiator::computeAbsorbedEnergy(HSNBDataContainer &traceSet, CFint curCell)
{

    m_curPhoton=traceSet.headerData;

    // Photon energy has been spent
    if (m_curPhoton->energyResiduum<=m_tolerance) {
        m_curPhoton->energyResiduum=0.0;
        return 0.0;
    }

    //  No distance has been crossed yet within the mesh
    if (m_curPhoton->nbCrossedCells==0){
        return 0.0;
    }

//    m_thermoData.printState(curCell);
//    traceSet.print();

    CFLog(DEBUG_MAX, "HSNBRadiator::computeAbsorbedEnergy => nbCells=" << traceSet.trace->m_nbCellsCrossed << "\n");


    //m_curTrace=traceSet.trace.get();

    // Loop over each wall intersection
    m_tempTransp = computeTransmission(traceSet)*m_curPhoton->wallTransmissivity;



    cf_assert(!isnan(m_tempTransp));

    CFLog(DEBUG_MED, "HSNBRadiator::computeAbsorbedEnergy => m_tempTransp=" << m_tempTransp << "\n");

    // Absorb energy in current cell if tolerance is met
    if (m_tempTransp < m_tolerance) {
        //std::cout << "P:" << m_rank << " HSNBRadiator::computeAbsorbedEnergy => Lower than tolerance " << m_tolerance  << " (transm=" << m_tempTransp << ") \n";
        m_tempAbsorbedEnergy = m_curPhoton->energyFraction*m_curPhoton->transm;

        CFLog(DEBUG_MED, "HSNBRadiator::computeAbsorbedEnergy => Completely absorbed photon. m_tempAbsorbedEnergy= " << m_tempAbsorbedEnergy << ", Residuum=" << m_curPhoton->energyResiduum
            << ", transm=" << m_curPhoton->transm<< "\n");

        if (static_cast<HSNBMechanismType>(m_curPhoton->mechType)==ATOMICLINE || curCell==0) {
            m_atomicSelfAbsorption++;
        }

        //We are done with this photon, no more energy left
        m_curPhoton->energyResiduum=0.0;
        m_curPhoton->transm = m_tempTransp;
        CFLog(DEBUG_MED, "HSNBRadiator::computeAbsorbedEnergy => We want to return m_tempAbsorbedEnergy=" << m_tempAbsorbedEnergy << "\n");

        return m_tempAbsorbedEnergy;

    }

    CFLog(DEBUG_MED, "HSNBRadiator::computeAbsorbedEnergy => Partial absorbtion computed as. photonEnergy= "
         << m_curPhoton->energyFraction << ", Residuum=" << m_curPhoton->energyResiduum
         << ", transm=" << m_curPhoton->transm
         << ", transp" <<m_tempTransp<<  "\n");



    m_tempAbsorbedEnergy=m_curPhoton->energyFraction*(m_curPhoton->transm-m_tempTransp);
    m_curPhoton->transm = m_tempTransp;

    //Photon energy gets cut by the amount absorbed.
    m_curPhoton->energyResiduum-=m_tempAbsorbedEnergy;
    //std::cout << "P:" << m_rank << " Residuum: " << m_curPhoton->energyResiduum << " transp: "<<m_curPhoton->transm << " transp: " << m_tempTransp<< "\n";

   CFLog(DEBUG_MED, "HSNBRadiator::computeAbsorbedEnergy => Partially absorbed. absorbedEnergy= " << m_tempAbsorbedEnergy << ", Residuum=" << m_curPhoton->energyResiduum
         << ", transm" <<m_curPhoton->transm<<  "\n");


    return m_tempAbsorbedEnergy;

}

void HSNBRadiator::addStateParams(HSNBDataContainer& photonTraceSet, CFuint cellID, CFreal raydistance)
{

    //Careful: The raytracer computes distances in [m] but the HSNB radiator computes absorption in [cm]
    //Thus we have to convert stepDistance by multiplying with a factor of 100 once it's added to the trace.
    m_rayDistance=raydistance*100;

//    m_nbSteps++;
//    m_avgStepDistance+=(rayDistance-m_avgStepDistance)/m_nbSteps;
//    CFLog(INFO, "HSNBRadiator::addStateParams => stepDistance=" << rayDistance << " / avg=" << m_avgStepDistance << "\n");

    m_thermoData.setState(cellID);
    m_tempLocalID=m_thermoData.getCurrentLocalCellID();

    m_tempWavenumber=1/photonTraceSet.headerData->wavelength;



    if (photonTraceSet.trace->isUninitialised()) {

        if (static_cast<HSNBMechanismType>(photonTraceSet.headerData->mechType)==THICKDIATOMIC) {
            m_tempMechIndex=photonTraceSet.headerData->mechanismIndex;
            photonTraceSet.trace->m_kappa0=m_diatomics[m_tempMechIndex].getKappa(m_thermoData, m_tempLocalID, m_tempWavenumber);
//            std::cout <<  "HSNBRadiator::addStateParams => m_kappa0=" << photonTraceSet.trace->m_kappa0 <<std::endl;

        }

        photonTraceSet.trace->setup(m_nbContinua,m_nbAtoms,m_nbSpecies,m_co2Exists);
        initTrace(*photonTraceSet.trace, photonTraceSet.headerData->mechanismIndex, m_rayDistance);
        //std::cout << "NEW TRACE INITIALISED, m_nbThickDiatomics:" << m_nbThickDiatomics <<", m_nbNonThickDiatomics"
        //          << m_nbNonThickDiatomics << ", m_thickDiatomics.size()" << photonTraceSet.trace->m_thickDiatomics.size()  <<
        //          ", m_thinDiatomics.size() " << photonTraceSet.trace->m_thinDiatomics.size() << std::endl;

    }


//    photonTraceSet.trace->cellID.push_back(m_tempLocalID);

    photonTraceSet.trace->addTrackingState();
    photonTraceSet.headerData->nbCrossedCells++;


//    Used to time in place computation
//    Stopwatch<WallTime> stpStateAdd;
//    stpStateAdd.restart();

    for (int i=0; i<m_thinDiatomicsIndices.size(); i++) {
        m_tempNewStateMechIndex=m_thinDiatomicsIndices[i];
        m_diatomics[m_tempNewStateMechIndex].addStateParams(m_thermoData,photonTraceSet.trace->m_thinDiatomics[i], m_rayDistance, m_tempLocalID, m_tempWavenumber);
    }

//    std::cout << "HSNBRadiator::addStateParams => THIN ADD =" << stpStateAdd.read() << std::endl;
//    stpStateAdd.restart();

    for (int i=0; i<m_thickDiatomicsIndices.size(); i++) {
        m_tempNewStateMechIndex=m_thickDiatomicsIndices[i];
        m_diatomics[m_tempNewStateMechIndex].addStateParams(m_thermoData,photonTraceSet.trace->m_thickDiatomics[i], m_rayDistance, m_tempLocalID, m_tempWavenumber);
    }

//    std::cout << "HSNBRadiator::addStateParams => THICK ADD =" << stpStateAdd.read() << std::endl;
//    stpStateAdd.restart();

    for (int i=0; i<m_continua.size(); i++) {
        m_continua[i].addStateParams(m_thermoData,photonTraceSet.trace->m_continua[i], m_rayDistance, m_tempLocalID, m_tempWavenumber);

    }

//    std::cout << "HSNBRadiator::addStateParams => CONTINUA ADD =" << stpStateAdd.read() << std::endl;
//    stpStateAdd.restart();

    for (int i=0; i<m_atoms.size(); i++) {
        m_atoms[i].addStateParams(m_thermoData, photonTraceSet.trace->m_atoms[i], m_rayDistance, m_tempLocalID, m_tempWavenumber);
    }

//    std::cout << "HSNBRadiator::addStateParams => ATOMS ADD =" << stpStateAdd.read() << std::endl;
//    stpStateAdd.restart();

    //In case CO2 is present in the mixture
    if (m_co2Exists) {
//        photonTraceSet.trace->m_co2Exists=true;
        m_co2->addStateParams(m_thermoData, photonTraceSet.trace->m_co2, m_rayDistance, m_tempLocalID, m_tempWavenumber);
    }

//    std::cout << "HSNBRadiator::addStateParams => CO2 ADD =" << stpStateAdd.read() << std::endl;
//    stpStateAdd.stop();
//    photonTraceSet.trace->print();




}

void HSNBRadiator::setAbsorptionTolerance(const CFreal tol)
{
  cf_assert(tol>0.0);

  m_tolerance=tol;

}

void HSNBRadiator::initTrace(HSNBPhotonTrace& trace, CFuint mechID, CFreal startDistance)
{
    //distance0=distance0-ds
    trace.m_distance0=startDistance;
    trace.m_co2Exists=m_co2Exists;

    for (int i=0; i<m_diatomics.size(); i++) {
        if (m_diatomics[i].getMechanismType()==THIN) {
            trace.addThinDiatomicState((i==mechID));
        }
        else {
            trace.addThickDiatomicState((i==mechID));
        }
    }

    //std::cout << "P" << m_rank <<": HSNBRadiator::initTrace: ThinNb=" << trace.m_thinDiatomics.size() << ", ThickNb: " << trace.m_thickDiatomics.size() << " \n";
}

CFreal HSNBRadiator::computeTransmission(HSNBDataContainer &traceSet)
{
    //CFLog(DEBUG_MED, "HSNBRadiator::computeTransmission \n");


    m_tempEmittingMechanism=traceSet.headerData->mechanismIndex;
    //double sig = (1/traceSet.headerData->wavelength);
    CFLog(DEBUG_MED, "HSNBRadiator::computeTransmission => sig="<< 1/traceSet.headerData->wavelength << "\n");
    // Compute the total transmissivity for wavenumber sig

    //CFreal tempTau;
    m_tempTau = 0.0;
//    CFreal oldTau;
    for (int j = 0; j< m_thinDiatomicsIndices.size(); ++j) {
        m_tempCurMechIndex=m_thinDiatomicsIndices[j];

//        tempTau = m_diatomics[m_tempCurMechIndex].opticalThickness(traceSet.trace->m_thinDiatomics[j]);
//        CFLog(DEBUG_MED, "HSNBRadiator::computeTransmission => Thin::tau[" << j << "]=" << tempTau << " \n");

//        if (tempTau!=0.0) {
//            if (isnan(tempTau)) {
//                std::cout << "HSNBRadiator::computeTransmission => Tau is singular for thin diatomics= " <<  tempTau << "\n";
//                std::cout << "m_curPhoton->prevProcessId=" << m_curPhoton->prevProcessId << "\n";
//                std::cout << "My rank: " << m_rank << "\n";
//                traceSet.print(true);

//                while (true) {

//                }
//            }
//        }

//        oldTau=m_tempTau;

        m_tempTau += m_diatomics[m_tempCurMechIndex].opticalThickness(traceSet.trace->m_thinDiatomics[j]);
//        CFLog(INFO, "SnbDiatomicSystem::opticalThickness::Thin => kappa=" <<  m_tempTau-oldTau << " \n");

    }



    for (int j = 0; j< m_thickDiatomicsIndices.size(); ++j) {
        m_tempCurMechIndex=m_thickDiatomicsIndices[j];

//        tempTau = m_diatomics[m_tempCurMechIndex].opticalThickness(traceSet.trace->m_thickDiatomics[j]);
//        CFLog(DEBUG_MED, "HSNBRadiator::computeTransmission => Thick::tau[" << j << "]=" << tempTau << " \n");

//        if (tempTau!=0.0) {
//            if (isnan(tempTau)) {
//                std::cout << "HSNBRadiator::computeTransmission => Tau is singular for thick diatomics= " <<  tempTau << "\n";
//                std::cout << "m_curPhoton->prevProcessId=" << m_curPhoton->prevProcessId << "\n";
//                std::cout << "My rank: " << m_rank << "\n";
//                traceSet.print(true);

//                while (true) {

//                }
//            }
//        }



        m_tempTau += m_diatomics[m_tempCurMechIndex].opticalThickness(traceSet.trace->m_thickDiatomics[j]);

//        CFLog(INFO, "SnbDiatomicSystem::opticalThickness::Thick => kappa=" <<  m_tempTau-oldTau << " \n");

    }

    for (int j = 0; j < m_continua.size(); ++j) {
//        tempTau = m_continua[j].opticalThickness(traceSet.trace->m_continua[j]);
//        CFLog(DEBUG_MED, "HSNBRadiator::computeTransmission => Continua::tau[" << j << "]=" << tempTau << " \n");

//        if (tempTau!=0.0) {
//            if (isnan(tempTau)) {
//                std::cout << "HSNBRadiator::computeTransmission => Tau is singular for continuum= " <<  tempTau << "\n";
//                std::cout << "m_curPhoton->prevProcessId=" << m_curPhoton->prevProcessId << "\n";
//                std::cout << "My rank: " << m_rank << "\n";
//                traceSet.print(true);

//                while (true) {

//                }
//            }
//        }

        m_tempTau += m_continua[j].opticalThickness(traceSet.trace->m_continua[j]);
    }

    for (int j = 0; j < m_atoms.size(); ++j){
//        tempTau = m_atoms[j].opticalThickness(traceSet.trace->m_atoms[j]);
//        CFLog(DEBUG_MED, "HSNBRadiator::computeTransmission => Atoms::tau[" << j << "]=" << tempTau << " \n");

//        if (tempTau!=0.0) {
//            if (isnan(tempTau)) {
//                std::cout << "HSNBRadiator::computeTransmission => Tau is singular for atoms= " <<  tempTau << "\n";
//                std::cout << "atom index=" << j << "\n";
//                std::cout << "m_curPhoton->prevProcessId=" << m_curPhoton->prevProcessId << "\n";
//                std::cout << "My rank: " << m_rank << "\n";
//                traceSet.print(true);

//                while (true) {

//                }
//            }
//        }

        m_tempTau += m_atoms[j].opticalThickness(traceSet.trace->m_atoms[j]);
    }

    //"Add" absorption by CO2 systems which is computed only "stepwise", meaning that the
    // overall tau value has to be mulitplied in each each photon tracing step
    if (m_co2Exists) {
        //FALSCH muss mit gesamt-tau multi. werden
        m_tempTau += m_co2->opticalThickness(traceSet.trace->m_co2);
    }

//    if (m_co2Exists) {
//        m_tempTau+=m_co2->opticalThickness(traceSet.trace->m_co2);
//    }

//    CFLog(DEBUG_MED, "HSNBRadiator::computeTransmission => tau=" << m_tempTau << " TraceSet=" << traceSet.trace->m_nbCellsCrossed << " \n");

//    if (m_tempTau!=0.0) {
//        if (isnan(m_tempTau)) {
//            std::cout << "HSNBRadiator::computeTransmission => m_tempTau is singular= " <<  m_tempTau << "\n";
//            std::cout << "m_curPhoton->energyFraction=" << m_curPhoton->energyFraction << "\n";
//            std::cout << "m_curPhoton->energyResiduum=" << m_curPhoton->energyResiduum << "\n";
//            std::cout << "m_curPhoton->prevProcessId=" << m_curPhoton->prevProcessId << "\n";
//            std::cout << "My rank: " << m_rank << "\n";
//        }
//    }

    CFLog(DEBUG_MED, "HSNBRadiator::computeTransmission => TAU="<< m_tempTau<<  "\n");

    m_tempTau = std::exp(-m_tempTau);

    // Handle self absorption for thick molecular systems
    if (static_cast<HSNBMechanismType>(traceSet.headerData->mechType)!=THICKDIATOMIC) {
        return m_tempTau;
    }

    // Compute the tau for this system

  //  CFLog(DEBUG_MED, "HSNBRadiator::computeTransmission => Compute transmission for thick diatomic system. \n");

    m_tempTauk = m_diatomics[m_tempEmittingMechanism].tau(traceSet.trace->m_thickDiatomics[traceSet.trace->m_mechLocalID]);
    if (m_tempTauk == 1.0 || m_tempTauk == 0.0) {
        return m_tempTau;
    }

    // Build a new line of sight to account for ds
    m_tempds = 0.01*traceSet.trace->m_distance0;
    m_tempShortenedStartDistance=traceSet.trace->m_distance0-m_tempds;


    // Compute the tau for slightly shortened LOS
    m_Taukds=m_diatomics[m_tempEmittingMechanism].tauShortened(traceSet.trace->m_thickDiatomics[traceSet.trace->m_mechLocalID], traceSet.trace->m_kappa0, m_tempShortenedStartDistance);

//    //Don't forget to reset the pathlength
//    traceSet.trace->m_distance[0]+=ds;

  //  CFLog(DEBUG_MED, "HSNBRadiator::computeTransmission => Computed shortened tau. m_Taukds="<< m_Taukds << "\n");
    if (m_Taukds <= m_tempTauk) {
        return m_tempTau;
    }

    // Update the total transmissivity
    if (traceSet.trace->m_kappa0 > 0.0) {
        m_tempTau *= (m_Taukds - m_tempTauk)/(m_tempTauk*(traceSet.trace->m_kappa0)*m_tempds);
    }

    return m_tempTau;
}

bool HSNBRadiator::wavenumberIsValid(CFreal sig)
{
    //Wavenumber is either in range or no custom range is specified(for negative values of m_sigMin, use full available dataset in that case)
    if (sig<0) {
        return false;
    }
    else {
//        std::cout << "sig: " << sig <<  ", sigMax=" << m_sigMax << ", m_sigMin " << m_sigMin << std::endl;
        return (((sig>=m_sigMin) && (sig<=m_sigMax)) || m_sigMin<0);
    }
}

bool HSNBRadiator::co2Exists() const
{
    return m_co2Exists;
}

CFreal HSNBRadiator::getAbsorption(CFreal lambda, RealVector &s_o)
{
//    const double tolerance = 0.001;
//    double transm = 1.0, transp, walltrans = 1.0, wall_emis;

//    //REPLACE iw with a "currentWall"
//    int wall_pos = path.wallPosition(iw);

//    //LAST CELL

//    //std::cout << "cell " << ic << std::endl;
//    // Cell absorption
//    transp = computeTransmission(path, CURRENTCELL) * walltrans;

//    // Absorb energy in current cell
//    if (transp < tolerance) {
//        //SET FLAG: photon is extinct
//    }

//    transm = transp;

//    // Wall absorption
//    //std::cout << "wall " << iw << std::endl;
//    //XN: wall_emis=currentWallEmissivity
////    wall_emis = path.wallEmissivity(iw);
//    m_face_absorption[path.wallID(iw) - mp_mesh->nInternalFaces()] +=
//            path.energy() * transp * wall_emis;

//    // Update the wall transmissivity
//    walltrans *= (1.0 - wall_emis);
//    transm    *= (1.0 - wall_emis);

//    //Absorption: 1-tau,
    //    return (1-transm);
}

//Get Tau as described in Duarte2016
CFreal HSNBRadiator::getAbsorption(HSNBDataContainer &traceSet)
{
    m_tempTransp = computeTransmission(traceSet);
//    std::cout << "m_tempTransp=" << m_tempTransp << "\n";
    return -std::log(m_tempTransp);
}

CFreal HSNBRadiator::tau(HSNBDataContainer &traceSet)
{
    return -std::log(computeTransmission(traceSet));
}

void HSNBRadiator::absorb(HSNBPhotonTrace &trace, HSNBMechanismType mechType)
{
//    // Loop over each cell in the path along the way to the next wall
//    const double tolerance = 0.001;
//    double transm = 1.0, transp, walltrans = 1.0, wall_emis;

//    // Loop over each wall intersection
//    int ic = 0;
//    for (int iw = 0; iw < path.nWalls(); ++iw) {
//        int wall_pos = path.wallPosition(iw);
//        for (; ic < wall_pos; ++ic) {
//            //std::cout << "cell " << ic << std::endl;
//            // Cell absorption
//            transp = computeTransmission(path, ic+1) * walltrans;

//            // Absorb energy in current cell
//            if (transp < tolerance) {
//                m_cell_absorption[path.cellID(ic)] += path.energy()*transm;
//                //of << path.cellID(ic) << " " << path.spectralLocation() << " " << path.mechanism() << " " << path.energy() << "\n";
//                cell_extinctions++;
//                if (path.mechanism() >(m_continua.size()+m_diatomics.size()) || iw==0 || ic==0) {
//                    counter++;
//                }
//                return;
//            }
//            m_cell_absorption[path.cellID(ic)] += path.energy()*(transm-transp);
//            transm = transp;
//        }


//        // Wall absorption
//        //std::cout << "wall " << iw << std::endl;
//        wall_emis = path.wallEmissivity(iw);
//        m_face_absorption[path.wallID(iw) - mp_mesh->nInternalFaces()] +=
//            path.energy() * transp * wall_emis;

//        // Update the wall transmissivity
//        walltrans *= (1.0 - wall_emis);
//        transm    *= (1.0 - wall_emis);

//        if (transm < tolerance)
//            return;
//    }
}

CFreal HSNBRadiator::updateLineEmission(bool printDebug)
{

    CFreal sumEmissionOverLines=0.0;

    //Not necessary, correct state is precondition
    //m_thermoData.setState(cellID);

    int start = 0;
    for (int i = 0; i < m_atoms.size(); ++i) {
        // Compute emission (also updates line data)
        m_atoms[i].thinLineEmission(m_thermoData, &m_line_emis[start], sumEmissionOverLines, printDebug);
        start += m_atoms[i].nLines();
    }


    return sumEmissionOverLines;
}

bool HSNBRadiator::mechanismIsDiatomic(CFuint mechanismID)
{
    return mechanismID<m_nbDiatomics;
}

bool HSNBRadiator::mechanismIsContinous(CFuint mechanismID)
{
    return ((mechanismID>=m_nbDiatomics) && (mechanismID<m_nbDiatomics+m_nbContinua));
}

void HSNBRadiator::loadProcessData()
{

    boost::filesystem::path speciesDirectory;

    // is : ifstream
    // Load the atoms
    std::string line;
    std::vector<std::string> tokens;

    // Open the file
    std::ifstream& procFile = m_inFileHandle->open(m_processFile);


    // Read in the species data from the table header
    getline(procFile, line);
    String::tokenize(String::trim(line, " \t"), tokens, " \t");

    m_nbSpecies=tokens.size()-1;
    //Start with 1 (0 is the row title "Species:")
    for (int i = 1; i < tokens.size(); ++i) {
       m_thermoData.addSpecies(tokens[i]);
        //std::cout << tokens[i] << " " << m_thermoData.speciesIndex(tokens[i])<< " species" << std::endl;
    }

    CFLog(INFO, "HSNBRadiator::loadProcessData => m_sigMin="<< m_sigMin << " m_sigMax=" << m_sigMax << "]\n");

    m_thermoData.setSpectralRange(m_sigMin, m_sigMax);
    CFLog(INFO, "HSNB PATH " << m_HSNBPath.string() << "\n");

    //2nd row
    std::getline(procFile, line);

    m_co2Exists=m_thermoData.speciesExists("CO2");

    if (m_co2Exists) {
        if (!line.empty()) {
            String::tokenize(line, tokens, " \t\r\n");
            m_co2 = new SnbCO2System(m_HSNBPath.string(),tokens[4], m_thermoData);
            CFLog(INFO, "HSNBRadiator::loadProcessData => CO2 is present in the mixture. Path treatment: ["<< tokens[4] << "]\n");
        }
        else {
//            CFLog(INFO, "HSNBRadiator::loadProcessData: CO2 has been detected to be part of the mixture but your" <<
//                  " HSNB configuration does not specify a treatment for the CO2 species (Curtis-Godson [CG] or Linquist-Simmons [LS])." <<
//                  " The simulation will default to a CG treatment for CO2. \n");
//            m_co2 = new SnbCO2System(m_HSNBPath.string(),"CG", m_thermoData);

              CFLog(INFO, "HSNBRadiator::loadProcessData: CO2 has been detected to be part of the mixture but your" <<
                    " HSNB configuration does not specify a treatment for the CO2 species (Curtis-Godson [CG] or Linquist-Simmons [LS])." <<
                    " The simulation will ignore CO2 contributions. \n");
              m_co2Exists=false;
        }
    }

    CFLog(INFO, line << "\n");

    std::getline(procFile, line);
    CFLog(INFO, line <<"\n");

    while (line.empty()) {
        std::getline(procFile, line);
    }
    //Row to define atomic species
    std::getline(procFile, line);
    CFLog(INFO, "HSNBRadiator::loadProcessData => Atomic species: " << line << "\n");


    if (!line.empty()) {
        String::tokenize(line, tokens, " \t\r\n");

        for (int i = 0; i < tokens.size(); ++i) {
            m_atoms.push_back(AtomicLines(m_HSNBPath.string(), tokens[i], m_thermoData));
        }

        std::getline(procFile, line);
    }


    //Skip to Diatomic Header
    while (line.empty()) {
        std::getline(procFile, line);
    }

    std::vector<SnbDiatomicSystem>::iterator it;
    // Rows to define diatomic species
    // Load diatomics
    std::getline(procFile, line);
    m_nbNonThickDiatomics=0;
    m_nbThickDiatomics=0;
    while (!line.empty()) {
        tokens.clear();
        String::tokenize(line, tokens, " \t\r\n");

        m_diatomics.push_back(
            SnbDiatomicSystem(SpeciesLoadData(m_HSNBPath.string(),tokens[0]), tokens[1], tokens.size() > 2 ? tokens[2] : "",m_thermoData));
        it=m_diatomics.end()-1;


        if (m_diatomics.back().getMechanismType()==THIN){
            m_thinDiatomicsIndices.push_back(m_diatomics.size()-1);
            m_nbNonThickDiatomics++;
        }
        else {
            m_thickDiatomicsIndices.push_back(m_diatomics.size()-1);
            m_nbThickDiatomics++;
        }

        getline(procFile, line);
    }

    //Skip to Continuum Header
    while (line.empty()) {
        std::getline(procFile, line);
    }

    // Load continua
    getline(procFile, line);
    while (!line.empty()) {
        tokens.clear();
        String::tokenize(line, tokens, " \t\r\n");

        if (tokens[1] == "BF")
            m_continua.push_back(SnbContinuumSystem(SpeciesLoadData(m_HSNBPath.string(),tokens[0]), m_thermoData, BOUNDFREE));
        else
            m_continua.push_back(SnbContinuumSystem(SpeciesLoadData(m_HSNBPath.string(),tokens[0]), m_thermoData, FREEFREE));
        getline(procFile, line);
    }


    CFuint curMechIndex;
    for (int i=0; i<m_thinDiatomicsIndices.size(); i++) {
        curMechIndex=m_thinDiatomicsIndices[i];
        CFLog(DEBUG_MIN, "HSNBRadiator::loadProcessData: m_thinDiatomics["<<i<<"].NUP=" << m_diatomics[curMechIndex].getMechanismType() << " \n");
    }

    for (int i=0; i<m_thickDiatomicsIndices.size(); i++) {
        curMechIndex=m_thickDiatomicsIndices[i];
        CFLog(DEBUG_MIN, "HSNBRadiator::loadProcessData: m_thickDiatomics["<<i<<"].NUP=" << m_diatomics[curMechIndex].getMechanismType() << " \n");
    }

    for (int i=0; i<m_diatomics.size(); i++) {
        CFLog(DEBUG_MIN, "HSNBRadiator::loadProcessData: m_diatomics["<<i<<"].NUP=" << m_diatomics[i].getMechanismType() << " \n");
    }


//    for (int i=0; i<m_diatomics.size(); i++) {
//        cout << "m_diatomics(" << i << ")="<< m_diatomics[i].systemName() << ", lowBand=" << m_diatomics[i].lowBand() << endl;
//    }
//    for (int i=0; i<m_continua.size(); i++) {
//        cout << "m_diatomics(" << i << ")="<< m_diatomics[i].systemName() << ", lowBand=" << m_diatomics[i].lowBand() << endl;
//    }
//    while (true){

//    }


    m_inFileHandle->close();



}

CFreal HSNBRadiator::emissivePowerDebug(CFuint stateID)
{
    //NOTE: NOT SURE IF MAYBE THE CURRENT CELL HAS TO BE SET ELSEHOW
    const CFuint cellID = stateID;

    m_thermoData.printState(stateID);
    //std::cout << "Radiator State " << cellID << std::endl;

    m_nbDiatomics=m_diatomics.size();
    m_nbContinua=m_continua.size();
    m_nbAtoms= m_atoms.size();

    m_nbMechanisms=m_nbDiatomics+m_nbContinua+(m_nbAtoms>0)+m_co2Exists;

    m_thermoData.setState(cellID);
    m_emissionByMechanisms.clear();
    m_emissionByMechanisms.resize(m_nbMechanisms, 0.0);

    m_totalEmissionCurrentCell=0.0;



    for (int k = 0; k < m_nbDiatomics; ++k){
        m_emissionByMechanisms[k]=m_diatomics[k].emittedPower(cellID,m_thermoData);
        m_totalEmissionCurrentCell+=m_emissionByMechanisms[k];
        CFLog(INFO, "HSNBRadiator::emissivePowerDebug => m_emissionByMechanisms[" << k << "]="<<m_emissionByMechanisms[k]<<" \n");
    }

    for (int k = m_nbDiatomics; k < m_nbContinua+m_nbDiatomics; ++k){
        m_emissionByMechanisms[k]=m_continua[k-m_nbDiatomics].emittedPower(m_thermoData,cellID);
        m_totalEmissionCurrentCell+=m_emissionByMechanisms[k];
    }

    if (m_atoms.size() > 0) {
        m_emissionByMechanisms[m_nbDiatomics+m_nbContinua]=updateLineEmission();
        m_totalEmissionCurrentCell+=m_emissionByMechanisms[m_nbDiatomics+m_nbContinua];
    }

    if (m_co2Exists) {
        m_emissionByMechanisms[m_nbMechanisms-1]=m_co2->emittedPower(m_tempLocalID,m_thermoData);
        m_totalEmissionCurrentCell+=m_emissionByMechanisms[m_nbMechanisms-1];
    }

    CFreal currentCellSpecificPower=m_totalEmissionCurrentCell;


    // P=4*pi*V*sum_mechanisms(emission)
    return currentCellSpecificPower;
}




/////////////////////////////////////////////////////////////////////////////

}
}
