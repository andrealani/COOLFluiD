#include <fstream>
#include <iostream>

#include "Common/PE.hh"
#include "Common/BadValueException.hh"
#include "Common/CFPrintContainer.hh"

#include "MathTools/MathConsts.hh"

#include "Environment/ObjectProvider.hh"
#include "Environment/CFEnv.hh"
#include "Environment/FileHandlerOutput.hh"
#include "Environment/FileHandlerInput.hh"
#include "Environment/DirPaths.hh"
#include "Environment/SingleBehaviorFactory.hh"

#include "Framework/PathAppender.hh"
#include "Framework/DataProcessing.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"
#include "Framework/PhysicalChemicalLibrary.hh"
#include "Framework/SocketBundleSetter.hh"

#include "FiniteVolume/CellCenterFVM.hh"

#include "RadiativeTransfer/RadiativeTransfer.hh"
#include "RadiativeTransfer/Solvers/FiniteVolumeDOM/RadiativeTransferFVDOM.hh"
#include "RadiativeTransfer/RadiationLibrary/Radiator.hh"
#include "RadiativeTransfer/RadiationLibrary/RadiationPhysicsHandler.hh"

/////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Numerics::FiniteVolume;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace RadiativeTransfer {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<RadiativeTransferFVDOM, DataProcessingData, RadiativeTransferModule>
radiativeTransferFVDOMProvider("RadiativeTransferFVDOM");

//////////////////////////////////////////////////////////////////////////////

RadiativeTransferFVDOM::RadiativeTransferFVDOM(const std::string& name) :
  DataProcessingCom(name),
  socket_states("states"),
  socket_gstates("gstates"),
  socket_nstates("nstates"),
  socket_volumes("volumes"),
  socket_nodes("nodes"),
  socket_isOutward("isOutward"),
  socket_normals("normals"),
  socket_faceCenters("faceCenters"),
  socket_CellID("CellID"),
  socket_divq("divq"),
  socket_qx("qx"),
  socket_qy("qy"),
  socket_qz("qz"),
  socket_TempProfile("TempProfile"),
  socket_alpha_avbin("alpha_avbin"),
  socket_B_bin("B_bin"),
  socket_qradFluxWall("qradFluxWall"),
  m_library(CFNULL),
  m_radiation(new RadiationPhysicsHandler("RadiationPhysicsHandler")),
  m_mapGeoToTrs(CFNULL),
  m_isWallFace(),
  m_mapWallTRSToOffsetFaceID(),
  m_geoBuilder(), 
  m_wallFaceBuilder(),
  m_normal(),
  m_weight(),
  m_fieldSource(),
  m_fieldAbsor(),
  m_fieldAbSrcV(),
  m_fieldAbV(),
  m_In(),
  m_II(),
  m_opacities(),
  m_radSource(),
  m_Ttable(),
  m_Ptable(),
  m_dotProdInFace(),
  m_sdone(),
  m_cdoneIdx(),
  m_wallTrsNames(),
  m_dirs(),
  m_advanceOrder(),
  m_q(),
  m_qrAv(),
  m_divqAv(),
  m_nbDirTypes()
{
  addConfigOptionsTo(this);
  
  m_nbDirs = 8;
  this->setParameter("nDirs",&m_nbDirs);
  
  m_dirName = Environment::DirPaths::getInstance().getWorkingDir();
  setParameter("DirName", &m_dirName);
  
  m_binTabName = "air-100.dat";
  setParameter("BinTabName", &m_binTabName); 
  
  m_outTabName = "air-100.out";
  setParameter("OutTabName", &m_outTabName);
  
  m_writeToFile = true;
  setParameter("WriteToFile", &m_writeToFile);
    
  m_useExponentialMethod = true;
  setParameter("UseExponentialMethod", &m_useExponentialMethod);
  
  m_radialData = false;
  setParameter("RadialData", &m_radialData);
  
  m_oldAlgo = true;
  setParameter("OldAlgorithm", &m_oldAlgo);
  
  m_Nr = 100;
  setParameter("NbRadialPoints", &m_Nr);
  
  m_constantP = -1.; // negative by default to force to give pressure if needed
  setParameter("ConstantP", &m_constantP);
  
  m_Tmin = 1000.;
  setParameter("Tmin", &m_Tmin);
  
  m_Tmax = 12000.;
  setParameter("Tmax", &m_Tmax);

  m_deltaT = 0.0071;
  setParameter("DeltaT", &m_deltaT);
  
  m_PID = 0;
  setParameter("PID", &m_PID);
  
  m_TID = 4;
  setParameter("TID", &m_TID);
  
  m_nbThreads = 1;
  setParameter("NbThreads", &m_nbThreads);
  
  m_threadID = 0;
  setParameter("ThreadID", &m_threadID);
  
  m_loopOverBins = true;
  setParameter("LoopOverBins", &m_loopOverBins);
  
  m_emptyRun = false;
  setParameter("EmptyRun", &m_emptyRun);
  
  m_readOpacityTables = true;
  setParameter("ReadOpacityTables", &m_readOpacityTables);
}
      
//////////////////////////////////////////////////////////////////////////////

RadiativeTransferFVDOM::~RadiativeTransferFVDOM()
{
}

//////////////////////////////////////////////////////////////////////////////

void RadiativeTransferFVDOM::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFuint >
    ("nDirs","Number of directions. Only allowed 8, 24, 48 and 80");
  options.addConfigOption< boost::filesystem::path >
    ("DirName","Name of the directories where the .dat file is located.");
  options.addConfigOption< string >
    ("BinTabName","Name of the .dat file");
  options.addConfigOption< bool >
    ("WriteToFile","Writing the table in ASCII");   
  options.addConfigOption< string >
    ("OutTabName","Name of the output file");    
  options.addConfigOption< bool >
    ("UseExponentialMethod","Exponential method for radiation. Explained in ICCFD7-1003");
  options.addConfigOption< bool >
    ("RadialData","radial q and divQ for the sphere case");
  options.addConfigOption< bool >
    ("OldAlgorithm","old algorithm (very inefficient) just kept for comparison purposes");
  options.addConfigOption< CFuint >
    ("NbRadialPoints","Number of Radial points in the sphere case");
  options.addConfigOption< CFreal >("ConstantP","Constant pressure temperature");
  options.addConfigOption< CFreal >("Tmin","Minimum temperature");
  options.addConfigOption< CFreal >("Tmax","Maximum temperature");
  options.addConfigOption< CFreal >("DeltaT","Temperature interval");
  options.addConfigOption< CFuint >("PID","ID of pressure in the state vector");
  options.addConfigOption< CFuint >("TID","ID of temperature in the state vector");
  options.addConfigOption< CFuint >("NbThreads","Number of threads/CPUs in which the algorithm has to be split.");
  options.addConfigOption< CFuint >("ThreadID","ID of the current thread within the parallel algorithm."); 
  options.addConfigOption< bool >("LoopOverBins","Loop over bins and then over directions (do the opposite if =false).");
  options.addConfigOption< bool >("EmptyRun","Run without actually solving anything, just for testing purposes.");
  options.addConfigOption< bool >("ReadOpacityTables","Read the opacity tables instead if using the radiation library.");
}
      
//////////////////////////////////////////////////////////////////////////////

void RadiativeTransferFVDOM::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;
  
  CFLog(VERBOSE, "RadiativeTransferFVDOM::configure() => START\n");
  
  DataProcessingCom::configure(args);
  
  cf_assert(m_radiation.isNotNull());
  configureNested ( m_radiation.getPtr(), args );
  
  CFLog(VERBOSE, "RadiativeTransferFVDOM::configure() => END\n");
}
      
//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSource> >
RadiativeTransferFVDOM::providesSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSource> > result;

  result.push_back(&socket_CellID);
  result.push_back(&socket_divq);
  result.push_back(&socket_qx);
  result.push_back(&socket_qy);
  result.push_back(&socket_qz);
  result.push_back(&socket_TempProfile);
  result.push_back(&socket_alpha_avbin);
  result.push_back(&socket_B_bin);
  result.push_back(&socket_qradFluxWall);
  
  return result;
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
RadiativeTransferFVDOM::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;
  
  result.push_back(&socket_states);
  result.push_back(&socket_gstates);
  result.push_back(&socket_nstates);
  result.push_back(&socket_volumes);
  result.push_back(&socket_nodes);
  result.push_back(&socket_isOutward);
  result.push_back(&socket_normals);  
  result.push_back(&socket_faceCenters);  
  
  return result;
}

//////////////////////////////////////////////////////////////////////////////

void RadiativeTransferFVDOM::setup()
{
  CFAUTOTRACE;
  
  DataProcessingCom::setup();
  
  SocketBundle sockets;
  sockets.states      = socket_states;
  sockets.gstates     = socket_gstates;
  sockets.nstates     = socket_nstates;
  sockets.volumes     = socket_volumes;
  sockets.nodes       = socket_nodes;
  sockets.isOutward   = socket_isOutward;
  sockets.normals     = socket_normals;
  sockets.faceCenters = socket_faceCenters;
  
  // source sockets cannot be copied
  sockets.alpha_avbin = socket_alpha_avbin.getDataHandle();
  sockets.B_bin       = socket_B_bin.getDataHandle();
  
  // sockets.faceCenters = socket_faceCenters;
  // sockets.faceAreas   = socket_faceAreas;
  
  m_radiation->setupDataSockets(sockets);
  m_radiation->setup();
  m_radiation->configureTRS();
  // m_radiation->setupAxiFlag(m_isAxi);
  // vector<string> boundaryTrsNames;
  // m_radiation->getBoundaryTRSnames(boundaryTrsNames);
  
  // cell builder
  m_geoBuilder.getGeoBuilder()->setDataSockets(socket_states, socket_gstates, socket_nodes);
  m_geoBuilder.setup();
  CellTrsGeoBuilder::GeoData& geoData = m_geoBuilder.getDataGE();
  geoData.trs = MeshDataStack::getActive()->getTrs("InnerCells");
  
  // wall face builder
  m_wallFaceBuilder.setup();
  m_wallFaceBuilder.getGeoBuilder()->setDataSockets(socket_states, socket_gstates, socket_nodes);
  m_wallFaceBuilder.getDataGE().isBFace = true;
  
  m_mapGeoToTrs = MeshDataStack::getActive()->getMapGeoToTrs("MapFacesToTrs");
  
  CFLog(VERBOSE, "RadiativeTransferFVDOM::setup() => start\n");
  CFLog(VERBOSE, "RadiativeTransferFVDOM::setup() => [threadID / nbThreads] = [" 
	<< m_threadID << " / " << m_nbThreads << "]\n");
  
  cf_assert(m_PID < PhysicalModelStack::getActive()->getNbEq());
  cf_assert(m_TID < PhysicalModelStack::getActive()->getNbEq());
  cf_assert(m_PID != m_TID);
  
  const std::string nsp = getMethodData().getNamespace();
  cf_assert(PE::GetPE().GetProcessorCount(nsp) == 1);
  
  // Setting up the file containing the binary table with opacities
  CFLog(VERBOSE, "RadiativeTransferFVDOM::setup() => m_dirName = "<< m_dirName <<"\n"); 
  CFLog(VERBOSE, "RadiativeTransferFVDOM::setup() => m_binTabName = "<< m_binTabName <<"\n");
  
  m_inFileHandle = Environment::SingleBehaviorFactory<Environment::FileHandlerInput>::getInstance().create();
  m_binTableFile = m_dirName / boost::filesystem::path(m_binTabName);
  
  CFLog(VERBOSE, "RadiativeTransferFVDOM::setup() => m_binTableFile = "<< m_binTableFile <<"\n");
  
  // Check to force an allowed number of directions
  if ((m_nbDirs != 8) && (m_nbDirs != 24) && (m_nbDirs != 48) && (m_nbDirs != 80)) {
    CFLog(WARN, "RadiativeTransferFVDOM::setup() => This ndirs is not allowed. 8 directions is chosen \n");
    m_nbDirs = 8;
  } 
  
  if (m_readOpacityTables) {
    // Reading the table
    readOpacities();
  }
  
  const CFuint DIM = 3;
  
  // Selecting the # of direction types depending on the option
  switch(m_nbDirs) {
  case 8:     
    m_nbDirTypes = 1;
    break;
  case 24:
    m_nbDirTypes = 1;    
    break;
  case 48:
    m_nbDirTypes = 2;
    break;
  case 80:
    m_nbDirTypes = 3;
    break;
  default:		//Case nDirs == 8
    m_nbDirTypes = 1;    
    break;
  }
  
  //Resizing the vectors
  m_weight.resize(m_nbDirs);
  
  // Get number of cells
  Common::SafePtr<Common::ConnectivityTable<CFuint> > cells =
    MeshDataStack::getActive()->getConnectivity("cellStates_InnerCells");

  const CFuint nbCells = cells->nbRows();
    
  m_sdone.resize(nbCells);
  m_cdoneIdx.reserve(nbCells);
  
  if(m_useExponentialMethod){
    m_fieldSource.resize(nbCells);
    m_fieldAbsor.resize(nbCells);
  }
  else{
    m_fieldSource.resize(nbCells);
    m_fieldAbSrcV.resize(nbCells);
    m_fieldAbV.resize(nbCells);
  }
  m_In.resize(nbCells);
  m_II.resize(nbCells);
  
  for (CFuint iCell = 0; iCell < nbCells; iCell++){
    m_sdone[iCell] = false;
  }
  
  // m_nbThreads, m_threadID
  cf_assert(m_nbDirs > 0);
  cf_assert(m_nbBins > 0);
  
  // set the start/end bins for this process
  if (m_nbThreads == 1) { 
    m_startEndDir.first  = 0;
    m_startEndBin.first  = 0;
    m_startEndDir.second = m_nbDirs-1;
    m_startEndBin.second = m_nbBins-1;
  }
  else {
    const CFuint nbBinDir     = m_nbBins*m_nbDirs; 
    const CFuint minNbThreadsPerProc = nbBinDir/m_nbThreads;
    const CFuint maxNbThreadsPerProc = minNbThreadsPerProc + nbBinDir%m_nbThreads;
    cf_assert(minNbThreadsPerProc > 0);
    
    // same direction has same meshdata structure, therefore if you have 
    // m_nbThreads <= m_nbDirs it is more scalable to split by direction 
    const CFuint startThread = m_threadID*minNbThreadsPerProc;
    const CFuint nbThreadsPerProc = (m_threadID < m_nbThreads-1) ? 
      minNbThreadsPerProc : maxNbThreadsPerProc; 
    const CFuint endThread = startThread + nbThreadsPerProc;
    
    CFLog(VERBOSE, "RadiativeTransferFVDOM::setup() => nbBins, nbDirs   = [" << m_nbBins << ", " << m_nbDirs << "]\n");
    CFLog(VERBOSE, "RadiativeTransferFVDOM::setup() => startThread      = [" << startThread << "]\n");
    CFLog(VERBOSE, "RadiativeTransferFVDOM::setup() => endThread        = [" << endThread << "]\n");
    CFLog(VERBOSE, "RadiativeTransferFVDOM::setup() => nbThreadsPerProc = [" << nbThreadsPerProc << "]\n");
    
    if (m_loopOverBins) {
      // suppose you sweep all entries while walking row-wise through 
      // a matrix (b,d) of size m_nbBins*m_nbDirs, consider the corresponding 
      // (b_start,d_start) and (b_end,d_end) points
      m_startEndDir.first = startThread%m_nbDirs; 
      m_startEndBin.first = startThread/m_nbDirs;
      m_startEndDir.second = (endThread-1)%m_nbDirs;
      m_startEndBin.second = (endThread-1)/m_nbDirs;
    }
    else {
      // suppose you sweep all entries while walking row-wise through 
      // a matrix (d,b) of size m_nbDirs*m_nbBins, consider the corresponding 
      // (d_start,b_start) and (d_end,b_end) points
      m_startEndDir.first = startThread/m_nbBins; 
      m_startEndBin.first = startThread%m_nbBins;
      m_startEndDir.second = (endThread-1)/m_nbBins;
      m_startEndBin.second = (endThread-1)%m_nbBins;
    }
  }
  
  const CFuint startBin = m_startEndBin.first;
  const CFuint endBin   = m_startEndBin.second+1;
  cf_assert(endBin <= m_nbBins);
  const CFuint startDir = m_startEndDir.first;
  const CFuint endDir   = m_startEndDir.second+1;
  cf_assert(endDir <= m_nbDirs);

  if (m_loopOverBins) {
    CFLog(INFO, "RadiativeTransferFVDOM::setup() => start/end Bin = [" << startBin << ", " << endBin << "]\n");
    CFLog(INFO, "RadiativeTransferFVDOM::setup() => start/end Dir = [" << startDir << ", " << endDir << "]\n");
    CFLog(VERBOSE, "RadiativeTransferFVDOM::setup() => full (Bin, Dir) list: \n");
    
    for(CFuint ib = startBin; ib < endBin; ++ib) {
      const CFuint dStart = (ib != startBin) ? 0 : startDir;
      const CFuint dEnd   = (ib != m_startEndBin.second) ? m_nbDirs : endDir;
      for (CFuint d = dStart; d < dEnd; ++d) {
	CFLog(VERBOSE, "(" << ib << ", " << d <<"), ");
      }
      CFLog(VERBOSE, "\n");
    }
    CFLog(VERBOSE, "\n");
  }
  else {
    CFLog(INFO, "RadiativeTransferFVDOM::setup() => start/end Dir = [" << startDir << ", " << endDir << "]\n");
    CFLog(INFO, "RadiativeTransferFVDOM::setup() => start/end Bin = [" << startBin << ", " << endBin << "]\n");
    CFLog(VERBOSE, "RadiativeTransferFVDOM::setup() => full (Dir, Bin) list: \n");
    
    for(CFuint d = startDir; d < endDir; ++d) {
      const CFuint bStart = (d != startDir) ? 0 : startBin;
      const CFuint bEnd   = (d != m_startEndDir.second) ? m_nbBins : endBin;
      for (CFuint ib = bStart; ib < bEnd; ++ib) {
	CFLog(VERBOSE, "(" << d << ", " << ib <<"), ");
      }
      CFLog(VERBOSE, "\n");
    }
    CFLog(VERBOSE, "\n");
  }
  
  m_dirs.resize(m_nbDirs, 3);
  m_advanceOrder.resize(m_nbDirs);
  m_q.resize(nbCells, DIM);
  
  cf_assert(endDir <= m_nbDirs);
  // resize only the rows corresponding to considered directions 
  for (CFuint i = startDir; i< endDir; i++) {
    m_advanceOrder[i].resize(nbCells);
  }
  
  m_normal.resize(DIM, 0.); 
  
  DataHandle<CFreal> CellID = socket_CellID.getDataHandle();
  DataHandle<CFreal> divQ = socket_divq.getDataHandle();
  DataHandle<CFreal> qx = socket_qx.getDataHandle();
  DataHandle<CFreal> qy = socket_qy.getDataHandle();
  DataHandle<CFreal> qz = socket_qz.getDataHandle();
  DataHandle<CFreal> TempProfile = socket_TempProfile.getDataHandle();
  divQ.resize(nbCells);
  CellID.resize(nbCells);
  TempProfile.resize(nbCells);
  qx.resize(nbCells);
  qy.resize(nbCells);
  qz.resize(nbCells);
  divQ   = 0.0;
  qx     = 0.0;
  qy     = 0.0;
  qz     = 0.0;
  CellID = 0.0;
  TempProfile = 0.0;
  
  // resize the bins storage
  socket_alpha_avbin.getDataHandle().resize(nbCells*m_nbBins);
  socket_B_bin.getDataHandle().resize(nbCells*m_nbBins);
  
  // the following returns a list of all the .ApplyTRS with .TypeTRS=Wall in the .CFcase
  // for which radiative heat flux has to be computed 
  m_radiation->getWallTRSnames(m_wallTrsNames);
  // FaceTrsGeoBuilder::GeoData& wallFacesData = m_wallFaceBuilder.getDataGE();

  const CFuint totalNbFaces = MeshDataStack::getActive()->Statistics().getNbFaces();
  cf_assert(socket_normals.getDataHandle().size()/3);
  m_isWallFace.resize(totalNbFaces, false);
  m_dotProdInFace.resize(totalNbFaces);
  
  CFuint nbFaces = 0; // total number of boundary faces belonging to TRS of type "Wall"  
  for(CFuint j=0; j< m_wallTrsNames.size(); ++j) {
    const string wallTRSName = m_wallTrsNames[j];
    SafePtr<TopologicalRegionSet> wallFaces = MeshDataStack::getActive()->getTrs(wallTRSName);
    CFLog(VERBOSE, "RadiativeTransferFVDOM::setup() => TRS["<< wallTRSName << "] is Wall\n");
    const CFuint nbLocalFaces = wallFaces->getLocalNbGeoEnts();
    for (CFuint f = 0; f < nbLocalFaces; ++f) {
      m_isWallFace[wallFaces->getLocalGeoID(f)] = true;
    }
    // store the offset for the corresponding wall TRS name
    m_mapWallTRSToOffsetFaceID[wallTRSName] = nbFaces;
    nbFaces += nbLocalFaces;
  }
  // preallocation of memory for qradFluxWall
  socket_qradFluxWall.getDataHandle().resize(nbFaces);
  
  //Averages for the Sphere case
  if (m_radialData) {
    m_qrAv.resize(m_Nr);
    m_divqAv.resize(m_Nr);
    m_qrAv   = 0;
    m_divqAv = 0;  
  }
  
  Stopwatch<WallTime> stp;
  
  stp.start();
  getDirections();
  CFLog(INFO, "RadiativeTransferFVDOM::setup() => getDirections() took " << stp.read() << "s\n");
  
  stp.start();
  
  if (!m_emptyRun) {
    // only get advance order for the considered directions
    for (CFuint d = startDir; d < endDir; d++){
      cf_assert(m_advanceOrder[d].size() == nbCells);
      getAdvanceOrder(d, m_advanceOrder[d]);
    }
  }
    
  CFLog(INFO, "RadiativeTransferFVDOM::setup() => getAdvanceOrder() took " << stp.read() << "s\n");
    
  CFLog(VERBOSE, "RadiativeTransferFVDOM::setup() => end\n");
}
      
//////////////////////////////////////////////////////////////////////////////

void RadiativeTransferFVDOM::getDirections()
{
  CFLog(DEBUG_MIN, "RadiativeTransferFVDOM::getDirections() => start\n");
  
  const CFreal pi = MathTools::MathConsts::CFrealPi();
  
  CFLog(VERBOSE, "RadiativeTransferFVDOM::getDirections() => Number of Directions = " << m_nbDirs << "\n");
  
  const CFreal overSq3 = 1./std::sqrt(3.);
  switch(m_nbDirs) {
  case 8:
    m_weight[0] = 4.*pi/m_nbDirs;
    m_dirs(0,0) = overSq3;
    m_dirs(0,1) = overSq3;
    m_dirs(0,2) = overSq3;
    break;
  case 24:
    m_weight[0] = 4.*pi/m_nbDirs;
    m_dirs(0,0) = 0.2958759;
    m_dirs(0,1) = 0.2958759;
    m_dirs(0,2) = 0.9082483;
    break;
  case 48:
    m_weight[0] = 0.1609517;
    m_dirs(0,0) = 0.1838670;
    m_dirs(0,1) = 0.1838670;
    m_dirs(0,2) = 0.9656013;
    m_weight[1] = 0.3626469;
    m_dirs(1,0) = 0.1838670;
    m_dirs(1,1) = 0.6950514;
    m_dirs(1,2) = 0.6950514;
    break;
  case 80:
    m_weight[0] = 0.1712359;
    m_dirs(0,0) = 0.1422555;
    m_dirs(0,1) = 0.1422555;
    m_dirs(0,2) = 0.9795543;
    m_weight[1] = 0.0992284;
    m_dirs(1,0) = 0.1422555;
    m_dirs(1,1) = overSq3;
    m_dirs(1,2) = 0.8040087;
    m_weight[2] = 0.4617179;
    m_dirs(2,0) = overSq3;
    m_dirs(2,1) = overSq3;
    m_dirs(2,2) = overSq3;
    break;
  default:	// nDirs = 8
    m_weight[0] = 4.*pi/m_nbDirs;
    m_dirs(0,0) = overSq3;
    m_dirs(0,1) = overSq3;
    m_dirs(0,2) = overSq3;
    break;
  }
  
  CFuint d = m_nbDirTypes - 1; //Note that it has been changed, because the counters start at 0
  for (CFuint dirType = 0; dirType < m_nbDirTypes; dirType++){
    for (CFuint p = 0; p <= 2; p++){
      const CFuint l = p;	    //Note that it's different because the counter starts at 0
      const CFuint m = (p+1) % 3; //Note a % b is the remainder of the division a/b
      const CFuint n = (p+2) % 3;
      
      if (p == 0 || m_dirs(dirType,0) != m_dirs(dirType,1) ||
	  m_dirs(dirType,1) != m_dirs(dirType,2) || m_dirs(dirType,2) != m_dirs(dirType,0)) {
        CFLog(VERBOSE, "Case1::dirTypes = " << dirType <<"\n");
	CFLog(DEBUG_MIN, "l = " << l << "m = " << m << "n = " << n  <<"\n");
	for (int i = 0; i <= 1; i++) {
	  for (int j = 0; j <= 1; j++) {
	    for (int k = 0; k <= 1; k++) {
	      if ( p+i+j+k != 0) {
		//Note that this is different because the counters are different
		d += 1;
		m_weight[d] = m_weight[dirType];
		m_dirs(d,0) = std::pow(-1.,i)*m_dirs(dirType,l);
		m_dirs(d,1) = std::pow(-1.,j)*m_dirs(dirType,m);
		m_dirs(d,2) = std::pow(-1.,k)*m_dirs(dirType,n);
		CFLog(DEBUG_MIN, "l = " << l << " m = " << m << " n = " << n  <<"\n");
		CFLog(DEBUG_MIN, "d = " << d <<"\n");
		CFLog(DEBUG_MIN, "dirs[" << d <<"] = ("<<  m_dirs(d,0) <<", " << m_dirs(d,1) <<", "<<m_dirs(d,2)<<")\n");
	      }
	    }
	  }
	}
      }     
      if (m_dirs(dirType,0) != m_dirs(dirType,1) && m_dirs(dirType,1) != m_dirs(dirType,2) 
	  && m_dirs(dirType,2) != m_dirs(dirType,0)) {
	CFLog(VERBOSE, "Case2::dirTypes = " << dirType <<"\n");
	CFLog(DEBUG_MIN, "l = " << l << "m = " << m << "n = " << n  <<"\n");
	for (int i = 0; i <= 1; i++) {
	  for (int j = 0; j <= 1; j++) {
	    for (int k = 0; k <= 1; k++) {
	      //Note that this is different because the counters are different
	      d += 1;
	      m_weight[d] = m_weight[dirType];
	      m_dirs(d,0) = std::pow(-1.,i)*m_dirs(dirType,l);
	      m_dirs(d,1) = std::pow(-1.,j)*m_dirs(dirType,m);
	      m_dirs(d,2) = std::pow(-1.,k)*m_dirs(dirType,n);
	      CFLog(DEBUG_MIN, "l = " << l << " m = " << m << " n = " << n  <<"\n");
	      CFLog(DEBUG_MIN, "d = " << d <<"\n");
	      CFLog(DEBUG_MIN, "dirs[" << d <<"] = ("<<  m_dirs(d,0) <<", " << m_dirs(d,1) <<", "<<m_dirs(d,2)<<")\n");
	    }
	  }
	}
      }          
    }
  }
  
  // Printing the Directions for debugging
  for (CFuint dir = 0; dir < m_nbDirTypes; dir++) {
    CFLog(DEBUG_MIN, "Direction[" << dir <<"] = (" << m_dirs(dir,0) <<", " << m_dirs(dir,1) <<", " << m_dirs(dir,2) <<")\n");
  }
  
  CFLog(DEBUG_MIN, "RadiativeTransferFVDOM::getDirections() => end\n");
}

//////////////////////////////////////////////////////////////////////////////
      
void RadiativeTransferFVDOM::execute()
{
  CFAUTOTRACE;
  
  CFLog(INFO, "RadiativeTransferFVDOM::execute() => start\n");
  
  Stopwatch<WallTime> stp;
  
  stp.start();
  
  // compute the spectra in steps (to avoid to have to store too much data at once): 
  // this calls Radiator::setupSpectra(wavMin, wavMax)
  const CFuint nbLoops = m_radiation->getNumberLoops();
  for(CFuint i= 0; i < nbLoops; ++i) {
    m_radiation->setupWavStride(i);
  }
  
  CFLog(INFO, "RadiativeTransferFVDOM::execute() => radiation library took " << stp.read() << "s\n");
  
  stp.start();
  
  if (!m_emptyRun) {
    DataHandle<CFreal> divQ = socket_divq.getDataHandle();
    
    // only one CPU allow for namespace => the mesh has not been partitioned
    cf_assert(PE::GetPE().GetProcessorCount(getMethodData().getNamespace()) == 1);
    
    // Compute the order of advance
    // Call the function to get the directions
    m_q  = 0.0;
    divQ = 0.0;
    m_II = 0.0;
    
    const CFuint startBin = m_startEndBin.first;
    const CFuint endBin   = m_startEndBin.second+1;
    cf_assert(endBin <= m_nbBins);
    const CFuint startDir = m_startEndDir.first;
    const CFuint endDir   = m_startEndDir.second+1;
    cf_assert(endDir <= m_nbDirs);
    
    if (m_loopOverBins) {
      loopOverBins(startBin, endBin, startDir, endDir);
    }
    else {
      loopOverDirs(startBin, endBin, startDir, endDir);
    }
    
    DataHandle<CFreal> volumes = socket_volumes.getDataHandle();
   
    DataHandle<CFreal> qx   = socket_qx.getDataHandle();
    DataHandle<CFreal> qy   = socket_qy.getDataHandle();
    DataHandle<CFreal> qz   = socket_qz.getDataHandle();
    
    // const string fileName = "divq-" + StringOps::to_str(PE::GetPE().GetRank("Default"));
    // ofstream fout(fileName.c_str());
    
    SafePtr<TopologicalRegionSet> cells = m_geoBuilder.getDataGE().trs;
    const CFuint nbCells = cells->getLocalNbGeoEnts();
    for (CFuint iCell = 0; iCell < nbCells; iCell++) {
      divQ[iCell] /= volumes[iCell]; //converting area from m^3 into cm^3
      
      // if (iCell == 1000){fout << "iCell1000 => " << divQ[iCell] << "\n";}
      
      qx[iCell] = m_q(iCell,XX);
      qy[iCell] = m_q(iCell,YY);
      qz[iCell] = m_q(iCell,ZZ);
    }
        
    if (m_radialData){
      writeRadialData();
    } 
  }
    
  CFLog(INFO, "RadiativeTransferFVDOM::execute() => took " << stp.read() << "s \n");
}
      
//////////////////////////////////////////////////////////////////////////////

void RadiativeTransferFVDOM::getFieldOpacities(CFuint ib)
{
  CFLog(VERBOSE, "RadiativeTransferFVDOM::getFieldOpacities() => start\n");
  
  m_fieldSource = 0.;
  if(m_useExponentialMethod) {
    m_fieldAbsor  = 0.;
  }
  else {
    m_fieldAbSrcV = 0.;
    m_fieldAbV    = 0.;
  }
  
  DataHandle<CFreal> TempProfile    = socket_TempProfile.getDataHandle();
  DataHandle<CFreal> volumes        = socket_volumes.getDataHandle();
  DataHandle<State*, GLOBAL> states = socket_states.getDataHandle();
  
  RadiativeTransferFVDOM::Interpolator interp;
  
  const CFuint nbCells = states.size();
  cf_assert(nbCells > 0);
  const CFuint totalNbEqs = PhysicalModelStack::getActive()->getNbEq();
  
  for (CFuint iCell = 0; iCell < nbCells; iCell++) {
    const State *currState = states[iCell];
    //Get the field pressure and T commented because now we impose a temperature profile
    CFreal p = 0.;
    CFreal T = 0.;
    
    if (m_constantP < 0.) {
      cf_assert(m_PID < currState->size());
      cf_assert(m_TID < currState->size());
      p = (*currState)[m_PID];
      T = (*currState)[m_TID];
    }
    else {
      cf_assert(m_constantP > 0.);
      cf_assert(m_Tmax > 0.);
      cf_assert(m_Tmin > 0.);
      cf_assert(m_deltaT > 0.);
      
      const CFreal r  = currState->getCoordinates().norm2();
      p = m_constantP; //Pressure in Pa
      const CFreal Tmax   = m_Tmax;
      const CFreal Tmin   = m_Tmin;
      const CFreal deltaT = m_deltaT;
      const CFreal A      = std::pow(r*0.01/deltaT, 2);
      const CFreal rmax   = 1.5;
      const CFreal Amax   = std::pow(rmax*0.01/deltaT, 2);
      T = Tmax - (Tmax - Tmin)*(1 - std::exp(-A))/(1 - std::exp(-Amax));
    }
    
    TempProfile[iCell] = T;
    
    const CFreal patm   = p/101325.; //converting from Pa to atm
    CFreal val1 = 0;
    CFreal val2 = 0;
    
    interp.tableInterpolate(m_nbBins, m_nbTemp, m_nbPress, 
			    &m_Ttable[0], &m_Ptable[0],
			    &m_opacities[0], &m_radSource[0],
			    T, patm, ib, val1, val2); 
    
    if(m_useExponentialMethod){
      if (val1 <= 1e-30 || val2 <= 1e-30 ){
	m_fieldSource[iCell] = 1e-30;
	m_fieldAbsor[iCell]  = 1e-30;
      }
      else {
	m_fieldSource[iCell] = val2/val1;
	m_fieldAbsor[iCell]  = val1;
      }
    }
    else{
      if (val1 <= 1e-30 || val2 <= 1e-30 ){
	m_fieldSource[iCell] = 1e-30;
	m_fieldAbV[iCell]    = 1e-30*volumes[iCell]; // Volumen converted from m^3 into cm^3
      }
      else {
	m_fieldSource[iCell] = val2/val1;
	m_fieldAbV[iCell]    = val1*volumes[iCell];
      }      
      m_fieldAbSrcV[iCell]   = m_fieldSource[iCell]*m_fieldAbV[iCell];
    }
  }
  
  CFLog(VERBOSE, "RadiativeTransferFVDOM::getFieldOpacities() => end\n");
}

//////////////////////////////////////////////////////////////////////////////

void RadiativeTransferFVDOM::getAdvanceOrder(const CFuint d, 
					     LocalArray<CFint>::TYPE& advanceOrder)
{
  CFLog(VERBOSE, "RadiativeTransferFVDOM::getAdvanceOrder() => start\n");
  
  // The order of advance calculation begins here
  DataHandle<CFreal> CellID = socket_CellID.getDataHandle();
  const CFuint nbCells = CellID.size();
  cf_assert(nbCells > 0);
  const CFuint DIM = PhysicalModelStack::getActive()->getDim();
  cf_assert(DIM == DIM_3D);
  
  CFLog(INFO, "RadiativeTransferFVDOM::getAdvanceOrder() => Direction number [" << d <<"]\n");
  
  // precompute the dot products for all faces and directions (a part from the sign)
  computeDotProdInFace(d, m_dotProdInFace);
  
  CFuint mLast = 0;
  CFuint m = 0;

  CFuint stage = 1;
  // Initializing the sdone and cdone for each direction
  m_sdone.assign(nbCells, false);
  
  m_cdoneIdx.clear();
  cf_assert(m_cdoneIdx.size() == 0);
  cf_assert(m_cdoneIdx.capacity() == nbCells);
  
  SafePtr<ConnectivityTable<CFuint> > cellFaces = MeshDataStack::getActive()->getConnectivity("cellFaces");
  DataHandle<CFint> isOutward = socket_isOutward.getDataHandle();
    
  while (m < nbCells) { //The loop over the cells begins
    mLast = m;	  //Checking to see if it counts all the cells
    for (CFuint iCell = 0; iCell < nbCells; iCell++) {
      CFLog(DEBUG_MAX, "RadiativeTransferFVDOM::getAdvanceOrder() => iCell = " << iCell <<"\n");
      if (m_sdone[iCell] == false) {
	const CFuint nbFaces = cellFaces->nbCols(iCell);
	for (CFuint iFace = 0; iFace < nbFaces; ++iFace) {
	  const CFuint faceID = (*cellFaces)(iCell, iFace);
	  const bool isBFace  = m_mapGeoToTrs->isBGeo(faceID);
	  const bool neighborIsSdone = (!isBFace) ? m_sdone[getNeighborCellID(faceID, iCell)] : true;
	  
	  // When dot product is < 0 and neighbor is undone, skip the cell and continue the loop
	  const CFreal factor = ((CFuint)(isOutward[faceID]) != iCell) ? -1. : 1.;
	  const CFreal dotMult = m_dotProdInFace[faceID]*factor;
	  if (dotMult < 0. && neighborIsSdone == false) {goto cell_loop;}
	}// end loop over the FACES
	
	// CFLog(DEBUG_MAX, "advanceOrder[" << d << "][" << m <<"] = " << iCell << "\n");
	
	advanceOrder[m] = iCell;
	CellID[iCell] = stage;
	m_cdoneIdx.push_back(iCell);
	m += 1;
      }// end if(Cell is not done)
      
    cell_loop:
      const bool dummy = true;
    }// end of the loop over the CELLS
    
    const string msg = "advanceOrder[" + StringOps::to_str(d) + "] = ";
#ifndef CF_HAVE_CUDA
    CFLog(DEBUG_MAX, CFPrintContainer<LocalArray<CFint>::TYPE>(msg, &advanceOrder) << "\n");
#else
    CFint *const aoPtr = advanceOrder.ptr(); 
    CFLog(DEBUG_MAX, msg << "\n");
    for (CFuint a = 0; a < advanceOrder.size(); ++a) {
      CFLog(DEBUG_MAX, aoPtr[a] << " ");
    }
    CFLog(DEBUG_MAX, "\n");
#endif
    
    if (m == mLast) {diagnoseProblem(d, m, mLast);}
    
    advanceOrder[m - 1] *= -1;
    
    CFLog(DEBUG_MAX, "At stage[" << stage << "] => m_cdoneIdx.size() = " << m_cdoneIdx.size() << "\n");
    for (CFuint id = 0; id < m_cdoneIdx.size(); ++id) {
      m_sdone[m_cdoneIdx[id]] = true;
    }
    m_cdoneIdx.clear();
    
    CFLog(VERBOSE, "RadiativeTransferFVDOM::getAdvanceOrder() => m  "<< m << " \n");
    CFLog(VERBOSE, "RadiativeTransferFVDOM::getAdvanceOrder() => End of the "<< stage << " stage\n");
    
    ++stage;
  }// end of the loop over the STAGES
  
  //Printing advanceOrder for debug purpuses
  const string msg = "RadiativeTransferFVDOM::getAdvanceOrder() => advanceOrder[" + StringOps::to_str(d) + "] = ";

#ifndef CF_HAVE_CUDA
    CFLog(DEBUG_MAX, CFPrintContainer<LocalArray<CFint>::TYPE>(msg, &advanceOrder) << "\n");
#else
    CFint *const aoPtr = advanceOrder.ptr(); 
    CFLog(DEBUG_MAX, msg << "\n");
    for (CFuint a = 0; a < advanceOrder.size(); ++a) {
      CFLog(DEBUG_MAX, aoPtr[a] << " ");
    }
    CFLog(DEBUG_MAX, "\n");
#endif
    
  CFLog(VERBOSE, "RadiativeTransferFVDOM::getAdvanceOrder() => end\n");
}
      
//////////////////////////////////////////////////////////////////////////////

void RadiativeTransferFVDOM::readOpacities()
{
  CFLog(VERBOSE, "RadiativeTransferFVDOM::readOpacities() => start\n");
  
  fstream& is = m_inFileHandle->openBinary(m_binTableFile);
  
  // the three first numbers are #bins, #Temps, #pressures
  vector<double> data(3);
  is.read((char*)&data[0], 3*sizeof(double));
  
  m_nbBins  = ((int) data[0]);
  m_nbTemp  = ((int) data[1]);
  m_nbPress = ((int) data[2]);
  
  vector<double> Pressures(m_nbPress);
  is.read((char*)&Pressures[0], m_nbPress*sizeof(double));
  
  vector<double> Temperatures(m_nbTemp);
  is.read((char*)&Temperatures[0], m_nbTemp*sizeof(double));
  
  // in the table there are 143 zeros 
  vector<double> Zeros(143);
  is.read((char*)&Zeros[0], 143*sizeof(double));
  
  // Reading the table: for each pressure
  // each bin: value1(Temp1) value2(Temp1) ... value1(TempN) value2(TempN)
  vector< vector<double> > bins(m_nbPress*m_nbBins,vector<double>(2*m_nbTemp)); //initialize  
  for (CFuint ib =0; ib < m_nbPress*m_nbBins; ++ib){
    is.read((char*)&bins[ib][0], 2*m_nbTemp*sizeof(double));
  }  
  
  double end;
  is.read((char*)&end, sizeof(double));
  
  // Verifying that the last number is zero, so the table is finished
  cf_assert(int(end) == 0);
  
  is.close();
  //}
  CFLog(VERBOSE, "RadiativeTransferFVDOM::readOpacities() => After closing the binary File\n");
  
  // Storing the opacities and the radiative source into memory
  m_opacities.resize(m_nbPress*m_nbBins*m_nbTemp);
  m_radSource.resize(m_nbPress*m_nbBins*m_nbTemp);
  
  // Setting up the temperature and pressure tables
  m_Ttable.resize(m_nbTemp);
  m_Ptable.resize(m_nbPress);
  
  for(CFuint ib = 0; ib < m_nbBins; ++ib) {
    for(CFuint ip = 0; ip < m_nbPress; ++ip) {
      for(CFuint it = 0; it < m_nbTemp; ++it) {
	m_opacities[it + ib*m_nbTemp + ip*m_nbBins*m_nbTemp] = bins[ip + ib*m_nbPress][2*it];
	m_radSource[it + ib*m_nbTemp + ip*m_nbBins*m_nbTemp] = bins[ip + ib*m_nbPress][2*it + 1];
      }
    }
  }
  for(CFuint it = 0; it < m_nbTemp; ++it){
    m_Ttable[it] = Temperatures[it];
  }
  
  for(CFuint ip = 0; ip < m_nbPress; ++ip){
    m_Ptable[ip] = Pressures[ip];
  }
  
  // Writing the table into a file
  if(m_writeToFile){
    CFLog(VERBOSE, "RadiativeTransferFVDOM::readOpacities() => Writing file \n");
    boost::filesystem::path file = m_dirName / boost::filesystem::path(m_outTabName);
    file = PathAppender::getInstance().appendParallel( file );
    
    SelfRegistPtr<Environment::FileHandlerOutput> fhandle = Environment::SingleBehaviorFactory<Environment::FileHandlerOutput>::getInstance().create();
    ofstream& fout = fhandle->open(file);
    
    fout <<"#Bins = "<< data[0] <<"\t#Temps = "<< data[1] <<"\t#Pressures = "<< data[2] << endl;
    
    fout << endl;
    fout <<"Pressures[atm]  = ";
    for(CFuint ip = 0; ip < m_nbPress; ++ip){
      fout << Pressures[ip] << " ";
    }
    fout << endl;
    fout << endl;
    fout <<"Temperatures[K] = ";
    for(CFuint it = 0; it < m_nbTemp; ++it){
      fout << Temperatures[it] << " ";
    }
    fout << endl;
    CFuint m = 0;
    fout.precision(18);
    for(CFuint ip = 0; ip < m_nbPress; ++ip){
      fout << endl;
      fout <<"Pressure = "<< Pressures[ip] << endl;
      fout <<"bin \t\t\t\t Temp \t\t\t\t val1 \t\t\t\t\t\t val2" << endl;
      CFuint beginLoop = m*ip;
      for(CFuint ib = 0; ib < m_nbBins; ++ib){
	for(CFuint it = 0; it < m_nbTemp; ++it){
	  fout << ib + 1 <<"\t\t\t\t" << Temperatures[it] <<"\t\t\t\t" 
	       << m_opacities[it + ib*m_nbTemp + ip*m_nbBins*m_nbTemp] 
	       <<"\t\t\t\t" << m_radSource[it + ib*m_nbTemp + ip*m_nbBins*m_nbTemp] << endl;
	}
      }
      m += m_nbBins; 
    }
    fhandle->close();     
  }
  
  CFLog(VERBOSE, "RadiativeTransferFVDOM::readOpacities() => end\n");
}

//////////////////////////////////////////////////////////////////////////////

void RadiativeTransferFVDOM::writeRadialData()
{
  CFLog(VERBOSE, "RadiativeTransferFVDOM::writeRadialData() = > Writing radial data for the spherical test case => start\n");
  
  boost::filesystem::path file = m_dirName / boost::filesystem::path("radialData.plt");
  file = PathAppender::getInstance().appendParallel( file );
  
  SelfRegistPtr<Environment::FileHandlerOutput> fhandle = 
    Environment::SingleBehaviorFactory<Environment::FileHandlerOutput>::getInstance().create();
  ofstream& outputFile = fhandle->open(file);
  
  outputFile << "TITLE  = RadiativeTransferFVDOM radial data for a sphere\n";
  outputFile << "VARIABLES = r  qr divq nbPoints\n";

  CellTrsGeoBuilder::GeoData& geoData = m_geoBuilder.getDataGE();
  Common::SafePtr<TopologicalRegionSet> cells = geoData.trs;
  
  const CFuint nbCells = cells->getLocalNbGeoEnts();
  CFuint nbPoints = 0;
  CFreal rCoord = 0.;
  const CFreal Radius = 1.5;
  DataHandle<CFreal> divQ = socket_divq.getDataHandle();
    
  for(CFuint ir = 0; ir < m_Nr; ir++){
    nbPoints = 0;
    rCoord = (ir + 0.5)*Radius/m_Nr; //middle point between ir and (ir + 1)
    
    for(CFuint iCell = 0; iCell < nbCells; iCell++){
      geoData.idx = iCell;
      GeometricEntity* currCell = m_geoBuilder.buildGE();
      
      const Node& coordinate = currCell->getState(0)->getCoordinates();
      const CFreal x = coordinate[XX];
      const CFreal y = coordinate[YY];
      const CFreal z = coordinate[ZZ];
      const CFreal rCell = std::sqrt(x*x + y*y + z*z);
      
      if(rCell >= ir*Radius/m_Nr && rCell < (ir + 1)*Radius/m_Nr){
	nbPoints++;
	m_divqAv[ir] += divQ[iCell];
	m_qrAv[ir]   += (m_q(iCell,XX)*x + m_q(iCell,YY)*y + m_q(iCell,ZZ)*z)/rCell; //*rCell*rCell; Multiply by r**2 for area-weighted average
      }
      m_geoBuilder.releaseGE();
    }
    if(nbPoints > 0){
      m_divqAv[ir] /= nbPoints;
      m_qrAv[ir]   /= nbPoints; //m_qrAv[ir]   /= nbPoints*rCoord*rCoord; //area-weighted average radial flux
      outputFile << rCoord << " " << m_qrAv[ir] << " " << m_divqAv[ir] << " " <<  nbPoints << "\n";
    }
  }
  fhandle->close();
  
  CFLog(VERBOSE, "RadiativeTransferFVDOM::writeRadialData() => Writing radial data for the spherical test case => end\n");
}

//////////////////////////////////////////////////////////////////////////////

void RadiativeTransferFVDOM::unsetup()
{
  CFAUTOTRACE;
  
  DataProcessingCom::unsetup();
}
      
//////////////////////////////////////////////////////////////////////////////

void RadiativeTransferFVDOM::getFieldOpacities(const CFuint ib, const CFuint iCell)
{
  m_fieldSource[iCell] = 0.;
  if(m_useExponentialMethod){
    m_fieldAbsor[iCell] = 0.;
  }
  else{
    m_fieldAbSrcV[iCell] = 0.;
    m_fieldAbV[iCell]    = 0.;
  }
  
  DataHandle<CFreal> TempProfile = socket_TempProfile.getDataHandle();
  DataHandle<CFreal> volumes     = socket_volumes.getDataHandle();
  
  const State *currState = socket_states.getDataHandle()[iCell]; 
  //Get the field pressure and T commented because now we impose a temperature profile
  CFreal p = 0.;
  CFreal T = 0.;
  if (m_constantP < 0.) {
    cf_assert(m_PID < currState->size());
    cf_assert(m_TID < currState->size());
    p = (*currState)[m_PID];
    T = (*currState)[m_TID];
  }
  else {
    cf_assert(m_constantP > 0.);
    cf_assert(m_Tmax > 0.);
    cf_assert(m_Tmin > 0.);
    cf_assert(m_deltaT > 0.);
    
    const CFreal r  = currState->getCoordinates().norm2();
    p = m_constantP; //Pressure in Pa
    const CFreal Tmax   = m_Tmax;
    const CFreal Tmin   = m_Tmin;
    const CFreal deltaT = m_deltaT;
    const CFreal A      = std::pow(r*0.01/deltaT, 2);
    const CFreal rmax   = 1.5;
    const CFreal Amax   = std::pow(rmax*0.01/deltaT, 2);
    T = Tmax - (Tmax - Tmin)*(1 - std::exp(-A))/(1 - std::exp(-Amax));
  }
  
  TempProfile[iCell] = T;
  
  const CFreal patm   = p/101325.; //converting from Pa to atm
  CFreal val1 = 0;
  CFreal val2 = 0;

  RadiativeTransferFVDOM::Interpolator interp;
  interp.tableInterpolate(m_nbBins, m_nbTemp, m_nbPress, 
			  &m_Ttable[0], &m_Ptable[0],
			  &m_opacities[0], &m_radSource[0],
			  T, patm, ib, val1, val2); 
  
  if(m_useExponentialMethod){
    if (val1 <= 1e-30 || val2 <= 1e-30 ){
      m_fieldSource[iCell] = 1e-30;
      m_fieldAbsor[iCell]  = 1e-30;
    }
    else {
      m_fieldSource[iCell] = val2/val1;
      m_fieldAbsor[iCell]  = val1;
    }
  }
  else{
    if (val1 <= 1e-30 || val2 <= 1e-30 ){
      m_fieldSource[iCell] = 1e-30;
      m_fieldAbV[iCell]    = 1e-30*volumes[iCell]; // Volume converted from m^3 into cm^3
    }
    else {
      m_fieldSource[iCell] = val2/val1;
      m_fieldAbV[iCell]    = val1*volumes[iCell];
    }      
    m_fieldAbSrcV[iCell]   = m_fieldSource[iCell]*m_fieldAbV[iCell];
  }
}
      
//////////////////////////////////////////////////////////////////////////////

void RadiativeTransferFVDOM::computeQ(const CFuint ib, 
				      const CFuint d)
{      
  CFLog(VERBOSE, "RadiativeTransferFVDOM::computeQ() in (bin, dir) = ("
	<< ib << ", " << d << ") => start\n");
  
  DataHandle<CFreal> volumes = socket_volumes.getDataHandle();
  DataHandle<CFreal> divQ = socket_divq.getDataHandle();
  CellTrsGeoBuilder::GeoData& geoData = m_geoBuilder.getDataGE();
  SafePtr<TopologicalRegionSet> cells = geoData.trs;
  const CFuint nbCells = cells->getLocalNbGeoEnts();
  cf_assert(m_advanceOrder[d].size() == nbCells);
  
  // precompute the dot products for all faces and directions (a part from the sign)
  if (m_loopOverBins) {
    computeDotProdInFace(d, m_dotProdInFace);
  }
  
  SafePtr<ConnectivityTable<CFuint> > cellFaces = MeshDataStack::getActive()->getConnectivity("cellFaces");
  DataHandle<CFint> isOutward = socket_isOutward.getDataHandle();
  
  for (CFuint m = 0; m < nbCells; m++) {
    CFreal inDirDotnANeg = 0.;
    CFreal Ic            = 0.;
    
    // allocate the cell entity
    const CFuint iCell = std::abs(m_advanceOrder[d][m]);
    
    // new algorithm (more parallelizable): opacities are computed cell by cell
    // for a given bin
    if (!m_oldAlgo) {getFieldOpacities(ib, iCell);} 
    
    const CFuint nbFaces = cellFaces->nbCols(iCell);
    
    if(m_useExponentialMethod){
      inDirDotnANeg = 0.;
      CFreal dirDotnANeg = 0;
      CFreal Lc      = 0;
      CFreal halfExp = 0;
      
      for (CFuint iFace = 0; iFace < nbFaces; ++iFace) {
	const CFuint faceID = (*cellFaces)(iCell, iFace);
	const CFreal factor = ((CFuint)(isOutward[faceID]) != iCell) ? -1. : 1.;
	const CFreal dirDotNA = m_dotProdInFace[faceID]*factor;
	
	if(dirDotNA < 0.) {
	  dirDotnANeg += dirDotNA;
	  
	  const bool isBFace = m_mapGeoToTrs->isBGeo(faceID);
	  if (!isBFace){
	    const CFuint neighborCellID = getNeighborCellID(faceID, iCell);
	    inDirDotnANeg += m_In[neighborCellID]*dirDotNA;
	  }
	  else {
	    const CFreal boundarySource = m_fieldSource[iCell];
	    inDirDotnANeg += boundarySource*dirDotNA;
	  }
	}
      } 
      Lc          = volumes[iCell]/(- dirDotnANeg); 
      halfExp     = std::exp(-0.5*Lc*m_fieldAbsor[iCell]);
      m_In[iCell] = (inDirDotnANeg/dirDotnANeg)*halfExp*halfExp + (1. - halfExp*halfExp)*m_fieldSource[iCell];
      Ic          = (inDirDotnANeg/dirDotnANeg)*halfExp + (1. - halfExp)*m_fieldSource[iCell];
    }
    else{
      CFreal dirDotnAPos = 0;
      for (CFuint iFace = 0; iFace < nbFaces; ++iFace) {
	const CFuint faceID = (*cellFaces)(iCell, iFace);
	const CFreal factor = ((CFuint)(isOutward[faceID]) != iCell) ? -1. : 1.;
	const CFreal dirDotNA = m_dotProdInFace[faceID]*factor;
	
	if (dirDotNA >= 0.){
	  dirDotnAPos += dirDotNA;
	}
	else {
	  const bool isBFace = m_mapGeoToTrs->isBGeo(faceID);
	  if (!isBFace){
	    const CFuint neighborCellID = getNeighborCellID(faceID, iCell);
	    inDirDotnANeg += m_In[neighborCellID]*dirDotNA;
	  }
	  else {
	    const CFreal boundarySource = m_fieldSource[iCell];
	    inDirDotnANeg += boundarySource*dirDotNA;
	  }
	}
      } 
      m_In[iCell] = (m_fieldAbSrcV[iCell] - inDirDotnANeg)/(m_fieldAbV[iCell] + dirDotnAPos);
      Ic = m_In[iCell];
    }
    
    m_q(iCell,XX) += Ic*m_dirs(d,0)*m_weight[d];
    m_q(iCell,YY) += Ic*m_dirs(d,1)*m_weight[d];
    m_q(iCell,ZZ) += Ic*m_dirs(d,2)*m_weight[d];
    
    CFreal inDirDotnA = inDirDotnANeg;
    for (CFuint iFace = 0; iFace < nbFaces; ++iFace) {
      const CFuint faceID = (*cellFaces)(iCell, iFace);
      const CFreal factor = ((CFuint)(isOutward[faceID]) != iCell) ? -1. : 1.;
      const CFreal dirDotNA = m_dotProdInFace[faceID]*factor;
      if (dirDotNA > 0.) {
	inDirDotnA += m_In[iCell]*dirDotNA;
      }
    }
    
    divQ[iCell] += inDirDotnA*m_weight[d];
    m_II[iCell]   += Ic*m_weight[d];
  }  
  
  CFLog(VERBOSE, "RadiativeTransferFVDOM::computeQ() in (bin, dir) = ("
	<< ib << ", " << d << ") => end\n");
}
      
//////////////////////////////////////////////////////////////////////////////

void RadiativeTransferFVDOM::diagnoseProblem(const CFuint d, 
					     const CFuint m, 
					     const CFuint mLast)
{
  //Check that it wrote a cell in the current stage
  std::cout << "No cell added to advance list in direction number = " << d <<". Problem with mesh.\n";
  
  const CFreal pi = MathTools::MathConsts::CFrealPi();
  
  /// Test to rotate the directions
  CFreal xAngleRotation = 0.;
  CFreal yAngleRotation = 0.;
  CFreal zAngleRotation = 0.;
  std::cout << "Try rotating in x. Introduce an angle in degrees\n";
  std::cout << "Theta = \n";
  std::cin >> xAngleRotation;
  std::cout << "Try rotating in y. Introduce an angle in degrees\n";
  std::cout << "Phi = \n";
  std::cin >> yAngleRotation;
  std::cout << "Try rotating in z. Introduce an angle in degrees\n";
  std::cout << "Psi = \n";
  std::cin >> zAngleRotation;
  
  xAngleRotation *= pi/180.;
  yAngleRotation *= pi/180.;
  zAngleRotation *= pi/180.;
  
  for(CFuint dirs = 0; dirs < m_nbDirs; dirs++){
    //Rotating over x
    const CFreal rot0 = m_dirs(dirs,0);
    const CFreal rot1 = m_dirs(dirs,1)*std::cos(xAngleRotation) - m_dirs(dirs,2)*std::sin(xAngleRotation);
    const CFreal rot2 = m_dirs(dirs,1)*std::sin(xAngleRotation) + m_dirs(dirs,2)*std::cos(xAngleRotation);
    //Rotating over y
    const CFreal rot3 = rot0*std::cos(yAngleRotation) + rot2*std::sin(yAngleRotation);
    const CFreal rot4 = rot1;
    const CFreal rot5 = -rot0*std::sin(yAngleRotation) + rot2*std::cos(yAngleRotation);
    //Rotating over z
    const CFreal rot6 = rot3*std::cos(zAngleRotation) - rot4*std::sin(zAngleRotation);
    const CFreal rot7 = rot3*std::sin(zAngleRotation) + rot4*std::cos(zAngleRotation);
    const CFreal rot8 = rot5;
    
    m_dirs(dirs,0) = rot6;
    m_dirs(dirs,1) = rot7;
    m_dirs(dirs,2) = rot8;
    CFLog(VERBOSE, "dirs[" << dirs <<"] = ("<<  m_dirs(dirs,0) <<", " << m_dirs(dirs,1) <<", "<<m_dirs(dirs,2)<<")\n");
  }
  CFLog(ERROR, "RadiativeTransferFVDOM::getAdvanceOrder() => No cell added to advance list. Problem with mesh\n");
  return; // goto directions_loop;
  cf_assert(m != mLast);
}
      
//////////////////////////////////////////////////////////////////////////////
  
void RadiativeTransferFVDOM::computeDotProdInFace
(const CFuint d, LocalArray<CFreal>::TYPE& dotProdInFace)
{
  const CFuint DIM = 3;
  DataHandle<CFreal> normals = socket_normals.getDataHandle();
  const CFuint totalNbFaces = normals.size()/DIM;
  cf_assert(dotProdInFace.size() == totalNbFaces);
  RealVector normalPtr(DIM, static_cast<CFreal*>(NULL));
  
  for (CFuint faceID = 0; faceID < totalNbFaces; ++faceID) {
    const CFuint startID = faceID*DIM;
    normalPtr.wrap(DIM, &normals[startID]);
    // the sign of each dot product will depend on the actual considered cell within the loop
    dotProdInFace[faceID] = getDirDotNA(d, normalPtr);
  }
}      
      

//////////////////////////////////////////////////////////////////////////////

void RadiativeTransferFVDOM::loopOverBins(const CFuint startBin, 
					  const CFuint endBin, 
					  const CFuint startDir,
					  const CFuint endDir)
{
  for(CFuint ib = startBin; ib < endBin; ++ib) {
    CFLog(INFO, "( bin: " << ib << " ), ( dir: ");
    // old algorithm: opacities computed for all cells at once for a given bin
    if (m_oldAlgo) {getFieldOpacities(ib);}
    const CFuint dStart = (ib != startBin) ? 0 : startDir;
    const CFuint dEnd = (ib != m_startEndBin.second)? m_nbDirs : endDir;
    for (CFuint d = dStart; d < dEnd; ++d) {
      // CFLog(INFO, "[bin, dir] = [" << ib << ", " << d << "]\n");
      CFLog(INFO, d << " ");
      computeQ(ib,d);
    }
    CFLog(INFO, ")\n");
  }
}
      
//////////////////////////////////////////////////////////////////////////////

void RadiativeTransferFVDOM::loopOverDirs(const CFuint startBin, 
					  const CFuint endBin, 
					  const CFuint startDir,
					  const CFuint endDir)
{
  
  
  for (CFuint d = startDir; d < endDir; ++d) {
    CFLog(INFO, "( dir: " << d << " ), ( bin: ");
    const CFuint bStart = (d != startDir) ? 0 : startBin;
    const CFuint bEnd   = (d != m_startEndDir.second) ? m_nbBins : endBin;
    // precompute dot products for all faces and directions (a part from the sign)
    computeDotProdInFace(d, m_dotProdInFace);
    
    for(CFuint ib = startBin; ib < endBin; ++ib) {
      // CFLog(INFO, "[dir, bin] = [" << d << ", " << ib << "]\n");
      // old algorithm: opacities computed for all cells at once for a given bin
      CFLog(INFO, ib << " ");
      if (m_oldAlgo) {getFieldOpacities(ib);}
      computeQ(ib,d);
    }
    CFLog(INFO, ")\n");
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace RadiativeTransfer

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

