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
  socket_faceAreas("faceAreas"),
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
  m_qrAv(),
  m_divqAv(),
  m_nbBins(1),
  m_nbBands(1),
  m_multiSpectralIdx(),
  m_nbDirTypes()
{
  addConfigOptionsTo(this);
  
  m_nbDirs = 8;
  this->setParameter("nDirs",&m_nbDirs);
  
  m_dirName = Environment::DirPaths::getInstance().getWorkingDir();
  setParameter("DirName", &m_dirName);
  
  m_binTabName = "";
  setParameter("BinTabName", &m_binTabName); 
  
  m_outTabName = "tableout.dat";
  setParameter("OutTabName", &m_outTabName);
  
  m_writeToFile = true;
  setParameter("WriteToFile", &m_writeToFile);
    
  m_useExponentialMethod = true;
  setParameter("UseExponentialMethod", &m_useExponentialMethod);
  
  m_radialData = false;
  setParameter("RadialData", &m_radialData);
  
  m_oldAlgo = false;
  setParameter("OldAlgorithm", &m_oldAlgo);
  
  m_binningPARADE = false;
  setParameter("BinningPARADE", &m_binningPARADE);
  
  m_Nr = 50;
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
    
  m_dirGenerator = "Default";
  setParameter("DirectionsGenerator", &m_dirGenerator);

  m_theta_max = 180.;
  setParameter("theta_max", &m_theta_max);

  m_nb_pts_polar = 16;
  setParameter("nb_pts_polar", &m_nb_pts_polar);

  m_nb_pts_azi = 16;
  setParameter("nb_pts_azi", &m_nb_pts_azi);

  m_rule_polar = "GL";
  setParameter("rule_polar", &m_rule_polar);

  m_rule_azi = "TRAP";
  setParameter("rule_azi", &m_rule_azi);
  
  m_directions = false;
  setParameter("directions", &m_directions);

  m_TGSData = false;
  setParameter("TGSData", &m_TGSData);

  m_wallEmissivity = 1;
  setParameter("wallEmissivity", &m_wallEmissivity);

// AL: to be removed once a cleaner solution is found 
  m_radNamespace = "Default";
  setParameter("RadNamespace", &m_radNamespace);
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
  options.addConfigOption< bool >
    ("BinningPARADE",
     "Computes the radiative properties based on the binning done with data from PARADE database");
  options.addConfigOption< CFuint >("NbRadialPoints","Number of Radial points in the sphere case");
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
  options.addConfigOption< string >("DirectionsGenerator","Name of the method for generating directions.");
  options.addConfigOption< CFreal >("theta_max","Maximum value of theta.");
  options.addConfigOption< CFuint >("nb_pts_polar","Number of polar points.");
  options.addConfigOption< CFuint >("nb_pts_azi","Number of azimut points.");
  options.addConfigOption< string >("rule_polar","Rule for polars computation.");
  options.addConfigOption< string >("rule_azi","Rule for azimut computation.");
  options.addConfigOption< bool >("directions","Option to write directions");
  options.addConfigOption< bool >("TGSData","axial q and divQ for the TGS case");
  options.addConfigOption< CFreal >("wallEmissivity","The value of wall emissivity");
  // AL: to be removed once a cleaner solution is found 
  options.addConfigOption< string >
    ("RadNamespace","Namespace grouping all ranks involved in parallel communication");
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
  result.push_back(&socket_faceAreas);
  
  return result;
}

//////////////////////////////////////////////////////////////////////////////

void RadiativeTransferFVDOM::setup()
{
  CFAUTOTRACE;
  
  CFLog(VERBOSE, "RadiativeTransferFVDOM::setup() => START\n");
  
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
  sockets.faceAreas   = socket_faceAreas;
  
  // source sockets cannot be copied
  sockets.alpha_avbin = socket_alpha_avbin.getDataHandle();
  sockets.B_bin       = socket_B_bin.getDataHandle();
 
  m_radiation->setupDataSockets(sockets);
  m_radiation->setup();
  m_radiation->configureTRS();
  m_radiation->setupAxiFlag(false);
  // vector<string> boundaryTrsNames;
  // m_radiation->getBoundaryTRSnames(boundaryTrsNames);
  
  // copy the content of the "initialNodalStates" (if available) into "nstates"
  const string nsp = MeshDataStack::getActive()->getPrimaryNamespace();
  const string dhName = nsp + "_initialNodalStates";
  if (MeshDataStack::getActive()->getDataStorage()->checkData(dhName)) {
    DataHandle<CFreal> initialNodalStates =
      MeshDataStack::getActive()->getDataStorage()->getData<CFreal>(dhName); 
  
    /// storage of nstates (states in nodes)
    DataHandle<RealVector> nstates = socket_nstates.getDataHandle();  
    const CFuint nbNodalVars = nstates[0].size();
    cf_assert(nstates.size()*nbNodalVars == initialNodalStates.size());
    for (CFuint i = 0; i < nstates.size(); ++i) {
      const CFuint startn = i*nbNodalVars;
      for (CFuint j = 0; j < nbNodalVars; ++j) {
	nstates[i][j] = initialNodalStates[startn+j]; 
      }
    }  
    
    MeshDataStack::getActive()->getDataStorage()->deleteData<CFreal>(dhName);
  }
  
  if (m_radiation->hasRadiationPhysics()) {
    SafePtr<Radiator> radiator = 
      m_radiation->getCellDistPtr(0)->getRadiatorPtr();
    if (radiator->getName() != "NullRadiator") {
      m_nbBins  = radiator->getNbBins();
      m_nbBands = radiator->getNbBands();
    }
  }
  
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
  cf_assert(PE::GetPE().GetProcessorCount(nsp) == 1);
  
  // Setting up the file containing the binary table with opacities
  CFLog(VERBOSE, "RadiativeTransferFVDOM::setup() => m_dirName = "<< m_dirName <<"\n"); 
  CFLog(VERBOSE, "RadiativeTransferFVDOM::setup() => m_binTabName = "<< m_binTabName <<"\n");
  
  m_inFileHandle = Environment::SingleBehaviorFactory<Environment::FileHandlerInput>::getInstance().create();
  m_binTableFile = m_dirName / boost::filesystem::path(m_binTabName);
  
  CFLog(VERBOSE, "RadiativeTransferFVDOM::setup() => m_binTableFile = "<< m_binTableFile <<"\n");
  
  // Check to force an allowed number of directions for Default
  if ((m_dirGenerator == "Default") && (m_nbDirs != 8) && (m_nbDirs != 24) && (m_nbDirs != 48) && (m_nbDirs != 80)) {
    CFLog(WARN, "RadiativeTransferFVDOM::setup() => This ndirs is not allowed for Default directions generator. 8 directions is chosen \n");
    m_nbDirs = 8;
  } 
  
  // Check to force an allowed number points for Munafo (GL method)
  if ((m_dirGenerator == "Munafo") && (m_rule_polar=="GL") && (m_nb_pts_polar != 2) && (m_nb_pts_polar != 4) && (m_nb_pts_polar != 8) && (m_nb_pts_polar != 16) && (m_nb_pts_polar != 19) && (m_nb_pts_polar != 31) && (m_nb_pts_polar != 32) && (m_nb_pts_polar != 64) && (m_nb_pts_polar != 96))  {
   CFLog(WARN, "RadiativeTransferFVDOM::setup() => This nb_pts_polar is not allowed for Monafo directions generator (GL). 16 points are chosen \n");
   m_nb_pts_polar = 16;
  }
  
  if (readOpacityTables()) {
    // Reading the table
    readOpacities();
  }
  
  const CFuint DIM = 3;
  
  if(m_dirGenerator == "Munafo"){
    m_nbDirs = m_nb_pts_polar*m_nb_pts_azi;
  }

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
  cf_assert(m_nbBands > 0);

  if(m_nbBins > 1 && m_nbBands == 1){
    m_multiSpectralIdx = m_nbBins;
  }
  else if(m_nbBins ==1 && m_nbBands > 1){
    m_multiSpectralIdx = m_nbBands;
  }
  else if(m_nbBins > 1 && m_nbBands > 1){
    m_multiSpectralIdx = m_nbBands*m_nbBins;
  }
  
  // set the start/end bins for this process
  if (m_nbThreads == 1) { 
    m_startEndDir.first  = 0;
    m_startEndBin.first  = 0;
    m_startEndDir.second = m_nbDirs-1;
    m_startEndBin.second = m_multiSpectralIdx-1;
  }
  else {
    const CFuint nbBinDir     = m_multiSpectralIdx*m_nbDirs; 
    const CFuint minNbThreadsPerProc = nbBinDir/m_nbThreads;
    const CFuint maxNbThreadsPerProc = minNbThreadsPerProc + nbBinDir%m_nbThreads;
    cf_assert(minNbThreadsPerProc > 0);
    
    // same direction has same meshdata structure, therefore if you have 
    // m_nbThreads <= m_nbDirs it is more scalable to split by direction 
    const CFuint startThread = m_threadID*minNbThreadsPerProc;
    const CFuint nbThreadsPerProc = (m_threadID < m_nbThreads-1) ? 
      minNbThreadsPerProc : maxNbThreadsPerProc; 
    const CFuint endThread = startThread + nbThreadsPerProc;
    
    CFLog(INFO, "RadiativeTransferFVDOM::setup() => nbBins, nbDirs, nbBands   = [" << m_nbBins << ", " << m_nbDirs << ", " << m_nbBands <<"]\n");
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
      m_startEndDir.first = startThread/m_multiSpectralIdx; 
      m_startEndBin.first = startThread%m_multiSpectralIdx;
      m_startEndDir.second = (endThread-1)/m_multiSpectralIdx;
      m_startEndBin.second = (endThread-1)%m_multiSpectralIdx;
    }
  }
  
  const CFuint startBin = m_startEndBin.first;
  const CFuint endBin   = m_startEndBin.second+1;
  cf_assert(endBin <= m_multiSpectralIdx);
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
      const CFuint bEnd   = (d != m_startEndDir.second) ? m_multiSpectralIdx : endBin;
      for (CFuint ib = bStart; ib < bEnd; ++ib) {
	CFLog(VERBOSE, "(" << d << ", " << ib <<"), ");
      }
      CFLog(VERBOSE, "\n");
    }
    CFLog(VERBOSE, "\n");
  }
  
  m_dirs.resize(m_nbDirs*3);
  cf_assert(endDir <= m_nbDirs);
  // 1D array (logically 2D) to store advanceOrder
  m_advanceOrder.resize(nbCells*(endDir-startDir));
  cf_assert(m_advanceOrder.size() > 0);
  
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
  
  // this is fundamental since += is used afterwards
  divQ   = 0.0;
  qx     = 0.0;
  qy     = 0.0;
  qz     = 0.0;
  CellID = 0.0;
  TempProfile = 0.0;
  
  // resize the bins storage
  socket_alpha_avbin.getDataHandle().resize(nbCells*m_multiSpectralIdx);
  socket_B_bin.getDataHandle().resize(nbCells*m_multiSpectralIdx);
  
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
  socket_qradFluxWall.getDataHandle() = 0.0;

  //Averages for the Sphere case
  if (m_radialData || m_TGSData) {
    m_qrAv.resize(m_Nr);
    m_divqAv.resize(m_Nr);
    m_qrAv   = 0.;
    m_divqAv = 0.;  
  }
  
  Stopwatch<WallTime> stp;
  
  stp.start();
  getDirections();
  CFLog(INFO, "RadiativeTransferFVDOM::setup() => getDirections() took " << stp.read() << "s\n");
  
  stp.start();
  
  if (!m_emptyRun) {
    // only get advance order for the considered directions
    CFuint countd = 0;
    for (CFuint d = startDir; d < endDir; ++d, ++countd){
      getAdvanceOrder(d, &m_advanceOrder[countd*nbCells]);
    }
  }
    
  CFLog(INFO, "RadiativeTransferFVDOM::setup() => getAdvanceOrder() took " << stp.read() << "s\n");
    
  CFLog(VERBOSE, "RadiativeTransferFVDOM::setup() => END\n");
}
      
//////////////////////////////////////////////////////////////////////////////

void RadiativeTransferFVDOM::getDirections()
{
  CFLog(DEBUG_MIN, "RadiativeTransferFVDOM::getDirections() => start\n");
  
  const CFreal pi = MathTools::MathConsts::CFrealPi();
  
  CFLog(VERBOSE, "RadiativeTransferFVDOM::getDirections() => Number of Directions = " << m_nbDirs << "\n");
  
  RealMatrix mdirs(m_nbDirs, 3, &m_dirs[0]);

  if(m_dirGenerator == "Default") 
  {

  CFLog(VERBOSE, "RadiativeTransferFVDOM::getDirections() => Directions generator is Default \n");
 
  const CFreal overSq3 = 1./std::sqrt(3.);
  switch(m_nbDirs) {
  case 8:
    m_weight[0] = 4.*pi/m_nbDirs;
    mdirs(0,0) = overSq3;
    mdirs(0,1) = overSq3;
    mdirs(0,2) = overSq3;
    break;
  case 24:
    m_weight[0] = 4.*pi/m_nbDirs;
    mdirs(0,0) = 0.2958759;
    mdirs(0,1) = 0.2958759;
    mdirs(0,2) = 0.9082483;
    break;
  case 48:
    m_weight[0] = 0.1609517;
    mdirs(0,0) = 0.1838670;
    mdirs(0,1) = 0.1838670;
    mdirs(0,2) = 0.9656013;
    m_weight[1] = 0.3626469;
    mdirs(1,0) = 0.1838670;
    mdirs(1,1) = 0.6950514;
    mdirs(1,2) = 0.6950514;
    break;
  case 80:
    m_weight[0] = 0.1712359;
    mdirs(0,0) = 0.1422555;
    mdirs(0,1) = 0.1422555;
    mdirs(0,2) = 0.9795543;
    m_weight[1] = 0.0992284;
    mdirs(1,0) = 0.1422555;
    mdirs(1,1) = overSq3;
    mdirs(1,2) = 0.8040087;
    m_weight[2] = 0.4617179;
    mdirs(2,0) = overSq3;
    mdirs(2,1) = overSq3;
    mdirs(2,2) = overSq3;
    break;
  default:	// nDirs = 8
    m_weight[0] = 4.*pi/m_nbDirs;
    mdirs(0,0) = overSq3;
    mdirs(0,1) = overSq3;
    mdirs(0,2) = overSq3;
    break;
  }
  
  CFuint d = m_nbDirTypes - 1; //Note that it has been changed, because the counters start at 0
  for (CFuint dirType = 0; dirType < m_nbDirTypes; dirType++){
    for (CFuint p = 0; p <= 2; p++){
      const CFuint l = p;	    //Note that it's different because the counter starts at 0
      const CFuint m = (p+1) % 3; //Note a % b is the remainder of the division a/b
      const CFuint n = (p+2) % 3;
      
      if (p == 0 || mdirs(dirType,0) != mdirs(dirType,1) ||
	  mdirs(dirType,1) != mdirs(dirType,2) || mdirs(dirType,2) != mdirs(dirType,0)) {
        CFLog(VERBOSE, "Case1::dirTypes = " << dirType <<"\n");
	CFLog(DEBUG_MIN, "l = " << l << "m = " << m << "n = " << n  <<"\n");
	for (int i = 0; i <= 1; i++) {
	  for (int j = 0; j <= 1; j++) {
	    for (int k = 0; k <= 1; k++) {
	      if ( p+i+j+k != 0) {
		//Note that this is different because the counters are different
		d += 1;
		m_weight[d] = m_weight[dirType];
		mdirs(d,0) = std::pow(-1.,i)*mdirs(dirType,l);
		mdirs(d,1) = std::pow(-1.,j)*mdirs(dirType,m);
		mdirs(d,2) = std::pow(-1.,k)*mdirs(dirType,n);
		CFLog(DEBUG_MIN, "l = " << l << " m = " << m << " n = " << n  <<"\n");
		CFLog(DEBUG_MIN, "d = " << d <<"\n");
		CFLog(DEBUG_MIN, "dirs[" << d <<"] = ("<<  mdirs(d,0) <<", " << mdirs(d,1) <<", "<<mdirs(d,2)<<")\n");
	      }
	    }
	  }
	}
      }     
      if (mdirs(dirType,0) != mdirs(dirType,1) && mdirs(dirType,1) != mdirs(dirType,2) 
	  && mdirs(dirType,2) != mdirs(dirType,0)) {
	CFLog(VERBOSE, "Case2::dirTypes = " << dirType <<"\n");
	CFLog(DEBUG_MIN, "l = " << l << "m = " << m << "n = " << n  <<"\n");
	for (int i = 0; i <= 1; i++) {
	  for (int j = 0; j <= 1; j++) {
	    for (int k = 0; k <= 1; k++) {
	      //Note that this is different because the counters are different
	      d += 1;
	      m_weight[d] = m_weight[dirType];
	      mdirs(d,0) = std::pow(-1.,i)*mdirs(dirType,l);
	      mdirs(d,1) = std::pow(-1.,j)*mdirs(dirType,m);
	      mdirs(d,2) = std::pow(-1.,k)*mdirs(dirType,n);
	      CFLog(DEBUG_MIN, "l = " << l << " m = " << m << " n = " << n  <<"\n");
	      CFLog(DEBUG_MIN, "d = " << d <<"\n");
	      CFLog(DEBUG_MIN, "dirs[" << d <<"] = ("<<  mdirs(d,0) <<", " << mdirs(d,1) <<", "<<mdirs(d,2)<<")\n");
	    }
	  }
	}
      }          
    }
  }
  }

  else if(m_dirGenerator == "Munafo")
  {

   CFLog(VERBOSE, "RadiativeTransferFVDOM::getDirections() => Directions generator is Munafo \n"); 

   //Vectors declaration
   std::vector<CFreal> theta_vec;
   std::vector<CFreal> cos_theta_vec;
   std::vector<CFreal> sin_theta_vec;
   std::vector<CFreal> phi_vec;
   std::vector<CFreal> cos_phi_vec;
   std::vector<CFreal> sin_phi_vec;
   std::vector<CFreal> w_theta;
   std::vector<CFreal> w_phi;
   std::vector<CFreal> w_total;
   
   //conversion to rad
   m_theta_max = m_theta_max*pi/180.;

   //Integration limits for polar angle (theta)
   CFreal a_theta = 0.;
   CFreal b_theta = m_theta_max;

   //Integration limits for azimuthal angle (phi)
   CFreal a_phi = 0.;
   CFreal b_phi = 2.*pi;

   //Number of discrete points
   CFreal nb_pts = m_nb_pts_polar * m_nb_pts_azi;

   CFreal cos_theta;
   CFreal sin_theta;
   
   //Vectors allocation
   theta_vec.resize(m_nb_pts_polar+1);
   cos_theta_vec.resize(m_nb_pts_polar+1);
   sin_theta_vec.resize(m_nb_pts_polar+1);
   phi_vec.resize(m_nb_pts_azi+1);
   cos_phi_vec.resize(m_nb_pts_azi+1);
   sin_phi_vec.resize(m_nb_pts_azi+1);
   w_theta.resize(m_nb_pts_polar+1);
   w_phi.resize(m_nb_pts_azi+1);
   w_total.resize(nb_pts);

   // Useless value to switch from fortran to c++ vectors
   theta_vec[0] = 0.;
   cos_theta_vec[0] = 0.;
   sin_theta_vec[0] = 0.;
   phi_vec[0] = 0.;
   cos_phi_vec[0] = 0.;
   sin_phi_vec[0] = 0.;
   w_theta[0] = 0.;
   w_phi[0] = 0.;

   //Polar angle computation
    if(m_rule_polar == "GL") {
      
      CFLog(VERBOSE, "RadiativeTransferFVDOM::getDirections() => The method for polar angle computation is Gauss-Legendre \n");

      a_theta = std::cos(a_theta);
      b_theta = std::cos(b_theta);
      
     switch(m_nb_pts_polar) {

      case 2:
      // Nodes
      theta_vec[2] = 0.5773502691896257645091488;
      // Weights
      w_theta[2] = 1.;
      break;

      case 4:
      // Nodes
      theta_vec[3] = 0.3399810435848562648026658;
      theta_vec[4] = 0.8611363115940525752239465;
      // Weights
      w_theta[3] = 0.6521451548625461426269361;
      w_theta[4] = 0.3478548451374538573730639;
      break;

      case 8:
      // Nodes     
      theta_vec[5] = 0.1834346424956498049394761;
      theta_vec[6] = 0.5255324099163289858177390;
      theta_vec[7] = 0.7966664774136267395915539;
      theta_vec[8] = 0.9602898564975362316835609;
      // Weights
      w_theta[5] = 0.3626837833783619829651504;
      w_theta[6] = 0.3137066458778872873379622;
      w_theta[7] = 0.2223810344533744705443560;
      w_theta[8] = 0.1012285362903762591525314;
      break;

      case 16:
      // Nodes
      theta_vec[9]  = 0.0950125098376374401853193;
      theta_vec[10] = 0.2816035507792589132304605;
      theta_vec[11] = 0.4580167776572273863424194;
      theta_vec[12] = 0.6178762444026437484466718;
      theta_vec[13] = 0.7554044083550030338951012;
      theta_vec[14] = 0.8656312023878317438804679;
      theta_vec[15] = 0.9445750230732325760779884;
      theta_vec[16] = 0.9894009349916499325961542;
      // Weights
      w_theta[9]  = 0.1894506104550684962853967;
      w_theta[10] = 0.1826034150449235888667637;
      w_theta[11] = 0.1691565193950025381893121;
      w_theta[12] = 0.1495959888165767320815017;
      w_theta[13] = 0.1246289712555338720524763;
      w_theta[14] = 0.0951585116824927848099251;
      w_theta[15] = 0.0622535239386478928628438;
      w_theta[16] = 0.0271524594117540948517806;
      break;

      case 19:
      // Nodes
      theta_vec[10] = 0.;
      theta_vec[11] = 0.1603586456402253758680961;
      theta_vec[12] = 0.3165640999636298319901173; 
      theta_vec[13] = 0.4645707413759609457172671;
      theta_vec[14] = 0.6005453046616810234696382;
      theta_vec[15] = 0.7209661773352293786170959;
      theta_vec[16] = 0.8227146565371428249789225;
      theta_vec[17] = 0.9031559036148179016426609;
      theta_vec[18] = 0.9602081521348300308527788;
      theta_vec[19] = 0.9924068438435844031890177;

      // Weights
      w_theta[10] = 0.1610544498487836959791636;
      w_theta[11] = 0.1589688433939543476499564;
      w_theta[12] = 0.1527660420658596667788554;
      w_theta[13] = 0.1426067021736066117757461;
      w_theta[14] = 0.1287539625393362276755158;
      w_theta[15] = 0.1115666455473339947160239;
      w_theta[16] = 0.0914900216224499994644621;
      w_theta[17] = 0.0690445427376412265807083;
      w_theta[18] = 0.0448142267656996003328382;
      w_theta[19] = 0.0194617882297264770363120;
      break;

      case 31:
      // Nodes
      theta_vec[16] = 0.;
      theta_vec[17] = 0.0995553121523415;
      theta_vec[18] = 0.1981211993355706;
      theta_vec[19] = 0.2947180699817016;
      theta_vec[20] = 0.3883859016082329;
      theta_vec[21] = 0.4781937820449025;
      theta_vec[22] = 0.5632491614071493;
      theta_vec[23] = 0.6427067229242603;
      theta_vec[24] = 0.7157767845868532;
      theta_vec[25] = 0.7817331484166249;
      theta_vec[26] = 0.8399203201462674;
      theta_vec[27] = 0.8897600299482711;
      theta_vec[28] = 0.9307569978966481;
      theta_vec[29] = 0.9625039250929497;
      theta_vec[30] = 0.9846859096651525;
      theta_vec[31] = 0.9970874818194770;

      // Weights
      w_theta[16] = 0.0997205447934265;
      w_theta[17] = 0.0992250112266723;
      w_theta[18] = 0.0977433353863287;
      w_theta[19] = 0.0952902429123195;
      w_theta[20] = 0.0918901138936415;
      w_theta[21] = 0.0875767406084779;
      w_theta[22] = 0.0823929917615893;
      w_theta[23] = 0.0763903865987766;
      w_theta[24] = 0.0696285832354104;
      w_theta[25] = 0.0621747865610284;
      w_theta[26] = 0.0541030824249169;
      w_theta[27] = 0.0454937075272011;
      w_theta[28] = 0.0364322739123855;
      w_theta[29] = 0.0270090191849794;
      w_theta[30] = 0.0173186207903106;
      w_theta[31] = 0.0074708315792488;
      break;

      case 32:
      // Nodes
      theta_vec[17] = 0.0483076656877383162348126;
      theta_vec[18] = 0.1444719615827964934851864;
      theta_vec[19] = 0.2392873622521370745446032;
      theta_vec[20] = 0.3318686022821276497799168;
      theta_vec[21] = 0.4213512761306353453641194;
      theta_vec[22] = 0.5068999089322293900237475;
      theta_vec[23] = 0.5877157572407623290407455;
      theta_vec[24] = 0.6630442669302152009751152;
      theta_vec[25] = 0.7321821187402896803874267;
      theta_vec[26] = 0.7944837959679424069630973;
      theta_vec[27] = 0.8493676137325699701336930;
      theta_vec[28] = 0.8963211557660521239653072;
      theta_vec[29] = 0.9349060759377396891709191;
      theta_vec[30] = 0.9647622555875064307738119;
      theta_vec[31] = 0.9856115115452683354001750;
      theta_vec[32] = 0.9972638618494815635449811;

      // Weights
      w_theta[17] = 0.0965400885147278005667648;
      w_theta[18] = 0.0956387200792748594190820;
      w_theta[19] = 0.0938443990808045656391802;
      w_theta[20] = 0.0911738786957638847128686;
      w_theta[21] = 0.0876520930044038111427715;
      w_theta[22] = 0.0833119242269467552221991;
      w_theta[23] = 0.0781938957870703064717409;
      w_theta[24] = 0.0723457941088485062253994;
      w_theta[25] = 0.0658222227763618468376501;
      w_theta[26] = 0.0586840934785355471452836;
      w_theta[27] = 0.0509980592623761761961632;
      w_theta[28] = 0.0428358980222266806568786;
      w_theta[29] = 0.0342738629130214331026877;
      w_theta[30] = 0.0253920653092620594557526;
      w_theta[31] = 0.0162743947309056706051706;
      w_theta[32] = 0.0070186100094700966004071;
      break;

      case 64:
      // Nodes
      theta_vec[33] = 0.0243502926634244325089554;
      theta_vec[34] = 0.0729931217877990394495429;
      theta_vec[35] = 0.1214628192961205544703765;
      theta_vec[36] = 0.1696444204239928180373136;
      theta_vec[37] = 0.2174236437400070841496487;
      theta_vec[38] = 0.2646871622087674163739642;
      theta_vec[39] = 0.3113228719902109561575127;
      theta_vec[40] = 0.3572201583376681159504426;
      theta_vec[41] = 0.4022701579639916036957668;
      theta_vec[42] = 0.4463660172534640879849477;
      theta_vec[43] = 0.4894031457070529574785263;
      theta_vec[44] = 0.5312794640198945456580139;
      theta_vec[45] = 0.5718956462026340342838781;
      theta_vec[46] = 0.6111553551723932502488530;
      theta_vec[47] = 0.6489654712546573398577612;
      theta_vec[48] = 0.6852363130542332425635584;
      theta_vec[49] = 0.7198818501716108268489402;
      theta_vec[50] = 0.7528199072605318966118638;
      theta_vec[51] = 0.7839723589433414076102205;
      theta_vec[52] = 0.8132653151227975597419233;
      theta_vec[53] = 0.8406292962525803627516915;
      theta_vec[54] = 0.8659993981540928197607834;
      theta_vec[55] = 0.8893154459951141058534040;
      theta_vec[56] = 0.9105221370785028057563807;
      theta_vec[57] = 0.9295691721319395758214902;
      theta_vec[58] = 0.9464113748584028160624815;
      theta_vec[59] = 0.9610087996520537189186141;
      theta_vec[60] = 0.9733268277899109637418535;
      theta_vec[61] = 0.9833362538846259569312993;
      theta_vec[62] = 0.9910133714767443207393824;
      theta_vec[63] = 0.9963401167719552793469245;
      theta_vec[64] = 0.9993050417357721394569056;

      // Weights
      w_theta[33] = 0.0486909570091397203833654;
      w_theta[34] = 0.0485754674415034269347991;
      w_theta[35] = 0.0483447622348029571697695;
      w_theta[36] = 0.0479993885964583077281262;
      w_theta[37] = 0.0475401657148303086622822;
      w_theta[38] = 0.0469681828162100173253263;
      w_theta[39] = 0.0462847965813144172959532;
      w_theta[40] = 0.0454916279274181444797710;
      w_theta[41] = 0.0445905581637565630601347;
      w_theta[42] = 0.0435837245293234533768279;
      w_theta[43] = 0.0424735151236535890073398;
      w_theta[44] = 0.0412625632426235286101563;
      w_theta[45] = 0.0399537411327203413866569;
      w_theta[46] = 0.0385501531786156291289625;
      w_theta[47] = 0.0370551285402400460404151;
      w_theta[48] = 0.0354722132568823838106931;
      w_theta[49] = 0.0338051618371416093915655;
      w_theta[50] = 0.0320579283548515535854675;
      w_theta[51] = 0.0302346570724024788679741;
      w_theta[52] = 0.0283396726142594832275113;
      w_theta[53] = 0.0263774697150546586716918;
      w_theta[54] = 0.0243527025687108733381776;
      w_theta[55] = 0.0222701738083832541592983;
      w_theta[56] = 0.0201348231535302093723403;
      w_theta[57] = 0.0179517157756973430850453;
      w_theta[58] = 0.0157260304760247193219660;
      w_theta[59] = 0.0134630478967186425980608;
      w_theta[60] = 0.0111681394601311288185905;
      w_theta[61] = 0.0088467598263639477230309;
      w_theta[62] = 0.0065044579689783628561174;
      w_theta[63] = 0.0041470332605624676352875;
      w_theta[64] = 0.0017832807216964329472961;  
      break;

      case 96:
      // Nodes
      theta_vec[49] = 0.0162767448496029695791346;
      theta_vec[50] = 0.0488129851360497311119582;
      theta_vec[51] = 0.0812974954644255589944713;
      theta_vec[52] = 0.1136958501106659209112081;
      theta_vec[53] = 0.1459737146548969419891073;
      theta_vec[54] = 0.1780968823676186027594026;
      theta_vec[55] = 0.2100313104605672036028472;
      theta_vec[56] = 0.2417431561638400123279319;
      theta_vec[57] = 0.2731988125910491414872722;
      theta_vec[58] = 0.3043649443544963530239298;
      theta_vec[59] = 0.3352085228926254226163256;
      theta_vec[60] = 0.3656968614723136350308956;
      theta_vec[61] = 0.3957976498289086032850002;
      theta_vec[62] = 0.4254789884073005453648192;
      theta_vec[63] = 0.4547094221677430086356761;
      theta_vec[64] = 0.4834579739205963597684056;
      theta_vec[65] = 0.5116941771546676735855097;
      theta_vec[66] = 0.5393881083243574362268026;
      theta_vec[67] = 0.5665104185613971684042502;
      theta_vec[68] = 0.5930323647775720806835558;
      theta_vec[69] = 0.6189258401254685703863693;
      theta_vec[70] = 0.6441634037849671067984124;
      theta_vec[71] = 0.6687183100439161539525572;
      theta_vec[72] = 0.6925645366421715613442458;
      theta_vec[73] = 0.7156768123489676262251441;
      theta_vec[74] = 0.7380306437444001328511657;
      theta_vec[75] = 0.7596023411766474987029704;
      theta_vec[76] = 0.7803690438674332176036045;
      theta_vec[77] = 0.8003087441391408172287961;
      theta_vec[78] = 0.8194003107379316755389996;
      theta_vec[79] = 0.8376235112281871214943028;
      theta_vec[80] = 0.8549590334346014554627870;
      theta_vec[81] = 0.8713885059092965028737748;
      theta_vec[82] = 0.8868945174024204160568774;
      theta_vec[83] = 0.9014606353158523413192327;
      theta_vec[84] = 0.9150714231208980742058845;
      theta_vec[85] = 0.9277124567223086909646905;
      theta_vec[86] = 0.9393703397527552169318574;
      theta_vec[87] = 0.9500327177844376357560989;
      theta_vec[88] = 0.9596882914487425393000680;
      theta_vec[89] = 0.9683268284632642121736594;
      theta_vec[90] = 0.9759391745851364664526010;
      theta_vec[91] = 0.9825172635630146774470458;
      theta_vec[92] = 0.9880541263296237994807628;
      theta_vec[93] = 0.9925439003237626245718923;
      theta_vec[94] = 0.9959818429872092906503991;
      theta_vec[95] = 0.9983643758631816777241494;
      theta_vec[96] = 0.9996895038832307668276901;

      // Weights
      w_theta[49] = 0.0325506144923631662419614;
      w_theta[50] = 0.0325161187138688359872055;
      w_theta[51] = 0.0324471637140642693640128;
      w_theta[52] = 0.0323438225685759284287748;
      w_theta[53] = 0.0322062047940302506686671;
      w_theta[54] = 0.0320344562319926632181390;
      w_theta[55] = 0.0318287588944110065347537;
      w_theta[56] = 0.0315893307707271685580207;
      w_theta[57] = 0.0313164255968613558127843;
      w_theta[58] = 0.0310103325863138374232498;
      w_theta[59] = 0.0306713761236691490142288;
      w_theta[60] = 0.0302999154208275937940888;
      w_theta[61] = 0.0298963441363283859843881;
      w_theta[62] = 0.0294610899581679059704363;
      w_theta[63] = 0.0289946141505552365426788;
      w_theta[64] = 0.0284974110650853856455995;
      w_theta[65] = 0.0279700076168483344398186;
      w_theta[66] = 0.0274129627260292428234211;
      w_theta[67] = 0.0268268667255917621980567;
      w_theta[68] = 0.0262123407356724139134580;
      w_theta[69] = 0.0255700360053493614987972;
      w_theta[70] = 0.0249006332224836102883822;
      w_theta[71] = 0.0242048417923646912822673;
      w_theta[72] = 0.0234833990859262198422359;
      w_theta[73] = 0.0227370696583293740013478;
      w_theta[74] = 0.0219666444387443491947564;
      w_theta[75] = 0.0211729398921912989876739;
      w_theta[76] = 0.0203567971543333245952452;
      w_theta[77] = 0.0195190811401450224100852;
      w_theta[78] = 0.0186606796274114673851568;
      w_theta[79] = 0.0177825023160452608376142;
      w_theta[80] = 0.0168854798642451724504775;
      w_theta[81] = 0.0159705629025622913806165;
      w_theta[82] = 0.0150387210269949380058763;
      w_theta[83] = 0.0140909417723148609158616;
      w_theta[84] = 0.0131282295669615726370637;
      w_theta[85] = 0.0121516046710883196351814;
      w_theta[86] = 0.0111621020998384985912133;
      w_theta[87] = 0.0101607705350084157575876;
      w_theta[88] = 0.0091486712307833866325846;
      w_theta[89] = 0.0081268769256987592173824;
      w_theta[90] = 0.0070964707911538652691442;
      w_theta[91] = 0.0060585455042359616833167;
      w_theta[92] = 0.0050142027429275176924702;
      w_theta[93] = 0.0039645543384446866737334;
      w_theta[94] = 0.0029107318179349464084106;
      w_theta[95] = 0.0018539607889469217323359;
      w_theta[96] = 0.0007967920655520124294381;
      break;
     }

     // Apply reflection

     if(m_nb_pts_polar % 2 == 0) {
       for(CFuint i=1; i <= m_nb_pts_polar/2; i++) {
          theta_vec[m_nb_pts_polar/2 - i + 1] = -theta_vec[m_nb_pts_polar/2 +i];
          w_theta[m_nb_pts_polar/2 - i + 1] = w_theta[m_nb_pts_polar/2 +i];
       }
      }

     else {
       CFuint m = m_nb_pts_polar/2 + m_nb_pts_polar % 2;
	 for(CFuint ii=1; ii <= (m_nb_pts_polar/2); ii++) {
          theta_vec[m-ii] = -theta_vec[m+ii];
          w_theta[m-ii] = w_theta[m+ii];
       }
      }
   
     for(CFuint iii=1; iii <= m_nb_pts_polar; iii++) {
       theta_vec[iii] = 0.5 * ((a_theta - b_theta) * theta_vec[iii] + (b_theta + a_theta));
       w_theta[iii] = w_theta[iii]*0.5*(a_theta-b_theta);
       cos_theta_vec[iii] = theta_vec[iii];
       sin_theta_vec[iii] = std::sqrt(1.-std::pow(cos_theta_vec[iii],2));
      }
    } // end of GL


  else  if(m_rule_polar == "TRAP") {
    
    CFLog(VERBOSE, "RadiativeTransferFVDOM::getDirections() => The method for polar angle computation is Trapezoidal \n");

    // weights  
    w_theta[1] = 0.5;    
    for(CFuint i=2; i <= (m_nb_pts_polar - 1); i++)  {
      w_theta[i] = 1.;
    }
    w_theta[m_nb_pts_polar] = 0.5;

    // Nodes
    CFreal dx = (b_theta - a_theta)/(m_nb_pts_polar -1);
    theta_vec[1] = a_theta;
    for(CFuint ii = 2; ii <= m_nb_pts_polar; ii++) {
      theta_vec[ii] = theta_vec[ii - 1] + dx;
    }

    for(CFuint iv=1; iv <= m_nb_pts_polar; iv++) {

    w_theta[iv] = w_theta[iv]*dx;
    cos_theta_vec[iv] = std::cos(theta_vec[iv]);
    sin_theta_vec[iv] = std::sin(theta_vec[iv]);
    }
  } // end of TRAP

  //Azimutal angle computation (TRAP rule)

    CFLog(VERBOSE, "RadiativeTransferFVDOM::getDirections() => The method for azimutal angle computation is Trapezoidal \n");

    // Weights
    w_phi[1] = 0.5;    
    for(CFuint i=2; i <= (m_nb_pts_azi - 1); i++)  {
      w_phi[i] = 1.;
    }
    w_phi[m_nb_pts_azi] = 0.5;

    // Nodes
    CFreal dx = (b_phi - a_phi)/(m_nb_pts_azi - 1);
    phi_vec[1] = a_phi;
    for(CFuint ii = 2; ii <= m_nb_pts_azi; ii++) {
      phi_vec[ii] = phi_vec[ii - 1] + dx;
    }

    for(CFuint i=1; i <= m_nb_pts_azi; i++) {
    w_phi[i] = w_phi[i]*dx;
    cos_phi_vec[i] = std::cos(phi_vec[i]);
    sin_phi_vec[i] = std::sin(phi_vec[i]);
    }

   // Compute direction cosines and weights
    CFuint d = 0;
    for(CFuint p = 1; p <= m_nb_pts_polar; p++) {
      cos_theta = cos_theta_vec[p];
      sin_theta = sin_theta_vec[p];
  
      for(CFuint q = 1; q <= m_nb_pts_azi; q++) {
        //cosines
	mdirs(d,0)=cos_theta;
        mdirs(d,1)=sin_theta*cos_phi_vec[q];
	mdirs(d,2)=sin_theta*sin_phi_vec[q];
 
        //weights
	if(m_rule_polar == "GL") {
          m_weight[d]=w_theta[p]*w_phi[q];
        }

        else {
	  m_weight[d]=w_theta[p]*w_phi[q]*sin_theta;
        }  
        
        ++d; 
      } 
    } 

    m_nbDirTypes = d; 
 }
  // Printing the Directions for debugging
  for (CFuint dir = 0; dir < m_nbDirTypes; dir++) {
    CFLog(DEBUG_MIN, "Direction[" << dir <<"] = (" << mdirs(dir,0) <<", " << mdirs(dir,1) <<", " << mdirs(dir,2) <<")\n");
  }
  
  if (m_directions) {
    writeDirections();
  }
  
  CFLog(DEBUG_MIN, "RadiativeTransferFVDOM::getDirections() => end\n");
}

//////////////////////////////////////////////////////////////////////////////
      
void RadiativeTransferFVDOM::execute()
{
  CFAUTOTRACE;
  
  CFLog(INFO, "RadiativeTransferFVDOM::execute() => START\n");
  
  Stopwatch<WallTime> stp;
  stp.start();
  
  // compute the spectra in steps (to avoid to have to store too much data at once): 
  // this calls Radiator::setupSpectra(wavMin, wavMax)
  const CFuint nbLoops = m_radiation->getNumberLoops();
  for(CFuint i= 0; i < nbLoops; ++i) {
    m_radiation->setupWavStride(i);
  }
  
  CFLog(INFO, "RadiativeTransferFVDOM::execute() => radiation library took " << stp.read() << "s\n");
   
  DataHandle<CFreal> alpha_avbin    = socket_alpha_avbin.getDataHandle();
  DataHandle<CFreal> B_bin          = socket_B_bin.getDataHandle();
  DataHandle<State*, GLOBAL> states = socket_states.getDataHandle();
  
  const CFuint nbCells = states.size();

  for(CFuint i=0;i< nbCells;i++){
    for(CFuint j=1;j< m_multiSpectralIdx;j++){
      CFLog(VERBOSE,"Vector alpha(" << j << "," << i << ") = " << alpha_avbin[j+m_multiSpectralIdx*i] << "\n");
    }
  }

  Common::SafePtr<Common::ConnectivityTable<CFuint> > cells =
    MeshDataStack::getActive()->getConnectivity("cellStates_InnerCells");

  stp.start();
  
  if (!m_emptyRun) {
    DataHandle<CFreal> divQ = socket_divq.getDataHandle();
    
    // only one CPU allow for namespace => the mesh has not been partitioned
    cf_assert(PE::GetPE().GetProcessorCount(getMethodData().getNamespace()) == 1);
    
    // Compute the order of advance
    // Call the function to get the directions
    divQ = 0.0;
    for (CFuint i = 0; i < m_II.size(); ++i) {m_II[i] = 0.0;}
    
    const CFuint startBin = m_startEndBin.first;
    const CFuint endBin   = m_startEndBin.second+1;
    cf_assert(endBin <= m_multiSpectralIdx);
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
    }
    
    if (m_radialData){
      writeRadialData();
    } 

    if (m_TGSData){
      writeTGSData();
    } 

    reduceHeatFlux();
  }
  
  CFLog(INFO, "RadiativeTransferFVDOM::execute() => took " << stp.read() << "s \n");
}
    
//////////////////////////////////////////////////////////////////////////////

void RadiativeTransferFVDOM::getFieldOpacitiesBinning(const CFuint ib)
{
  DataHandle<CFreal> alpha_avbin    = socket_alpha_avbin.getDataHandle();
  DataHandle<CFreal> B_bin          = socket_B_bin.getDataHandle();
  DataHandle<State*, GLOBAL> states = socket_states.getDataHandle();
  
  const CFuint nbCells = states.size();

  for(CFuint i=0;i< nbCells;i++){
    for(CFuint j=0;j< (m_multiSpectralIdx);j++){
      CFLog(VERBOSE,"Vector alpha(" << j << "," << i << ") = " << alpha_avbin[j+m_multiSpectralIdx*i] << "\n");
    }
  }
  
  for (CFuint iCell = 0; iCell < nbCells; iCell++) {
    const State *currState = states[iCell];
    m_fieldSource[iCell] = 0.;
    
    if(m_useExponentialMethod){
      m_fieldAbsor[iCell] = 0.;
    }
    else{
      m_fieldAbSrcV[iCell] = 0.;
      m_fieldAbV[iCell]    = 0.;
    }
    
    DataHandle<CFreal> volumes     = socket_volumes.getDataHandle();
    
    if(m_useExponentialMethod){
      if (B_bin[ib+m_multiSpectralIdx *iCell] <= 1e-30 || alpha_avbin[ib+m_multiSpectralIdx *iCell] <= 1e-30 ){
	m_fieldSource[iCell] = 1e-30;
	m_fieldAbsor[iCell]  = 1e-30;
      }
      else {
	m_fieldSource[iCell] = B_bin[ib+m_multiSpectralIdx *iCell];
	m_fieldAbsor[iCell]  = alpha_avbin[ib+m_multiSpectralIdx *iCell];
      }
    }
    else {
      if (B_bin[ib+m_multiSpectralIdx *iCell] <= 1e-30 || alpha_avbin[ib+m_multiSpectralIdx *iCell] <= 1e-30){
	m_fieldSource[iCell] = 1e-30;
	m_fieldAbV[iCell]    = 1e-30*volumes[iCell]; // Volume converted from m^3 into cm^3
      }
      else {
	m_fieldSource[iCell] = B_bin[ib+m_multiSpectralIdx *iCell];
	m_fieldAbV[iCell]    = alpha_avbin[ib+m_multiSpectralIdx*iCell]*volumes[iCell];
      }      
      m_fieldAbSrcV[iCell]   = m_fieldSource[iCell]*m_fieldAbV[iCell];
    }
  }
}

//////////////////////////////////////////////////////////////////////  
    
void RadiativeTransferFVDOM::getAdvanceOrder(const CFuint d, 
					     CFint *const advanceOrder)
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
	
	CFLog(DEBUG_MAX, "advanceOrder[" << d << "][" << m <<"] = " << iCell << "\n");
	
	advanceOrder[m] = iCell;
	CellID[iCell] = stage;
	m_cdoneIdx.push_back(iCell);
	m += 1;
      }// end if(Cell is not done)
      
    cell_loop:
      const bool dummy = true;
    }// end of the loop over the CELLS
    
    const string msg = "advanceOrder[" + StringOps::to_str(d) + "] = ";
    CFLog(DEBUG_MAX, msg << "\n");
    for (CFuint a = 0; a < nbCells; ++a) {
      CFLog(DEBUG_MAX, advanceOrder[a] << " ");
    }
    CFLog(DEBUG_MAX, "\n");

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
  
  CFLog(DEBUG_MAX, msg << "\n");
  for (CFuint a = 0; a < nbCells; ++a) {
    CFLog(DEBUG_MAX, advanceOrder[a] << " ");
  }
  CFLog(DEBUG_MAX, "\n");
    
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
  
  if (!m_binningPARADE) {m_nbBins = ((int) data[0]);}
  
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
  DataHandle<CFreal> qx = socket_qx.getDataHandle();
  DataHandle<CFreal> qy = socket_qy.getDataHandle();
  DataHandle<CFreal> qz = socket_qz.getDataHandle();
  
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
	m_qrAv[ir]   += (qx[iCell]*x + qy[iCell]*y + qz[iCell]*z)/rCell; //*rCell*rCell; Multiply by r**2 for area-weighted average
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
  CFLog(DEBUG_MIN, "RadiativeTransferFVDOM::getFieldOpacities(" << ib << ", " 
	<< iCell << ") => START\n");
  
  m_fieldSource[iCell] = 0.;
  if(m_useExponentialMethod){
    m_fieldAbsor[iCell] = 0.;
  }
  else{
    m_fieldAbSrcV[iCell] = 0.;
    m_fieldAbV[iCell]    = 0.;
  }
  
  DataHandle<CFreal> TempProfile    = socket_TempProfile.getDataHandle();
  DataHandle<CFreal> volumes        = socket_volumes.getDataHandle();
  DataHandle<CFreal> alpha_avbin    = socket_alpha_avbin.getDataHandle();
  DataHandle<CFreal> B_bin          = socket_B_bin.getDataHandle();
  DataHandle<State*, GLOBAL> states = socket_states.getDataHandle();
  
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
  
  if (!m_binningPARADE) {
    RadiativeTransferFVDOM::DeviceFunc interp;
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
  else {
    if(m_useExponentialMethod){
      if (((alpha_avbin[ib+m_multiSpectralIdx*iCell]) <= 1e-30) || ((B_bin[ib+m_multiSpectralIdx*iCell]) <= 1e-30) ){
	m_fieldSource[iCell] = 1e-30;
	m_fieldAbsor[iCell]  = 1e-30;
      }
      else {
	m_fieldSource[iCell] = B_bin[ib+m_multiSpectralIdx*iCell];
	m_fieldAbsor[iCell]  = alpha_avbin[ib+m_multiSpectralIdx*iCell];
      }
    }
    else{
      if (alpha_avbin[ib+m_multiSpectralIdx*iCell] <= 1e-30 || B_bin[ib+m_multiSpectralIdx*iCell] <= 1e-30 ){
	m_fieldSource[iCell] = 1e-30;
	m_fieldAbV[iCell]    = 1e-30*volumes[iCell]; // Volume converted from m^3 into cm^3
      }
      else {
	m_fieldSource[iCell] = B_bin[ib+m_multiSpectralIdx*iCell];
	m_fieldAbV[iCell]    = alpha_avbin[ib+m_multiSpectralIdx*iCell]*volumes[iCell];
      }      
      m_fieldAbSrcV[iCell]   = m_fieldSource[iCell]*m_fieldAbV[iCell];
    }
  }
  
  CFLog(DEBUG_MIN, "RadiativeTransferFVDOM::getFieldOpacities(" << ib << ", " << iCell << ") => END\n");
} 
    
//////////////////////////////////////////////////////////////////////////////

void RadiativeTransferFVDOM::computeQExponential(const CFuint ib,
						 const CFuint dStart,
						 const CFuint d)
{      
  CFLog(VERBOSE, 
	"RadiativeTransferFVDOM::computeQExponential() in (bin, dir) = ("
	<< ib << ", " << d << ") => start\n");
  DataHandle<CFreal> qradFluxWall = socket_qradFluxWall.getDataHandle();
  DataHandle<CFreal> faceAreas = socket_faceAreas.getDataHandle();
  DataHandle<CFreal> volumes = socket_volumes.getDataHandle();
  DataHandle<CFreal> divQ = socket_divq.getDataHandle();
  CellTrsGeoBuilder::GeoData& geoData = m_geoBuilder.getDataGE();
  SafePtr<TopologicalRegionSet> cells = geoData.trs;
  const CFuint nbCells = cells->getLocalNbGeoEnts();
  
  SafePtr<ConnectivityTable<CFuint> > cellFaces = 
    MeshDataStack::getActive()->getConnectivity("cellFaces");
  DataHandle<CFint> isOutward = socket_isOutward.getDataHandle();
  DataHandle<CFreal> qx = socket_qx.getDataHandle();
  DataHandle<CFreal> qy = socket_qy.getDataHandle();
  DataHandle<CFreal> qz = socket_qz.getDataHandle();
  
  //////// 
  /*const CFuint totalNbFaces = MeshDataStack::getActive()->Statistics().getNbFaces();
  vector<CFint> faceCell(totalNbFaces*2, -1);
  vector<CFuint> nbFacesInCell(nbCells);
  
  for (CFuint iCell = 0; iCell < nbCells; ++iCell) {
    const CFuint nbFaces = cellFaces->nbCols(iCell);
    nbFacesInCell[iCell] = nbFaces;
    for (CFuint iFace = 0; iFace < nbFaces; ++iFace) {
      const CFuint faceID2 = (*cellFaces)(iCell, iFace)*2;
      if (faceCell[faceID2] == -1) {
	faceCell[faceID2] = iCell;
      }
      else {
	faceCell[faceID2+1] = iCell;
      }
    }
    }*/
  //////
  
  // AL: allocation arrays for heat flux (use reserve!!!)
  std::vector< CFuint > wallfIdx;
  std::vector< CFreal > ddd;
  std::vector< CFreal > Ibq;

  const CFuint startCell = (d-dStart)*nbCells;
  for (CFuint m = 0; m < nbCells; m++) {
    CFreal inDirDotnANeg = 0.;
    CFreal Ic            = 0.;
    CFreal dirDotnANeg   = 0.;
    CFreal Lc            = 0.;
    CFreal halfExp       = 0.;
    CFreal POP_dirDotNA  = 0.;
      
    // allocate the cell entity
    cf_assert(startCell+m < m_advanceOrder.size());
    const CFuint iCell = std::abs(m_advanceOrder[startCell+m]);
    
    // new algorithm (more parallelizable): opacities are computed cell by cell
    // for a given bin
    if (!m_oldAlgo) {getFieldOpacities(ib, iCell);} 
    
    const CFuint nbFaces = cellFaces->nbCols(iCell);
    //    cf_assert(nbFaces == nbFacesInCell[iCell]);
    
    for (CFuint iFace = 0; iFace < nbFaces; ++iFace) {
      const CFuint faceID = (*cellFaces)(iCell, iFace);
      const CFreal factor = ((CFuint)(isOutward[faceID]) != iCell) ? -1. : 1.;
      const CFreal dirDotNA = m_dotProdInFace[faceID]*factor;
      
      // AL: check on wall faces and emissivities (use Radiator for this!!!)
      const CFint isWallFace = getWallFaceID(faceID);
      // isInner would always be true...
      const bool isComputingWallFace = qradFluxWall.size() > 0 && isWallFace != -1 && isInner(iCell, faceID) == 1;
      if(isComputingWallFace) {
	wallfIdx.push_back(isWallFace);
	ddd.push_back(dirDotNA*m_weight[d]/faceAreas[faceID]);
	Ibq.push_back(getFaceIbq(faceID));
      }

      if(dirDotNA < 0.) {
	dirDotnANeg += dirDotNA;
	
	/*const CFint fcellID = faceCell[faceID*2]; 
	  const CFint neighborCellID = (fcellID == iCell) ? faceCell[faceID*2+1] : fcellID;
	  const CFreal source = (neighborCellID >=0) ? m_In[neighborCellID] : m_fieldSource[iCell];
	  inDirDotnANeg += source*dirDotNA;*/
	
	const bool isBFace = m_mapGeoToTrs->isBGeo(faceID);
	if (!isBFace){
	  const CFuint neighborCellID = getNeighborCellID(faceID, iCell);
	  inDirDotnANeg += m_In[neighborCellID]*dirDotNA;
	}
	else { // it recognizes a wall as it was a boundary
         if(isComputingWallFace) { // is a wall
	   inDirDotnANeg += m_wallEmissivity*getFaceIbq(faceID)*dirDotNA/m_multiSpectralIdx; //divided by m_multiSpectralIdx
             //CFLog(INFO,"WALL FACE ------------> " << getWallFaceID(faceID) << ", SOURCE ----> " << m_wallEmissivity*getFaceIbq(faceID) <<"\n");
          }
          else { // is another boundary
	  const CFreal boundarySource = m_fieldSource[iCell];
	  inDirDotnANeg += boundarySource*dirDotNA;
	 }
        }
      }
      else if (dirDotNA > 0.) {
	POP_dirDotNA += dirDotNA;
      }
    } 
    
    Lc          = volumes[iCell]/(- dirDotnANeg); 
    halfExp     = std::exp(-0.5*Lc*m_fieldAbsor[iCell]);
    const CFreal InCell = (inDirDotnANeg/dirDotnANeg)*halfExp*halfExp + (1. - halfExp*halfExp)*m_fieldSource[iCell];
    Ic          = (inDirDotnANeg/dirDotnANeg)*halfExp + (1. - halfExp)*m_fieldSource[iCell];
    
    // AL: computation of heat fluxes
    if(wallfIdx.size() > 0) {
      for(CFuint fcount=0; fcount < wallfIdx.size(); fcount++){
	const CFuint IDX = wallfIdx[fcount];
	const CFreal DDD = ddd[fcount];
	const CFreal IBQ = Ibq[fcount];
	if(DDD > 0.0){
          qradFluxWall[IDX] += -m_wallEmissivity*InCell*DDD;
	}
	else{
	  if(ib == 0){ //AL: why for the first bin you do this????
            qradFluxWall[IDX] += m_wallEmissivity*IBQ*std::abs(DDD);
	    //CFLog(INFO,"IBQ = " << IBQ <<"\n");
	  }
	}
      }
      wallfIdx.clear();
      ddd.clear();
      Ibq.clear();
    }

    CFreal inDirDotnA = inDirDotnANeg;
    inDirDotnA += InCell*POP_dirDotNA;
    m_In[iCell] = InCell;
    const CFreal IcWeight = Ic*m_weight[d];
    const CFuint d3 = d*3;
    
    qx[iCell]   += m_dirs[d3]*IcWeight;
    qy[iCell]   += m_dirs[d3+1]*IcWeight;
    qz[iCell]   += m_dirs[d3+2]*IcWeight;
    divQ[iCell] += inDirDotnA*m_weight[d];
    // m_II[iCell] += Ic*m_weight[d]; // useless
    
    /*if (iCell==100 && d == 0) {
      printf ("IcWeight    : %6.6f \n", IcWeight);
      printf ("inDirDotnA  : %6.6f \n",inDirDotnA);
      printf ("InCell      : %6.6f \n", InCell);
      printf ("cellIDin    : %d  \n", iCell*m_nbDirs+d);
      const CFreal qxIcell = qx[iCell];
      printf ("qx[iCell]   : %6.6f  \n", qxIcell);
      const CFreal divqIcell = divQ[iCell];
      printf ("divq[iCell] : %6.6f  \n", divqIcell);
      const CFreal In0 = m_In[iCell];
      printf ("In[iCell]   : %6.6f  \n", In0);
      printf ("d3          : %d  \n", d3);
      printf ("mdirs[d3]   : %6.6f  \n", m_dirs[d3]);
      printf ("mdirs[d3+1] : %6.6f  \n", m_dirs[d3+1]);
        printf ("mdirs[d3+2] : %6.6f  \n", m_dirs[d3+2]);
	exit(1);
	}*/
  }
  
  CFLog(VERBOSE, 
	"RadiativeTransferFVDOM::computeQExponential() in (bin, dir) = ("
	<< ib << ", " << d << ") => end\n");
}
    
//////////////////////////////////////////////////////////////////////////////

void RadiativeTransferFVDOM::computeQNoExponential(const CFuint ib, 
						   const CFuint dStart,
						   const CFuint d)
{      
  CFLog(VERBOSE, 
	"RadiativeTransferFVDOM::computeQNoExponential() in (bin, dir) = ("
	<< ib << ", " << d << ") => start\n");
  
  DataHandle<CFreal> volumes = socket_volumes.getDataHandle();
  DataHandle<CFreal> divQ = socket_divq.getDataHandle();
  CellTrsGeoBuilder::GeoData& geoData = m_geoBuilder.getDataGE();
  SafePtr<TopologicalRegionSet> cells = geoData.trs;
  const CFuint nbCells = cells->getLocalNbGeoEnts();
  
  SafePtr<ConnectivityTable<CFuint> > cellFaces = 
    MeshDataStack::getActive()->getConnectivity("cellFaces");
  DataHandle<CFint> isOutward = socket_isOutward.getDataHandle();
  DataHandle<CFreal> qx = socket_qx.getDataHandle();
  DataHandle<CFreal> qy = socket_qy.getDataHandle();
  DataHandle<CFreal> qz = socket_qz.getDataHandle();
  
  const CFuint startCell = (d-dStart)*nbCells;
  for (CFuint m = 0; m < nbCells; m++) {
    CFreal inDirDotnANeg = 0.;
    CFreal Ic            = 0.;
    CFreal dirDotnAPos   = 0.;
    
    // allocate the cell entity
    cf_assert(startCell+m < m_advanceOrder.size());
    const CFuint iCell = std::abs(m_advanceOrder[startCell+m]);
    
    // new algorithm (more parallelizable): opacities are computed cell by cell
    // for a given bin
    if (!m_oldAlgo) {getFieldOpacities(ib, iCell);} 
    
    const CFuint nbFaces = cellFaces->nbCols(iCell);
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
    
    qx[iCell] += Ic*m_dirs[d*3]*m_weight[d];
    qy[iCell] += Ic*m_dirs[d*3+1]*m_weight[d];
    qz[iCell] += Ic*m_dirs[d*3+2]*m_weight[d];
    
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
    m_II[iCell] += Ic*m_weight[d];
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
  CFLog(ERROR, "RadiativeTransferFVDOM::diagnoseProblem() => No cell added to advance list in direction [" << d << "].\n");
  
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
  
  RealMatrix mdirs(m_nbDirs, 3, &m_dirs[0]);
  
  for(CFuint dirs = 0; dirs < m_nbDirs; dirs++){
    //Rotating over x
    const CFreal rot0 = mdirs(dirs,0);
    const CFreal rot1 = mdirs(dirs,1)*std::cos(xAngleRotation) - mdirs(dirs,2)*std::sin(xAngleRotation);
    const CFreal rot2 = mdirs(dirs,1)*std::sin(xAngleRotation) + mdirs(dirs,2)*std::cos(xAngleRotation);
    //Rotating over y
    const CFreal rot3 = rot0*std::cos(yAngleRotation) + rot2*std::sin(yAngleRotation);
    const CFreal rot4 = rot1;
    const CFreal rot5 = -rot0*std::sin(yAngleRotation) + rot2*std::cos(yAngleRotation);
    //Rotating over z
    const CFreal rot6 = rot3*std::cos(zAngleRotation) - rot4*std::sin(zAngleRotation);
    const CFreal rot7 = rot3*std::sin(zAngleRotation) + rot4*std::cos(zAngleRotation);
    const CFreal rot8 = rot5;
    
    mdirs(dirs,0) = rot6;
    mdirs(dirs,1) = rot7;
    mdirs(dirs,2) = rot8;
    CFLog(VERBOSE, "dirs[" << dirs <<"] = ("<<  mdirs(dirs,0) <<", " << mdirs(dirs,1) <<", "<<mdirs(dirs,2)<<")\n");
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
  
  RadiativeTransferFVDOM::DeviceFunc fun;
  for (CFuint faceID = 0; faceID < totalNbFaces; ++faceID) {
    const CFuint startID = faceID*DIM;
    // the sign of each dot product will depend on the actual considered cell within the loop
    dotProdInFace[faceID] = fun.getDirDotNA(d, &m_dirs[0], &normals[startID]);
  }
}      

//////////////////////////////////////////////////////////////////////////////

void RadiativeTransferFVDOM::writeDirections()
{
  CFLog(VERBOSE, "RadiativeTransferFVDOM::writeDirections() = > Writing directions => start\n");
  
  boost::filesystem::path file = m_dirName / boost::filesystem::path("directions.plt");
  file = PathAppender::getInstance().appendParallel( file );
  
  SelfRegistPtr<Environment::FileHandlerOutput> fhandle = 
    Environment::SingleBehaviorFactory<Environment::FileHandlerOutput>::getInstance().create();
  ofstream& outputFile = fhandle->open(file);
  
  outputFile << "TITLE  = RadiativeTransferFVDOM directions\n";
  outputFile << "VARIABLES = x y z l m n weight dir_ID\n";

  CFreal mu_1 = 0.0;
  RealVector mu_2;
  mu_2.resize(3);
  mu_2[0] = 0.0;
  mu_2[1] = 0.0;
  mu_2[2] = 0.0;
  for(CFuint c=0; c < m_nbDirs; c++) {
    outputFile << "0.0 0.0 0.0 " << m_dirs[c*3+0] << " " << m_dirs[c*3+1] << " " 
	       << m_dirs[c*3+2] << " " << m_weight[c] << " " << c << "\n";
    mu_1 += m_weight[c]; //zeroth moment
    mu_2[0] += m_weight[c]*m_dirs[c*3+0]; //first moment l
    mu_2[1] += m_weight[c]*m_dirs[c*3+1]; //first moment m
    mu_2[2] += m_weight[c]*m_dirs[c*3+2]; //first moment n
  }
  
  fhandle->close();

  CFLog(INFO, "RadiativeTransferFVDOM::writeDirections() = > Zeroth moment = " << mu_1 << "\n");
  CFLog(INFO, "RadiativeTransferFVDOM::writeDirections() = > First moment = ( " << mu_2[0] << ", " << mu_2[1] << ", " << mu_2[2] << " )\n");

  CFLog(VERBOSE, "RadiativeTransferFVDOM::writeDirections() = > Writing directions => end\n");
}      

//////////////////////////////////////////////////////////////////////////////

void RadiativeTransferFVDOM::loopOverBins(const CFuint startBin, 
					  const CFuint endBin, 
					  const CFuint startDir,
					  const CFuint endDir)
{
  CFLog(VERBOSE, "RadiativeTransferFVDOM::loopOverBins() => START\n");
  for(CFuint ib = startBin; ib < endBin; ++ib) {
    CFLog(INFO, "( bin: " << ib << " ), ( dir: ");
    // old algorithm: opacities computed for all cells at once for a given bin
    if (m_oldAlgo && !m_binningPARADE) {getFieldOpacities(ib);}
    if(m_binningPARADE){getFieldOpacitiesBinning(ib);}
    
    /*for (CFuint k = 0; k < m_fieldSource.size(); ++k)
      CFLog(INFO, "m_fieldSource[" <<k << "] => (" << m_fieldSource[k] << "\n");
    for (CFuint k = 0; k < m_fieldAbsor.size(); ++k)
      CFLog(INFO, "m_fieldAbsor[" <<k << "] => (" << m_fieldAbsor[k] << "\n");
    for (CFuint k = 0; k < m_fieldAbSrcV.size(); ++k)
      CFLog(INFO, "m_fieldAbSrcV[" <<k << "] => (" << m_fieldAbSrcV[k] << "\n");
    for (CFuint k = 0; k < m_fieldAbV.size(); ++k)
      CFLog(INFO, "m_fieldAbV[" <<k << "] => (" << m_fieldAbV[k] << "\n");
      exit(1);*/
    
    const CFuint dStart = (ib != startBin) ? 0 : startDir;
    const CFuint dEnd = (ib != m_startEndBin.second)? m_nbDirs : endDir;
    for (CFuint d = dStart; d < dEnd; ++d) {
      CFLog(INFO, d << " ");
      // precompute dot products for all faces and directions (a part from the sign)
      computeDotProdInFace(d, m_dotProdInFace);
      
      (m_useExponentialMethod) ? 
	computeQExponential(ib,dStart,d) : computeQNoExponential(ib,dStart,d);
    }
    CFLog(INFO, ")\n");
  }

  /*  for (CFuint k = 0; k < socket_divq.getDataHandle().size(); ++k) {
    CFLog(INFO, "divQ[" <<k << "] => (" << socket_divq.getDataHandle()[k] << "\n");
  }
  for (CFuint k = 0; k < socket_qx.getDataHandle().size(); ++k) {
    CFLog(INFO, "qx[" <<k << "] => (" << socket_qx.getDataHandle()[k] << "\n");
  }
  for (CFuint k = 0; k < socket_qy.getDataHandle().size(); ++k) {
    CFLog(INFO, "qy[" <<k << "] => (" << socket_qy.getDataHandle()[k] << "\n");
  }
  for (CFuint k = 0; k < socket_qz.getDataHandle().size(); ++k) {
    CFLog(INFO, "qz[" <<k << "] => (" << socket_qz.getDataHandle()[k] << "\n");
  }
  exit(1);*/
  CFLog(VERBOSE, "RadiativeTransferFVDOM::loopOverBins() => END\n");
}
      
//////////////////////////////////////////////////////////////////////////////

void RadiativeTransferFVDOM::loopOverDirs(const CFuint startBin, 
					  const CFuint endBin, 
					  const CFuint startDir,
					  const CFuint endDir)
{
  CFLog(VERBOSE, "RadiativeTransferFVDOM::loopOverDirs() => START\n");
  for (CFuint d = startDir; d < endDir; ++d) {
    CFLog(INFO, "( dir: " << d << " ), ( bin: ");
    const CFuint bStart = (d != startDir) ? 0 : startBin;
    const CFuint bEnd   = (d != m_startEndDir.second) ? m_nbBins : endBin;
    // precompute dot products for all faces and directions (a part from the sign)
    computeDotProdInFace(d, m_dotProdInFace);
    
    for(CFuint ib = startBin; ib < endBin; ++ib) {
      // old algorithm: opacities computed for all cells at once for a given bin
      CFLog(INFO, ib << " ");
      if (m_oldAlgo) {getFieldOpacities(ib);}
      
      (m_useExponentialMethod) ? 
	computeQExponential(ib,startDir,d) : computeQNoExponential(ib,startDir,d);
    }
    CFLog(INFO, ")\n");
  }
  CFLog(VERBOSE, "RadiativeTransferFVDOM::loopOverDirs() => END\n");
}
    
//////////////////////////////////////////////////////////////////////////////

void RadiativeTransferFVDOM::writeTGSData()
{
  // AL: THIS NEEDS TO BE PARALLELIZED. see exampe of reduceHeatFlux(), 
  // OTHERWISE IT CANNOT WORK IN PARALLEL LIKE THIS
  
  CFLog(VERBOSE, "RadiativeTransferFVDOM::writeTGSData() = > Writing stag line data for TGS testcase => start\n");
  
  boost::filesystem::path file = m_dirName / boost::filesystem::path("TGSData.plt");
  file = PathAppender::getInstance().appendParallel( file );
  
  SelfRegistPtr<Environment::FileHandlerOutput> fhandle = 
    Environment::SingleBehaviorFactory<Environment::FileHandlerOutput>::getInstance().create();
  ofstream& outputFile = fhandle->open(file);
  
  outputFile << "TITLE  = RadiativeTransferFVDOM data for TGS\n";
  outputFile << "VARIABLES = x  qx divq nbPoints\n";

  CellTrsGeoBuilder::GeoData& geoData = m_geoBuilder.getDataGE();
  Common::SafePtr<TopologicalRegionSet> cells = geoData.trs;
  
  const CFuint nbCells = cells->getLocalNbGeoEnts();
  CFuint nbPoints = 0;
  CFreal xCoord = 0.;
  const CFreal length = 0.06;
  DataHandle<CFreal> divQ = socket_divq.getDataHandle();
  DataHandle<CFreal> qx = socket_qx.getDataHandle();
  //DataHandle<CFreal> qy = socket_qy.getDataHandle();
  //DataHandle<CFreal> qz = socket_qz.getDataHandle();
  
  for(CFuint ir = 0; ir < m_Nr; ir++){
    nbPoints = 0;
    xCoord = -(ir + 0.5)*length/m_Nr; //middle point between ir and (ir + 1)
    
    for(CFuint iCell = 0; iCell < nbCells; iCell++){
      geoData.idx = iCell;
      GeometricEntity* currCell = m_geoBuilder.buildGE();
      
      const Node& coordinate = currCell->getState(0)->getCoordinates();
      const CFreal x = coordinate[XX];
      //const CFreal y = coordinate[YY];
      //const CFreal z = coordinate[ZZ];
      //const CFreal rCell = std::sqrt(x*x + y*y + z*z);
      
      if(abs(x) >= ir*length/m_Nr && abs(x) < (ir + 1)*length/m_Nr){
	nbPoints++;
	m_divqAv[ir] += divQ[iCell];
	m_qrAv[ir]   += qx[iCell]; //*rCell*rCell; Multiply by r**2 for area-weighted average
      }
      m_geoBuilder.releaseGE();
    }
    if(nbPoints > 0){
      m_divqAv[ir] /= nbPoints;
      m_qrAv[ir]   /= nbPoints; //m_qrAv[ir]   /= nbPoints*rCoord*rCoord; //area-weighted average radial flux
      outputFile << xCoord << " " << m_qrAv[ir] << " " << m_divqAv[ir] << " " <<  nbPoints << "\n";
    }
  }
  fhandle->close();
  
  CFLog(VERBOSE, "RadiativeTransferFVDOM::writeTGSData() => Writing stag line data for TGS testcase => end\n");
}

//////////////////////////////////////////////////////////////////////////////

void RadiativeTransferFVDOM::reduceHeatFlux()
{
  if (PE::GetPE().GetProcessorCount(m_radNamespace) > 1) {
    CFLog(VERBOSE, "RadiativeTransferFVDOM::reduceHeatFlux() => start\n");
    
    // AL: the following has to be moved to a postprocess() or something function
    //     to be called after execute() unless the ConcurrentCoupler can be modified to
    //     take care of it
    DataHandle< CFreal > qradFluxWall = socket_qradFluxWall.getDataHandle(); 
    const CFuint qsize = qradFluxWall.size();
    if (qsize > 0) { 
      PE::GetPE().setBarrier(m_radNamespace);
      
      vector<CFreal> localHeatFlux(qsize);
      for (CFuint i = 0; i < qsize; ++i) {
        localHeatFlux[i] = qradFluxWall[i];
      }
      
      MPIError::getInstance().check
	("MPI_Allreduce", "RadiativeTransferFVDOM::reduceHeatFlux()",
	 MPI_Allreduce(&localHeatFlux[0], &qradFluxWall[0], qsize, 
		       MPIStructDef::getMPIType(&localHeatFlux[0]), 
		       MPI_SUM, PE::GetPE().GetCommunicator(m_radNamespace)));
    }
    
    /*DataHandle<CFreal> divQ = socket_divq.getDataHandle();
    DataHandle<CFreal> qx   = socket_qx.getDataHandle();
    DataHandle<CFreal> qy   = socket_qy.getDataHandle();
    DataHandle<CFreal> qz   = socket_qz.getDataHandle();
    
    const CFuint qsize2 = divQ.size();
    vector<CFreal> localQ(qsize2);
    
    PE::GetPE().setBarrier(m_radNamespace);
    
    CFLog(VERBOSE, "RadiativeTransferFVDOM::reduceHeatFlux() => reducing divQ\n");
    for (CFuint i = 0; i < qsize2; ++i) {localQ[i] = divQ[i];}
    MPIError::getInstance().check
      ("MPI_Allreduce", "RadiativeTransferFVDOM::reduceHeatFlux()",
       MPI_Allreduce(&localQ[0], &divQ[0], qsize2, MPIStructDef::getMPIType(&localQ[0]), 
		     MPI_SUM, PE::GetPE().GetCommunicator(m_radNamespace)));
    
    PE::GetPE().setBarrier(m_radNamespace);
    
    CFLog(VERBOSE, "RadiativeTransferFVDOM::reduceHeatFlux() => reducing qx\n");
    for (CFuint i = 0; i < qsize2; ++i) {localQ[i] = qx[i];}
    MPIError::getInstance().check
      ("MPI_Allreduce", "RadiativeTransferFVDOM::reduceHeatFlux()",
       MPI_Allreduce(&localQ[0], &qx[0], qsize2, MPIStructDef::getMPIType(&localQ[0]), 
		     MPI_SUM, PE::GetPE().GetCommunicator(m_radNamespace)));
    
    PE::GetPE().setBarrier(m_radNamespace);
    
    CFLog(VERBOSE, "RadiativeTransferFVDOM::reduceHeatFlux() => reducing qy\n");
    for (CFuint i = 0; i < qsize2; ++i) {localQ[i] = qy[i];}
    MPIError::getInstance().check
      ("MPI_Allreduce", "RadiativeTransferFVDOM::reduceHeatFlux()",
       MPI_Allreduce(&localQ[0], &qy[0], qsize2, MPIStructDef::getMPIType(&localQ[0]), 
		     MPI_SUM, PE::GetPE().GetCommunicator(m_radNamespace)));
    
    PE::GetPE().setBarrier(m_radNamespace);
    
    CFLog(VERBOSE, "RadiativeTransferFVDOM::reduceHeatFlux() => reducing qz\n");
    for (CFuint i = 0; i < qsize2; ++i) {localQ[i] = qz[i];}
    MPIError::getInstance().check
      ("MPI_Allreduce", "RadiativeTransferFVDOM::reduceHeatFlux()",
       MPI_Allreduce(&localQ[0], &qz[0], qsize2, MPIStructDef::getMPIType(&localQ[0]), 
		     MPI_SUM, PE::GetPE().GetCommunicator(m_radNamespace)));
    */
    
    CFLog(VERBOSE, "RadiativeTransferFVDOM::reduceHeatFlux() => end\n");
  }
}

//////////////////////////////////////////////////////////////////////////////

/*void RadiativeTransferFVDOM::computeWallHeatFlux()
{
  boost::filesystem::path file = m_dirName / boost::filesystem::path("wallHeatFlux.out");
  file = PathAppender::getInstance().appendParallel( file );
  
  SelfRegistPtr<Environment::FileHandlerOutput> fhandle = 
    Environment::SingleBehaviorFactory<Environment::FileHandlerOutput>::getInstance().create();
  ofstream& outputFile = fhandle->open(file);

  outputFile << "TITLE  = RadiativeTransferFVDOM wallHeatFlux\n";
  outputFile << "VARIABLES = x y z q\n";

  CFLog(VERBOSE, "RadiativeTransferFVDOM::computeWallHeatFlux() => START\n");

  DataHandle< CFreal > qradFluxWall = socket_qradFluxWall.getDataHandle();
  DataHandle<CFreal> faceCenters = socket_faceCenters.getDataHandle();
  DataHandle<CFreal> faceAreas = socket_faceAreas.getDataHandle();
  DataHandle<CFreal> divQ = socket_divq.getDataHandle();
  DataHandle<CFreal> qx = socket_qx.getDataHandle();
  DataHandle<CFreal> qy = socket_qy.getDataHandle();
  DataHandle<CFreal> qz = socket_qz.getDataHandle();

  vector< string > wallTrsNames;  //to be moved to the .hh
  m_radiation->getWallTRSnames(wallTrsNames); //to be moved to the setup()
  const CFuint nbWallTrs = wallTrsNames.size();
  
  CFLog(INFO, "RadiativeTransferFVDOM::computeWallHeatFlux() => Number of walls for heat flux computation = [" << nbWallTrs << "]\n");

  RealVector distance(2);
  //RealVector divq_w(2); //wall values
  RealVector qx_w(2);
  RealVector qy_w(2);
  RealVector qz_w(2);
  CFreal tot_dist;
  RealVector q_tot(nbWallTrs); //storing the total heat flux on each wall
  q_tot = 0.;

  const CFuint m_dim = 3; //3D mesh
  RealVector faceCoord(m_dim);

  for (CFuint i=0; i<nbWallTrs; ++i){ //for each wall
    Framework::FaceTrsGeoBuilder::GeoData& facesData = m_wallFaceBuilder.getDataGE();
    SafePtr<TopologicalRegionSet> wallFaces = MeshDataStack::getActive()->getTrs(wallTrsNames[i]);
    const CFuint nbFacesWall = wallFaces->getLocalNbGeoEnts();
    facesData.trs = wallFaces;

    CFLog(INFO, "RadiativeTransferFVDOM::computeWallHeatFlux() => Wall TRS number [" << i << "] has [" << nbFacesWall << "] faces\n");

    for(CFuint f=0; f<nbFacesWall; ++f){
      facesData.idx = f;
      Framework::GeometricEntity *const face = m_wallFaceBuilder.buildGE();
      CFuint faceID = face->getID();
      const TopologicalRegionSet& faceTrs = *m_mapGeoToTrs->getTrs(faceID);
      const CFuint faceIdx = m_mapGeoToTrs->getIdxInTrs(faceID); //taking the local idx for the faceID within the TRS
      const CFuint faceGeoID = m_radiation->getCurrentWallGeoID();    

      m_wallFaceBuilder.releaseGE();

      for(CFuint ii=0; ii<m_dim; ++ii) {
	faceCoord[ii] = faceCenters[m_dim*faceGeoID + ii];
      }

      const CFreal faceArea = faceAreas[faceID];
      
      CellTrsGeoBuilder::GeoData& geoData = m_geoBuilder.getDataGE();
      Common::SafePtr<TopologicalRegionSet> cells = geoData.trs;
      
      CFreal tot_dist = 0.;
      //exploring the neighborhood of the current face
      for(CFuint ii=0; ii<2; ++ii) {
	CFuint cellID = faceTrs.getStateID(faceIdx, ii); //neighbor's ID
	//const CFuint cellID2 = faceTrs.getStateID(faceIdx, 1);
	
	//divq_w[ii] = divQ[cellID];
	qx_w[ii] = qx[cellID];
	qy_w[ii] = qy[cellID];
	qz_w[ii] = qz[cellID];
	
	geoData.idx = cellID;
	GeometricEntity* currCell = m_geoBuilder.buildGE();
	
	Node& coordinate = currCell->getState(0)->getCoordinates(); //coordinate[XX,YY,ZZ]
	distance[ii] = std::sqrt((faceCoord[0]-coordinate[XX])*(faceCoord[0]-coordinate[XX])+
				 (faceCoord[1]-coordinate[YY])*(faceCoord[1]-coordinate[YY])+
				 (faceCoord[2]-coordinate[ZZ])*(faceCoord[2]-coordinate[ZZ]));
	tot_dist += distance[ii];
	
	m_geoBuilder.releaseGE();
      }
      
      //CFreal divq_face = 0.;
      CFreal qx_face = 0.;
      CFreal qy_face = 0.;
      CFreal qz_face = 0.;

      for(CFuint ii=0; ii<2; ++ii) {
        cf_assert(tot_dist > 0.);
	//divq_face += divq_w[ii]*distance[ii]/tot_dist; //distance-weighted average
        qx_face += qx_w[ii]*distance[ii]/tot_dist;
        qy_face += qy_w[ii]*distance[ii]/tot_dist;
        qz_face += qz_w[ii]*distance[ii]/tot_dist;
      }
      
      const CFreal q_face = std::sqrt(qx_face*qx_face+qy_face*qy_face+qz_face*qz_face);
      q_tot[i] += q_face*faceArea;
      
      //CFuint idx = getWallFaceID(faceID); / == faceIdx
      //CFLog(INFO, " " << diff << "\n"); //to check the index
      qradFluxWall[faceIdx] = q_face;
      //qxFluxWall[idx] = qx_face;
      //qyFluxWall[idx] = qy_face;
      //qzFluxWall[idx] = qz_face;
      
      outputFile << faceCoord[0] << " " << faceCoord[1] << " " << faceCoord[2] << " " << q_face  <<"\n";
      //CFLog(INFO, divq_w[0] << " " << divq_w[1] << " " << q_face  <<"\n");
      
    }
    CFLog(INFO, "RadiativeTransferFVDOM::computeWallHeatFlux() => Total heat flux on the wall number [" << i  << "] = " << q_tot[i] << "\n");
  }
  
  fhandle->close();

  CFLog(VERBOSE, "RadiativeTransferFVDOM::computeWallHeatFlux() => END\n");
}
*/

//////////////////////////////////////////////////////////////////////////////

    } // namespace RadiativeTransfer

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

