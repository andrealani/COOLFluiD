#include <cmath>
#include <algorithm>
#include <boost/progress.hpp>
#include <boost/random.hpp>
#include "MathTools/MathChecks.hh"
#include "MathTools/MathConsts.hh"
#include "RadiativeTransferSanna.hh"
#include "RadiativeTransferSanna/RadiativeTransferMonteCarlo.hh"
#include "Framework/DataProcessingData.hh"
#include "Framework/ConvectiveVarSet.hh"
#include "Framework/DataSocketSource.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/GeometricEntity.hh"
#include "Framework/GeometricEntityPool.hh"
#include "Framework/FaceTrsGeoBuilder.hh"
#include "Framework/CellTrsGeoBuilder.hh"
#include "Framework/RadiationLibrary.hh"
#include "Framework/PhysicalConsts.hh"
#include "Common/CFPrintContainer.hh"
#include "Common/MPI/MPIError.hh"
#include "FiniteVolume/CellCenterFVM.hh"
#include "MathTools/MathFunctions.hh"
#define CELLK 1.0

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Numerics::FiniteVolume;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

namespace RadiativeTransfer {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<RadiativeTransferMonteCarlo, DataProcessingData, RadiativeTransfer> 
RadiativeTransferMonteCarloFVMCCProvider("RadiativeTransferMonteCarloFVMCC");

//////////////////////////////////////////////////////////////////////////////

void RadiativeTransferMonteCarlo::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< string >("RadiationLibrary","Name of the radiation library.");
  options.addConfigOption< vector<string> >("wallTrsNames","wall Trs Names.");
  options.addConfigOption< vector<string> >("boundaryTrsNames","boundary Trs Names.");
  options.addConfigOption< vector<string> >("symmetryTrsNames","symmetry Trs Names.");
  options.addConfigOption< CFuint >("numberOfRays","number of rays sent by each elemet.");
  //options.addConfigOption< CFuint >("NbPoints","Number of points.");
  options.addConfigOption< CFuint >("TestcaseID","Identifier for the testcase (1- 3D cylinder, 2- slab, 3 axisymmetric cylinder).");
  options.addConfigOption< CFuint >("MaxNbVisitedCells","Maximum number of visited cells.");
  options.addConfigOption< CFreal >("WallEmissivity","Emissivity at the wall.");
  options.addConfigOption< CFreal >("WallAbsorption","Absorption at the wall.");
  options.addConfigOption< CFuint, Config::DynamicOption<> >("FreezePattern","Freeze the statistical pattern.");
  options.addConfigOption< bool >("Axi","True if it is an axisymmetric simulation.");
  options.addConfigOption< bool >("PlanckFunction","Flag to tell iof to use Planck function at the wall");
  //options.addConfigOption< CFuint >("nCy","Set it for axisymmetric testcase, it is the number of cell in y direction");
  //options.addConfigOption< CFuint >("nCx","Set it for axisymmetric testcase, it is the number of cell in x direction");
  options.addConfigOption< CFuint >("ReducedSpectralSize","Number of points to use for reducing the spectra from the radiation library");
  options.addConfigOption< CFuint >("TemperatureID","Variable ID corresponding to the temperature");
  options.addConfigOption< CFreal >("FreeStreamTemperature","Temperature (typically free stream) for which qrad is assigned to 0");
  options.addConfigOption< CFuint >("NumberOfSteps","Number of post-prossing divisions");
}

//////////////////////////////////////////////////////////////////////////////

RadiativeTransferMonteCarlo::RadiativeTransferMonteCarlo(const std::string& name):
  DataProcessingCom(name),
  _particleTracking(),
  //m_MTRand(time(NULL)),
  socket_qrad("qrad"),
  socket_qradFluxWall("qradFluxWall"),
  socket_normals("normals"),
  socket_nodes("nodes"),
  socket_states("states"),
  socket_gstates("gstates"),
  socket_nstates("nstates"),
  socket_volumes("volumes"),
  socket_faceAreas("faceAreas"),
  socket_rankPartitionFaces("rankPartitionFaces"),
  socket_isOutward("isOutward"),
  m_nstatesProxy(CFNULL),
  m_nodeIdToStateId(),
  m_radCoeff(),
  m_radCoeffReduced(),
  absorptionCoeff(),
  absorptionCoeff2(),
  m_deltaWav(),
  m_deltaWavs(),
  _actualPoint(),
  _intersectionPointOld(),
  _internalPoint(),
  _intersectionPoint(),
  _reflectedDirection(),
  _intersection(),
  _nn(),
  m_midPointE(),
  m_midPointA(),
  m_diffAE(),
  m_projARandom(),
  _RHSGas(),
  _RHSWall(),
  m_mapWallFaceIDtoCellID(),
  m_axiVolumes(),
  m_wallSurfaces(),
  _dim(),
  _boundaryIDs(),
  _wallIDs(),
  _symmetryIDs(),
  _partitionIDs(),
  _myP(),
  _nP(),
  _comm(),
  _ProcessorOfPartitionGhostCell(),
  _Global2LocalInCellIdMap(),
  _CellIDmap(),
  _nWallFaces(),
  m_spectralIdx(0),
  m_iWavRange(0),
  m_cellCompletionIdx(),
  m_faceCompletionIdx(),
  _EmittingCellIDs(),
  _AbsorbingCellIDs(),
  _AbsorbingCellRanks(),
  _AbsorbingCellType(),
  _EmittingWallFaceIDs(),
  _AbsorbingWallFaceIDs(),
  _AbsorbingWallFaceRanks(),
  _AbsorbingWallFaceType(),
  m_sendCounts(),
  m_recvCounts(),
  m_sendDispls(),
  m_recvDispls(),
  m_isOverlapCell(),
  m_isWallFaceAbsorbing(),
  m_overlapCellRanks(),
  _CellG2LIDmap(),
  m_cellBuilder(),
  m_wallFaceBuilder(),
  m_radLibrary()
{
  addConfigOptionsTo(this);
  
  m_radLibraryName = "Null";
  setParameter("RadiationLibrary",&m_radLibraryName);
  
  setParameter("wallTrsNames",&_wallTrsNames);
  
  setParameter("boundaryTrsNames",&_boundaryTrsNames);
  
  setParameter("symmetryTrsNames",&_symmetryTrsNames);
  
  setParameter("numberOfRays",&_NRAYS);
  
  //_nbPoints = 16;
  //setParameter("NbPoints",&_nbPoints);
  
  _testcaseID = 0;
  setParameter("TestcaseID",&_testcaseID);

  _maxVisitedCells = 10000;
  setParameter("MaxNbVisitedCells",&_maxVisitedCells);

  m_freezePattern = 0;
  setParameter("FreezePattern", &m_freezePattern);

  m_nbSteps=30;
  setParameter("NumberOfSteps",&m_nbSteps);
  //_nCx = 64;
  //setParameter("nCx",&_nCx);

  //_nCy = 64;
  //setParameter("nCy",&_nCy);
  
  _Axi = false;
  setParameter("Axi",&_Axi);
  
  m_usePlanck = false;
  setParameter("PlanckFunction",&m_usePlanck);
  
  m_wallEm = 1.;
  setParameter("WallEmissivity",&m_wallEm);
  
  m_wallAbs = 1.;
  setParameter("WallAbsorption",&m_wallAbs);
  
  m_reducedSpectraSize = 0;
  setParameter("ReducedSpectralSize",&m_reducedSpectraSize);
  
  m_temperatureID = 0;
  setParameter("TemperatureID",&m_temperatureID);
  
  m_freeStreamT = 0.;
  setParameter("FreeStreamTemperature",&m_freeStreamT);
}

/////////////////////////////////////////////////////////////////////////////

RadiativeTransferMonteCarlo::~RadiativeTransferMonteCarlo()
{
}

/////////////////////////////////////////////////////////////////////////////

void RadiativeTransferMonteCarlo::configure(Config::ConfigArgs& args)
{
  CFAUTOTRACE;
  
  DataProcessingCom::configure(args);

  m_radLibrary = Environment::Factory<RadiationLibrary>::getInstance().
      getProvider(m_radLibraryName)->create(m_radLibraryName);
  cf_assert(m_radLibrary.isNotNull());
  configureNested ( m_radLibrary.getPtr(), args );
}

//////////////////////////////////////////////////////////////////////////////

vector<SafePtr<BaseDataSocketSink> > RadiativeTransferMonteCarlo::needsSockets(){

  vector<Common::SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_nodes);
  result.push_back(&socket_states);
  result.push_back(&socket_nstates);
  result.push_back(&socket_gstates);
  result.push_back(&socket_normals);
  result.push_back(&socket_volumes);
  result.push_back(&socket_faceAreas);
  result.push_back(&socket_rankPartitionFaces);
  result.push_back(&socket_isOutward);
  
  return result;
}

//////////////////////////////////////////////////////////////////////////////

vector<SafePtr<BaseDataSocketSource> > RadiativeTransferMonteCarlo::providesSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSource> > result;
  result.push_back(&socket_qrad);
  result.push_back(&socket_qradFluxWall);
  
  return result;
}

/////////////////////////////////////////////////////////////////////////////

void RadiativeTransferMonteCarlo::setup()
{ 
  // set up the radiation library
  m_radLibrary->setup();

  DataHandle<State*, GLOBAL> states = socket_states.getDataHandle();
  m_nodeIdToStateId.resize(states.size());
  for (CFuint i = 0; i < m_nodeIdToStateId.size(); ++i) {
    m_nodeIdToStateId[i] = i;
  }
  
  m_nstatesProxy.reset
      (new DofDataHandleIterator<CFreal, State, GLOBAL>(states,&m_nodeIdToStateId));
  
  // resizing of the array storing absorption and emission coefficients
  m_radCoeff.resize(states.size(), m_radLibrary->getWavelengthStride()*3);
  m_radCoeff = 0.;
  
  if (m_reducedSpectraSize > 0) {
    // number of intervals is user-defined m_reducedSpectraSize
    m_deltaWavs.resize(m_reducedSpectraSize, 1.);
    
    // AL: don't change the order of instructions here: it is meant to be like that!
    // count first and last point in the spectra as separate points
    //   m_reducedSpectraSize += 2;
    
 //  cf_always_assert(m_reducedSpectraSize < m_radLibrary->getWavelengthStride());
    // resizing of the array storing absorbtion and emission coefficients
    m_radCoeffReduced.resize(states.size(), m_reducedSpectraSize*3);
    m_radCoeffReduced = 0.;

    m_wavReduced.resize(states.size(),m_reducedSpectraSize);
    m_emReduced.resize(states.size(),m_reducedSpectraSize);
    m_amReduced.resize(states.size(),m_reducedSpectraSize);


    absorptionCoeff.resize(states.size(),m_reducedSpectraSize*3);
    absorptionCoeff=0.;

    absorptionCoeff2.resize(states.size());
  }
  else {
    m_deltaWavs.resize(m_radLibrary->getWavelengthStride(),1.);
  }
  
  if (_testcaseID > 0) {cf_always_assert(m_deltaWavs.size() == 1);}
  cf_assert(m_deltaWavs.size() > 0);

  // MPI parameters
  _myP = PE::GetPE().GetRank();
  _nP = PE::GetPE().GetProcessorCount();
  _comm = PE::GetPE().GetCommunicator();
  MPIError::getInstance().init(_comm, _myP);

  // set dimensions
  _dim = PhysicalModelStack::getActive()->getDim();
  if (m_temperatureID == 0) {
    CFLog(INFO, "RadiativeTransferMonteCarlo::setup() => TemperatureID set to default\n");
    m_temperatureID = _dim + 1;
  }

  _actualPoint.resize(_dim);
  _intersectionPointOld.resize(_dim);
  _internalPoint.resize(_dim);
  _intersectionPoint.resize(_dim);
  _reflectedDirection.resize(_dim);
  _intersection.resize(_dim);
  _nn.resize(_dim);
  
  m_midPointE.resize(_dim);
  m_midPointA.resize(_dim);
  m_diffAE.resize(_dim);
  m_projARandom.resize(_dim);
  
  // initialize ParticleTracking
  _particleTracking.setDataSockets(socket_states, socket_gstates, socket_nodes, socket_nstates, socket_normals,socket_isOutward);
  
  // setting up the cell builder
  m_cellBuilder.setup();
  m_cellBuilder.getGeoBuilder()->setDataSockets(socket_states, socket_gstates, socket_nodes);
  CellTrsGeoBuilder::GeoData& cellData = m_cellBuilder.getDataGE();
  SafePtr<TopologicalRegionSet> cells = MeshDataStack::getActive()->getTrs("InnerCells");
  cellData.trs = cells;
  
  // setting up the wall face builder
  m_wallFaceBuilder.setup();
  m_wallFaceBuilder.getGeoBuilder()->setDataSockets(socket_states, socket_gstates, socket_nodes);
  m_wallFaceBuilder.getDataGE().isBFace = true;
  
  const CFuint nCells = cellData.trs->getLocalNbGeoEnts();
  m_cellCompletionIdx.resize(nCells);
  
  _RHSGas.resize(nCells, 0.);
  
  // preallocation of memory for _RHSWall
  FaceTrsGeoBuilder::GeoData& WallFacesData = m_wallFaceBuilder.getDataGE();
  for(CFuint j=0; j<_wallTrsNames.size(); ++j){
    SafePtr<TopologicalRegionSet> WallFaces = MeshDataStack::getActive()->getTrs(_wallTrsNames[j]);
    WallFacesData.trs = WallFaces;
    if (_RHSWall.size() + WallFaces->getLocalNbGeoEnts() > 0) {
      _RHSWall.resize(_RHSWall.size() + WallFaces->getLocalNbGeoEnts(), 0.);
    }
  }
  
  m_mapWallFaceIDtoCellID.reserve(_RHSWall.size());
  
  if (_Axi) {
    m_axiVolumes.resize(states.size(), 0.);
  }
  m_wallSurfaces.resize(_RHSWall.size(), 0.);
  
  // resize the storage of qrad in each cell
  socket_qrad.getDataHandle().resize(nCells);
  if (_RHSWall.size() > 0) {
    socket_qradFluxWall.getDataHandle().resize(_RHSWall.size());
  }
  
  _EmittingCellIDs.reserve(nCells);

  const CFuint nray1cells = (_NRAYS + 1)*nCells;
  _AbsorbingCellIDs.resize(nray1cells,-1);
  _AbsorbingCellRanks.resize(nray1cells,-1);
  _AbsorbingCellType.resize(nray1cells,-1);
  
  // other tasks
  buildCellIDradiusmap();
  BoundaryKnowledge();
  
  // boundary faces are counted before this point inside BoundaryKnowledge()
  // it is possible that some partition domains do not have wall faces
  if (_nWallFaces > 0) {
    m_faceCompletionIdx.resize(_nWallFaces,0);
    _EmittingWallFaceIDs.reserve(_nWallFaces);
    const CFuint nray1wall = (_NRAYS + 1)*_nWallFaces;
    _AbsorbingWallFaceIDs.resize(nray1wall,-1);
    _AbsorbingWallFaceRanks.resize(nray1wall,-1);
    _AbsorbingWallFaceType.resize(nray1wall,-1);
  }
  
  cf_assert(_RHSWall.size() == _nWallFaces);
  
  m_sendCounts.resize(_nP);
  m_recvCounts.resize(_nP);
  m_sendDispls.resize(_nP);
  m_recvDispls.resize(_nP);
  cf_assert(m_sendCounts.size() == _nP);
  cf_assert(m_recvCounts.size() == _nP);
  cf_assert(m_sendDispls.size() == _nP);
  cf_assert(m_recvDispls.size() == _nP);
  
  buildCellIDmap();

  flagOverlapCells();
  
  const CFuint nbFaces = MeshDataStack::getActive()->Statistics().getNbFaces();
  m_isWallFaceAbsorbing.resize(nbFaces, false);
  
  // for testcases 1, 2 or 3, the following is precomputed
  if (_testcaseID == 1 || _testcaseID == 2 || _testcaseID == 3) {
    MonteCarlo();
  }
}

/////////////////////////////////////////////////////////////////////////////

void RadiativeTransferMonteCarlo::executeOnTrs()
{ 
  CFLog(VERBOSE, "RadiativeTransferMonteCarlo::executeOnTrs() START\n");
  
  // for testcases == 0, the following must be recomputed everytime
  if (_testcaseID == 0) {
    Stopwatch<WallTime> s;
    s.restart();

    // re-initialize the heat fluxes to 0
    _RHSGas.assign(_RHSGas.size(), 0.);
    if (_RHSWall.size() > 0) {
      _RHSWall.assign(_RHSWall.size(), 0.);
    }
    for(CFuint i=0;i<absorptionCoeff2.size();++i){
      absorptionCoeff2[i] = new typeAbsMap;
    }

    const CFuint nbWavLoops = m_radLibrary->getWavLoopSize();
    CFLog(VERBOSE, " RadiativeTransferMonteCarlo::executeOnTrs() => nbWavLoops = " << nbWavLoops << "\n");
    for (m_iWavRange = 0; m_iWavRange < nbWavLoops; ++m_iWavRange) {
      m_radLibrary->computeProperties(m_nstatesProxy.get(), m_radCoeff, m_iWavRange);

      // if the user defines a value of ReducedSpectralSize>0, then the spectra are reduced before running MC
      if (m_radCoeffReduced.size() > 0) myReduceSpectra(m_radCoeff, m_radCoeffReduced);
      
      CFuint wavStride = m_radLibrary->getWavelengthStride();
      if (m_radCoeffReduced.size() > 0) {
        // m_reducedSpectraSize-1 is considered because what counts is the number of wavelength intervals,
        // not the actual number of wavelengths
        wavStride = m_reducedSpectraSize;
      }
      
      cf_assert(wavStride <= m_radLibrary->getWavelengthStride());
      CFLog(VERBOSE, "RadiativeTransferMonteCarlo::executeOnTrs() => wavStride = " << wavStride << "\n");
      for (m_spectralIdx = 0; m_spectralIdx  < wavStride; ++m_spectralIdx) {
        CFLog(DEBUG_MAX, "RadiativeTransferMonteCarlo::executeOnTrs() => SpectralID = " << m_spectralIdx + wavStride*m_iWavRange << "\n");
        if (m_freezePattern == 0) {
          MonteCarlo();
        }

        CFLog(INFO, "spectralIdx = " << m_spectralIdx << " < " << m_deltaWavs.size() << "\n");
        computeHeatFlux();
      }
    }
    for(CFuint i=0;i<absorptionCoeff2.size();++i){
     delete absorptionCoeff2[i];
    }
    
    CFLog(DEBUG_MAX, "MonteCarlo() took " << s << "s\n");
  }
  else {
    // uncoupled cases
    computeHeatFlux();
    switch(_testcaseID) {
      case(1):
        getRHFcircularSections();
        break;
      case(2):
        getRHFslab();
        break;
      case(3):
        getRHFaxyCylinder();
        break;
      default:
        cout << "ERROR: RadiativeTransferMonteCarlo::executeOnTrs() => TestcaseID undefined!" << endl; abort();
        break;
    }
  }
  DataHandle<State*, GLOBAL> states = socket_states.getDataHandle();
  DataHandle<CFreal> gasRadiativeHeatSource = socket_qrad.getDataHandle();
  DataHandle<CFreal> volumes = socket_volumes.getDataHandle();
  cf_assert(gasRadiativeHeatSource.size() == _RHSGas.size());
  for(CFuint i=0; i<_RHSGas.size(); ++i){
    const CFreal volume = (_Axi) ? m_axiVolumes[i] : volumes[i];
    //CFreal volume=volumes[i];
    cf_assert(volume > 0.);
    gasRadiativeHeatSource[i] = _RHSGas[i]/volume; // this is divergence of qrad
    
    // this can be used to switch off radiation in the free stream
    if (m_freeStreamT > 0.) {
      const CFreal temp = (*states[i])[m_temperatureID];
      if ((temp > m_freeStreamT - 0.01) && (temp < m_freeStreamT + 0.01)) {
        gasRadiativeHeatSource[i] = 0.0;
      }
    }
  }
  
  
  DataHandle<CFreal> surfaces = socket_faceAreas.getDataHandle();
  DataHandle<CFreal> wallRadiativeHeatSource = socket_qradFluxWall.getDataHandle();
  cf_assert(wallRadiativeHeatSource.size() == _RHSWall.size());
  cf_assert(wallRadiativeHeatSource.size() == _nWallFaces);
  const CFuint nbWallTrs = _wallTrsNames.size();
  CFuint ff = 0;
  for(CFuint i=0; i<nbWallTrs; ++i){
    SafePtr<TopologicalRegionSet> WallFaces = MeshDataStack::getActive()->getTrs(_wallTrsNames[i]);
    const CFuint nbFacesWall = WallFaces->getLocalNbGeoEnts();
    for(CFuint f=0; f<nbFacesWall; ++f, ++ff){
      cf_assert(m_wallSurfaces[ff] > 0.);
      wallRadiativeHeatSource[ff] = _RHSWall[ff]/m_wallSurfaces[ff];
      // CFLog(INFO, " AFTER faceID = " << WallFaces->getLocalGeoID(f) << ", index = " << ff << "\n");
    }
  }
  
  // sanity check to see if all faces are absorbing some portion of radiative flux
  CFuint counter = 0;
  for (CFuint i = 0; i < _wallIDs.size(); ++i) {
    if (!m_isWallFaceAbsorbing[_wallIDs[i]]) {
      CFLog(INFO, "ERROR: Wall face with ID " <<_wallIDs[i] << " is NOT ABSORBING ANYTHING.\n");
      counter++;
    }
  }
  
  if (counter > 0) {
    CFLog(INFO, "ERROR: Number of wall faces not absorbing is " << counter << "/" << _nWallFaces << "\n");
  }
  
  CFLog(VERBOSE, "RadiativeTransferMonteCarlo::executeOnTrs() END\n");
}

/////////////////////////////////////////////////////////////////////////////

void RadiativeTransferMonteCarlo::BoundaryKnowledge()
{
  setFaceIDs(_wallTrsNames, _wallIDs, true);
  _nWallFaces = _wallIDs.size();
  
  setFaceIDs(_symmetryTrsNames, _symmetryIDs, false);
  
  // boundary faces but wall
  setFaceIDs(_boundaryTrsNames, _boundaryIDs, false);
  
  // partition faces
  SafePtr<TopologicalRegionSet> PartitionFaces = MeshDataStack::getActive()->getTrs("PartitionFaces");
  FaceTrsGeoBuilder::GeoData& facesData = m_wallFaceBuilder.getDataGE();
  facesData.trs = PartitionFaces;
  const CFuint nPF = PartitionFaces->getLocalNbGeoEnts();
  vector<pair<CFuint,CFuint> > PartitionIDs;
  PartitionIDs.reserve(nPF);
  for(CFuint i=0; i<nPF; ++i){
    facesData.idx = i;
    GeometricEntity *const face = m_wallFaceBuilder.buildGE();
    PartitionIDs.push_back(make_pair(face->getID(),i));
    m_wallFaceBuilder.releaseGE();
  }
  sort(PartitionIDs.begin(),PartitionIDs.end(), PairSort<CFuint,CFuint>());
  
  DataHandle<CFuint> rankPartitionFaces = socket_rankPartitionFaces.getDataHandle();
  assert(rankPartitionFaces.size() == nPF);
  
  _partitionIDs.reserve(nPF);
  for(CFuint i=0; i<nPF; ++i){
    _partitionIDs.push_back(PartitionIDs[i].first);
    const CFuint index = PartitionIDs[i].second;
    facesData.idx = index;
    GeometricEntity *const face = m_wallFaceBuilder.buildGE();
    // [i] or [index]
    // store the global ID of the right state of the current partition face (key) and the rank of the processor containing that state
    // AL: we are assuming that there is only one processor containing that state and this should be the case since that
    //     state should not belong to the overlap
    _ProcessorOfPartitionGhostCell.insert(face->getState(1)->getGlobalID(),rankPartitionFaces[index]);
    // cout<<" PartitionIDs[i].first = "<<PartitionIDs[i].first<<" face->getState(1)->getGlobalID() = "<<face->getState(1)->getGlobalID()<<" face->getState(0)->getGlobalID() = "<<face->getState(0)->getGlobalID()<<" rankPartitionFaces[index] = "<<rankPartitionFaces[index]<<" _myP = "<< _myP<<endl;
    assert(rankPartitionFaces[index] != _myP);
    // cout << face->getState(1)->getGlobalID() << " ----- " << rankPartitionFaces[index] << endl;
    m_wallFaceBuilder.releaseGE();
  }
  
  // store the mapping bewteen the global state IDs and the local state IDs
  DataHandle<State*, GLOBAL> states = socket_states.getDataHandle();
  for (CFuint i = 0; i < states.size(); ++i) {
    _Global2LocalInCellIdMap.insert(states[i]->getGlobalID(), states[i]->getLocalID());
  }
  _Global2LocalInCellIdMap.sortKeys();
}

/////////////////////////////////////////////////////////////////////////////

void RadiativeTransferMonteCarlo::computeCellRays(vector<vector<Ray> >& partitionBeams)
{ 
  CFLog(DEBUG_MAX, "RadiativeTransferMonteCarlo::computeCellRays()\n");
  _EmittingCellIDs.clearContent();
  CellTrsGeoBuilder::GeoData& cellData = m_cellBuilder.getDataGE();
  const CFuint nCells = cellData.trs->getLocalNbGeoEnts();
  const CFuint nray1cells = (_NRAYS + 1)*nCells;
  _AbsorbingCellIDs.assign(nray1cells,-1);
  _AbsorbingCellRanks.assign(nray1cells,-1);
  _AbsorbingCellType.assign(nray1cells,-1);

  m_rand.seed(time(NULL)*(_myP+1));
  
  CFuint o = 0; //number of photon absorbed by the emitting cell itself
  const string Trs = "InnerCells";
  EntityType entity;
  boost::progress_display* progressBar = NULL;
  if (_myP == 0) progressBar = new boost::progress_display(nCells);
  //stringstream ss;
  //ss<<"directions_"<<_myP<<".txt";
  //string ss2=ss.str();
  //ofstream debug(ss2.c_str());

  //  ofstream raytracingResults("raytracing.txt") ;
  for(CFuint c=0; c<nCells; ++c) {
    if (_myP == 0) {
      ++(*progressBar);
    }
    m_cellCompletionIdx[c] = 0;
    cellData.idx = c;
    GeometricEntity *const cell = m_cellBuilder.buildGE();
    //if (cell->getState(0)->isParUpdatable()){
      const CFuint cellID = cell->getState(0)->getGlobalID();

      _EmittingCellIDs.insert(cellID,c);
      CFreal minWav=0.,maxWav=0;
      if(m_usePlanck){
        minWav=absorptionCoeff(c,m_spectralIdx*3);
        maxWav=absorptionCoeff(c,m_spectralIdx*3+1);
      }
      for(CFuint r=0; r<_NRAYS; ++r){
        // build beam
        Ray beam;
        beam.tt=0.;
        CFuint ndirections = (_Axi) ? 3 : _dim;

        RealVector directions(ndirections);
        m_rand.sphereDirections(ndirections,directions);
        for(CFuint ii=0; ii<ndirections; ++ii){
          beam.direction[ii]=directions[ii];
          //cout<<directions[ii];
        }
        //cout<<endl;

        //     if (_Axi) {
        //     beam.direction[0]=-0.7747;
        //     beam.direction[1]=0.4293;
        //     beam.direction[2]=0.4643;
        //   }

        //calculate directions
        /*
      CFreal tetha, phi;
      phi=m_rand.uniformRand(0.,2.*MathConsts::CFrealPi());

      if(_dim==2){
        beam.direction[0]=cos(phi);
        beam.direction[1]=sin(phi);
      }
      else{
        tetha=m_rand.uniformRand(0.,MathConsts::CFrealPi());
        beam.direction[0]=sin(tetha)*cos(phi);
        beam.direction[1]=sin(tetha)*sin(phi);
         beam.direction[2]=cos(tetha);
      }

      cout<<endl<<"NEW BEAM"<<endl<<endl;
      cout<<"Beam emited on cell "<<c<<endl;
      cout<<"Direction: "<<beam.direction[0]<<' '<< beam.direction[1]<<' '<<beam.direction[2]<<endl;
      cout<<"Position: ";
      debug<<c<<' '<<beam.direction[0]<<' '<< beam.direction[1]<<' '<<beam.direction[2]<<endl;
      */
        //debug<<c<<' '<<beam.direction[0]<<' '<< beam.direction[1]<<' '<<beam.direction[2]<<endl;
        beam.startEntityID = cellID;
        beam.emittingEntityID = cellID;
        beam.emittingEntityType = -1;
        beam.emittingProcessorRank = _myP;
        beam.KS = - log(m_rand.uniformRand());
        beam.actualKS = 0;
        beam.wavelength=0.;
        if (m_usePlanck){
            beam.wavelength=m_rand.uniformRand(minWav,maxWav);
        }
        // rayTrace beam
        const CFuint AbsorbingEntityID = rayTracing(beam, 0, entity, Trs);
        //raytracingResults<<c<<' '<<AbsorbingEntityID<<' '<<entity<<endl;
        if(AbsorbingEntityID == cellID) o++;

        if(entity == WALL_FACE) {
          const CFuint idx = c*(_NRAYS +1) + m_cellCompletionIdx[c];
          assert(idx < _AbsorbingCellIDs.size());
          _AbsorbingCellIDs[idx]   = AbsorbingEntityID;
          _AbsorbingCellRanks[idx] = _myP;
          // here globalID of the internal cell adiacent to the face
          cf_assert(AbsorbingEntityID >= 0);
          _AbsorbingCellType[idx]  = m_mapWallFaceIDtoCellID.find(AbsorbingEntityID);
          m_cellCompletionIdx[c]++;
        }
        if (entity == NEGLIGIBLE){
          cout<<"particle is NEGLIGIBLE!!"<<endl;
        }
        if (entity == INTERNAL_CELL ){
          const CFuint idx = c*(_NRAYS + 1) + m_cellCompletionIdx[c];
          assert(idx < _AbsorbingCellIDs.size());
          _AbsorbingCellIDs[idx]   = AbsorbingEntityID;
          _AbsorbingCellRanks[idx] = _myP;
          _AbsorbingCellType[idx]  = -1;
          m_cellCompletionIdx[c]++;
        }

        if(entity == PARTITION_GHOSTCELL ){
          if (!_ProcessorOfPartitionGhostCell.exists(beam.startEntityID)) {
            CFLog(DEBUG_MAX, _myP << "RadiativeTransferMonteCarlo::computeCellRays() => beam.startEntityID " << beam.startEntityID << " doesn't exist\n");
            cout<<"error!"<<endl;
          }
          cf_assert(_ProcessorOfPartitionGhostCell.exists(beam.startEntityID));
          const CFuint rankP = _ProcessorOfPartitionGhostCell.find(beam.startEntityID);
          partitionBeams[rankP].push_back(beam);
        }
      }
    //}
    m_cellBuilder.releaseGE();
  }
  //debug.close();
  if (_myP == 0) delete progressBar;
  
  _EmittingCellIDs.sortKeys();

  // force synchronization at this point
  MPI_Barrier(_comm);


}

/////////////////////////////////////////////////////////////////////////////

void RadiativeTransferMonteCarlo::computeWallFaceRays(vector<vector<Ray> >& partitionBeams)
{     
  CFLog(DEBUG_MAX, "RadiativeTransferMonteCarlo::computeWallFaceRays()\n");
  
  _EmittingWallFaceIDs.clearContent();
  
  if (_nWallFaces > 0) {
    const CFuint nray1wall = (_NRAYS + 1)*_nWallFaces;
    _AbsorbingWallFaceIDs.assign(nray1wall,-1);
    _AbsorbingWallFaceRanks.assign(nray1wall,-1);
    _AbsorbingWallFaceType.assign(nray1wall,-1);
  }
  
  FaceTrsGeoBuilder::GeoData& WallFacesData = m_wallFaceBuilder.getDataGE();
  const CFuint nbWallTrs = _wallTrsNames.size();
  DataHandle<CFreal> normals = socket_normals.getDataHandle();
  
  // counter for the total number of wall faces, including all wall TRS's
  CFuint ff = 0;

  CFuint nbCell  = 0;
  CFuint nbFace  = 0;
  CFuint nbGhost = 0;
  CFuint nbNeg   = 0;

  for(CFuint j=0; j<nbWallTrs; ++j){
    SafePtr<TopologicalRegionSet> WallFaces = MeshDataStack::getActive()->getTrs(_wallTrsNames[j]);
    WallFacesData.trs = WallFaces;
    WallFacesData.isBFace = true;
    const string Trs = _wallTrsNames[j];
    
    const CFuint nbFacesWall = WallFaces->getLocalNbGeoEnts();
    for(CFuint f=0; f<nbFacesWall; ++f){
      cf_assert(ff < m_faceCompletionIdx.size());
      m_faceCompletionIdx[ff] = 0;
      WallFacesData.idx = f;
      GeometricEntity *const face = m_wallFaceBuilder.buildGE();
      const CFuint faceID = face->getID();
      const CFuint globalID = face->getState(0)->getGlobalID();
      const CFuint startID = faceID*_dim;
      _nn[0] = normals[startID];
      _nn[1] = normals[startID + 1];
      if(_dim == 3){
        _nn[2] = normals[startID + 2];
      }

      const CFreal invNLength = -1./_nn.norm2();     // normals always points outwards,
      _nn *= invNLength;                             // so have to reverse them
      
      vector<CFreal> nn2(_dim);
      for(CFuint i=0;i<_dim; ++i){
        nn2[i]=_nn[i];
      }
      //cout<<"FACE NORMAL: "<<nn2[0]<<' '<<nn2[1]<<' '<<nn2[2]<<endl;
      // CHANGE
      _EmittingWallFaceIDs.insert(globalID,ff);
      for(CFuint r=0; r<_NRAYS; ++r){

        // build beam
        Ray beam;
        if(!_Axi){
          vector<CFreal> directions(_dim);
          m_rand.hemiDirections(_dim,nn2,directions);
          for(CFuint d=0; d<_dim; ++d){
            beam.direction[d] = directions[d];
          }
         // cout<<"BEAM DIRECTION: "<<directions[0]<<' '<<directions[1]<<' '<<directions[2]<<endl;
        }

        beam.startEntityID = f; // use f instead of faceID: see the setupStartPoint function of particle tracking;

        beam.emittingEntityID = faceID; // AL: emitting entityID is never actually used !!
        beam.emittingEntityType = globalID;
        beam.emittingProcessorRank = _myP;
        beam.KS = - log(m_rand.uniformRand());
        beam.actualKS = 0.;
        beam.tt=0.;

        // WALL CELL RAYS ARE TURNED OFF!!
        EntityType entity = NEGLIGIBLE;
        const CFuint AbsorbingEntityID = 0;//rayTracing(beam, 2, entity, Trs);
        if (entity == WALL_FACE ){
          cf_assert(ff < m_faceCompletionIdx.size());
          const CFuint idx = ff*(_NRAYS + 1) + m_faceCompletionIdx[ff];
          _AbsorbingWallFaceIDs[idx]   = AbsorbingEntityID;
          _AbsorbingWallFaceRanks[idx] = _myP;
          // here globalID of the internal cell adiacent to the face
          cf_assert(AbsorbingEntityID >= 0);
          _AbsorbingWallFaceType[idx]  = m_mapWallFaceIDtoCellID.find(AbsorbingEntityID);
          m_faceCompletionIdx[ff]++;
          //cout<<"WALL FACE"<<endl;
          nbFace++;
        }

        if(entity == INTERNAL_CELL ){
          cf_assert(ff < m_faceCompletionIdx.size());
          const CFuint idx = ff*(_NRAYS + 1) + m_faceCompletionIdx[ff];
          _AbsorbingWallFaceIDs[idx]   = AbsorbingEntityID;
          _AbsorbingWallFaceRanks[idx] = _myP;
          _AbsorbingWallFaceType[idx]  = -1;
          m_faceCompletionIdx[ff]++;
          //cout<<"INTERNAL CELL"<<endl;
          nbCell++;
        }
        
        if(entity== PARTITION_GHOSTCELL ){
          cf_assert(_ProcessorOfPartitionGhostCell.exists(beam.startEntityID));
          if (_ProcessorOfPartitionGhostCell.exists(beam.startEntityID)) {
            const CFuint rankP = _ProcessorOfPartitionGhostCell.find(beam.startEntityID);
            partitionBeams[rankP].push_back(beam);
          }
          nbGhost++;
        }
        if(entity== NEGLIGIBLE ){
          nbNeg++;
          //cout<<"PARTICLE IS NEGLIGIBLE!!!"<<endl;
        }

      }
      m_wallFaceBuilder.releaseGE();
      ff++;
    }
  }
  //cout<<"Number of Neg: "<<nbNeg<<endl<<"Number of photons: "<<ff*_NRAYS<<endl;
  _EmittingWallFaceIDs.sortKeys();
  //cout<<"number of rays arriving to:"<<endl<<"Face: "<<nbFace<<endl<<"Cell: "<<nbCell<<endl<<"GhostFace: "<<nbGhost<<endl;
  
  // force synchronization at this point
  MPI_Barrier(_comm);
}

/////////////////////////////////////////////////////////////////////////////

void RadiativeTransferMonteCarlo::computePartitionFaceRays
(vector<vector<Ray> >& partitionBeams, 
 vector<vector<CFuint> >& returningPartitionBeams)
{   
  CFLog(DEBUG_MAX, "RadiativeTransferMonteCarlo::computePartitionFaceRays()\n");
  
  // create MPI ray structure
  MPI_Datatype MPI_Ray_ptr;
  build_MPI_strut_Ray(&MPI_Ray_ptr);
  
  // stop condition
  bool stopCondition = false;
  
  CFuint sendSize = 0;
  for(CFuint i=0; i<_nP; ++i){
    sendSize += partitionBeams[i].size();
  }
  
  // send and recv buffers are preallocated with size > 0
  vector<Ray> sendBuffer(max<CFuint>(sendSize,1));
  vector<Ray> recvBuffer(max<CFuint>(sendSize,1));
  
  // loop till total absorbtion or dispersion of photons achieved
  do {
    // communication operation
    
    sendBuffer.clear();
    for(CFuint i=0; i<_nP; ++i){
      const CFuint addedSize = partitionBeams[i].size();
      for(CFuint j=0; j<addedSize; ++j){
        sendBuffer.push_back(partitionBeams[i][j]);
      }
    }
    
    // set and verify stop condition
    CFuint localStop = 1;
    CFuint totalStop = 0;
    if(sendBuffer.size() == 0){ // if no one partition has to sent beam to any other (it means also that all photons have been absorbed or dispersed)
      localStop = 0;
    }
    MPI_Allreduce(&localStop, &totalStop, 1, MPI_UNSIGNED, MPI_MAX, _comm);
    if(totalStop == 0){
      stopCondition = true;
      break;
    }
    
    // build counts vectors
    for(CFuint p=0; p<_nP; p++){
      m_sendCounts[p] = partitionBeams[p].size();
    }
    MPI_Alltoall(&m_sendCounts[0], 1, MPI_UNSIGNED, &m_recvCounts[0], 1, MPI_UNSIGNED, _comm);
    
    // build displacements vectors
    m_recvDispls[0] = 0;
    m_sendDispls[0] = 0;
    
    CFuint totRecvCount = m_recvCounts[0];
    if(_nP>1){
      for(CFuint p=1; p<_nP; p++){
        m_recvDispls[p] = m_recvDispls[p-1] + m_recvCounts[p-1];
        m_sendDispls[p] = m_sendDispls[p-1] + m_sendCounts[p-1];
        totRecvCount += m_recvCounts[p];
      }
    }
    
    // cout << CFPrintContainer<vector<int> >("sendCounts  = ", &m_sendCounts) << endl;
    // cout << CFPrintContainer<vector<int> >("sendDispls  = ", &m_sendDispls) << endl;
    // cout << CFPrintContainer<vector<int> >("recvCounts  = ", &m_recvCounts) << endl;
    // cout << CFPrintContainer<vector<int> >("recvDispls  = ", &m_recvDispls) << endl;
    
    // build receive buffers
    const CFuint LastDisplacement = m_recvDispls[_nP-1] + m_recvCounts[_nP-1];
    cf_assert(totRecvCount == LastDisplacement);
    
    // the receiving buffers must be properly resized
    // vector<Ray> recvBuffer(LastDisplacement);
    // resize the receive buffer only if the new size is bigger
    if (LastDisplacement > recvBuffer.size()) {
      recvBuffer.resize(totRecvCount);
      // this leads to an inconsistent recvBuffer size which is >= the needed one
      // but avoids to have a 0 size array  which would cause out-of-bounds access
    }
    
    cf_assert(sendBuffer.capacity() > 0);
    cf_assert(recvBuffer.capacity() > 0);
    
    MPIError::getInstance().check("MPI_Alltoallv", "RadiativeTransferMonteCarlo::MonteCarlo()",
                                  MPI_Alltoallv(&sendBuffer[0], &m_sendCounts[0], &m_sendDispls[0], MPI_Ray_ptr,
                                                &recvBuffer[0], &m_recvCounts[0], &m_recvDispls[0], MPI_Ray_ptr, _comm));
    
    // clear content of partitionBeams in order to use it again
    for (CFuint p = 0; p < partitionBeams.size(); ++p) {
      partitionBeams[p].clear();
      cf_assert(partitionBeams[p].size() == 0);
    }
    // compute partitionFaces rays
    const string Trs = "InnerCells";
    for(CFuint i=0; i < LastDisplacement; ++i) {
      Ray& beam = recvBuffer[i];
      
      // rayTrace beam
      EntityType entity;
      const CFuint emitting_rankP = beam.emittingProcessorRank;
      //cout<<_myP<<" START PART FACE RAYS "<<endl;
      const CFuint AbsorbingEntityID = rayTracing(beam, 1, entity, Trs);
      //cout<<_myP<<" END PART FACE RAYS  RAYTRACING "<<endl;

      
      if (entity == WALL_FACE){
        returningPartitionBeams[emitting_rankP].push_back(beam.emittingEntityID);
        returningPartitionBeams[emitting_rankP].push_back(beam.emittingEntityType);
        returningPartitionBeams[emitting_rankP].push_back(AbsorbingEntityID);
        returningPartitionBeams[emitting_rankP].push_back(m_mapWallFaceIDtoCellID.find(AbsorbingEntityID));
        returningPartitionBeams[emitting_rankP].push_back(_myP);
        CFLog(DEBUG_MAX, "WALL_FACE\n");
      }
      
      if (entity == INTERNAL_CELL){
        returningPartitionBeams[emitting_rankP].push_back(beam.emittingEntityID);
        returningPartitionBeams[emitting_rankP].push_back(beam.emittingEntityType);
        returningPartitionBeams[emitting_rankP].push_back(AbsorbingEntityID);
        returningPartitionBeams[emitting_rankP].push_back(-1);
        returningPartitionBeams[emitting_rankP].push_back(_myP);
        CFLog(DEBUG_MAX, "INTERNAL_CELL\n");
      }

      if(entity== NEGLIGIBLE ){
        cout<<"PARTICLE IS NEGLIGIBLE!!!"<<endl;
      }
      
      /// to be fixed
      if (entity == PARTITION_GHOSTCELL) {
        if (!_ProcessorOfPartitionGhostCell.exists(beam.startEntityID)) {
          CFLog(INFO, "beam.startEntityID " << beam.startEntityID << " not found in _ProcessorOfPartitionGhostCell\n");
          cout<<"NOT FOUND!!"<<endl;
        }
        const CFuint rankP = _ProcessorOfPartitionGhostCell.find(beam.startEntityID);
        partitionBeams[rankP].push_back(beam);
        CFLog(DEBUG_MAX, "PARTITION_GHOSTCELL\n");
      }

      //cout<<_myP<<" END PARTITION FACE RAYS "<<endl;

    }

  } while(!stopCondition);

  // force synchronization at this point
  MPI_Barrier(_comm);
}

/////////////////////////////////////////////////////////////////////////////

void RadiativeTransferMonteCarlo::computeFinalCommunication
(std::vector<std::vector<CFuint> >& returningPartitionBeams)
{
  CFLog(DEBUG_MAX, "RadiativeTransferMonteCarlo::computeFinalCommunication()\n");
  
  // build send buffers
  CFuint sizeSendBuffer2 = 0;
  for(CFuint i=0; i<_nP; ++i){
    sizeSendBuffer2 += returningPartitionBeams[i].size();
  }
  
  // set and verify stop condition
  // if no one partition has to sent beam to any other (it means also that all photons has been absorbed or dispersed)
  CFuint totalStop = 0;
  CFuint localStop = (sizeSendBuffer2 == 0) ? 0 : 1;
  MPI_Allreduce(&localStop, &totalStop, 1, MPI_UNSIGNED, MPI_MAX, _comm);
  
  if (totalStop != 0) {
    vector<int> sendBuffer2(max<CFuint>(sizeSendBuffer2,1));
    vector<int> recvBuffer2(max<CFuint>(sizeSendBuffer2,1));
    
    sendBuffer2.clear();
    for(CFuint i=0; i<_nP; ++i){
      const CFuint addedSize = returningPartitionBeams[i].size();
      for(CFuint j=0; j<addedSize; ++j){
        sendBuffer2.push_back(returningPartitionBeams[i][j]);
      }
    }
    
    // build counts vectors
    for(CFuint p=0; p<_nP; p++){
      m_sendCounts[p] = returningPartitionBeams[p].size();
    }
    
    MPI_Alltoall(&m_sendCounts[0], 1, MPI_UNSIGNED, &m_recvCounts[0], 1, MPI_UNSIGNED, _comm);
    
    // build displacements vectors
    m_recvDispls[0] = 0;
    m_sendDispls[0] = 0;
    if(_nP>1){
      for(CFuint p=1; p<_nP; p++){
        m_recvDispls[p] = m_recvDispls[p-1] + m_recvCounts[p-1];
        m_sendDispls[p] = m_sendDispls[p-1] + m_sendCounts[p-1];
      }
    }
    
    const CFuint LastDisplacement2 = m_recvDispls[_nP-1] + m_recvCounts[_nP-1];
    
    // build receive buffers
    if (LastDisplacement2 > recvBuffer2.size()) {
      recvBuffer2.resize(LastDisplacement2);
    }
    assert(sendBuffer2.capacity()>0);
    assert(recvBuffer2.capacity()>0);
    
    // total exchange rays data
    MPI_Alltoallv(&sendBuffer2[0], &m_sendCounts[0], &m_sendDispls[0], MPI_INT, &recvBuffer2[0], &m_recvCounts[0], &m_recvDispls[0], MPI_INT, _comm);
    
    // clear content of partitionBeams in order to use it again
    for (CFuint p = 0; p < returningPartitionBeams.size(); ++p) {
      returningPartitionBeams[p].clear();
      cf_assert(returningPartitionBeams[p].size() == 0);
    }
    
    // store data
    const CFuint ld25 = LastDisplacement2/5; // blocks of 5-entries received
    for(CFuint i=0; i < ld25; ++i){
      
      // photon absorbed by a cell (emittingEntityType == -1)
      if (recvBuffer2[5*i+1] == -1) {
        const CFuint index = _EmittingCellIDs.find(recvBuffer2[5*i]);
        const CFuint index2 = index*(_NRAYS + 1) + m_cellCompletionIdx[index];
        _AbsorbingCellIDs  [index2] = recvBuffer2[5*i+2];
        _AbsorbingCellRanks[index2] = recvBuffer2[5*i+4];
        cf_assert(_AbsorbingCellRanks[index2] < (CFint)_nP);
        _AbsorbingCellType [index2] = recvBuffer2[5*i+3];
        m_cellCompletionIdx[index]++;
      }
      
      if (recvBuffer2[5*i+1] >= 0) {
        // photon absorbed by a wallFace (emittingEntityType >= 0)
        const CFuint index = _EmittingWallFaceIDs.find(recvBuffer2[5*i+1]);
        const CFuint index2 = index*(_NRAYS + 1) + m_faceCompletionIdx[index];
        _AbsorbingWallFaceIDs  [index2] = recvBuffer2[5*i+2];
        _AbsorbingWallFaceRanks[index2] = recvBuffer2[5*i+4]; // proc that sent the data
        cf_assert(_AbsorbingWallFaceRanks[index2] < (CFint)_nP);
        _AbsorbingWallFaceType [index2] = recvBuffer2[5*i+3];
        m_faceCompletionIdx[index]++;
      }
    }
  }
}

/////////////////////////////////////////////////////////////////////////////

void RadiativeTransferMonteCarlo::MonteCarlo()
{
  // initialize random number generator
  //time_t seconds;
  //seconds = time(NULL);
  //m_MTRand.seed(seconds);

  /*
  ofstream debug("hemiDirections.txt");
  vector<CFreal> normals(3);
  vector<CFreal> directions(3);
  normals[0]=0.;
  normals[1]=0.;
  normals[2]=1.;
  for(int i=0;i<10000;++i){
    directions=m_rand.hemiDirections(3,normals);
    debug<<directions[0]<<' '<<directions[1]<<' '<<directions[2]<<endl;
  }
  debug.close();
  */

  // vector storing beams that cross the partition faces
  vector<vector<Ray> > partitionBeams(_nP);
  // compute rays emitted by this partition domain //
  computeCellRays(partitionBeams);

  // WALL CELL RAYS ARE TURNED OFF!!
  // compute rays emitted by this partition domain //
  computeWallFaceRays(partitionBeams);
  // create vector of final total exchange
  vector<vector<CFuint> > returningPartitionBeams(_nP);

  // compute rays emitted by others partition domains
  computePartitionFaceRays(partitionBeams, returningPartitionBeams);
  // final comunication operation
  computeFinalCommunication(returningPartitionBeams);

}
RealVector RadiativeTransferMonteCarlo::convertCylindricalToCartesian(CFreal z,CFreal r, CFreal tetha){
  RealVector cartesian(3);
  cartesian[0]=r*cos(tetha);
  cartesian[1]=r*sin(tetha);
  cartesian[2]=z;
  return cartesian;
}
/////////////////////////////////////////////////////////////////////////////
RealVector RadiativeTransferMonteCarlo::covert2DAxito3Dcart(Ray ray, RealVector axiCoord){
  ray.tt=(axiCoord[0]-ray.startPoint[2])/ray.direction[2];
  return get3Dcoord(ray);
}


void RadiativeTransferMonteCarlo::convert3Dto2DAxi(Ray ray, RealVector axiCoord){
  axiCoord[0]=ray.startPoint[2]+ray.direction[2]*ray.tt; //z
  axiCoord[1]=sqrt(pow(ray.startPoint[0]+ray.direction[0]*ray.tt,2)+pow(ray.startPoint[1]+ray.direction[1]*ray.tt,2));//r
}

RealVector RadiativeTransferMonteCarlo::get3Dcoord(Ray ray){ //{x,y,z}
  RealVector cart(3);
  for (CFuint i=0; i<3;++i){
    cart[i]=ray.startPoint[i]+ray.direction[i]*ray.tt;
  }
  return cart;
}

/////////////////////////////////////////////////////////////////////////////
RealVector RadiativeTransferMonteCarlo::get2DAxicoord(Ray ray){ //{z,r}
  RealVector cart(3),axi(2);
  cart=get3Dcoord(ray); //{x,y,z}
  axi[0]=cart[2];
  axi[1]=sqrt( pow(cart[0],2.)+pow(cart[1],2.) );
  //cout<<"ray.t: "<<ray.tt<<" z: "<<r[0]<<" r: "<<r[1]<<endl;
  return axi;
}

/////////////////////////////////////////////////////////////////////////////
/*
RealVector RadiativeTransferMonteCarlo::nextAxiPoint(const RealVector previous2Daxi, const Ray ray){
  RealVector fNew(2),fOld(2),fMed(2),fMedAprox(2),fIni(2);
  CFreal curvatureRatio, threshold=1e-6, maxStep=0.5;
  static CFreal step=previous2Daxi[1]/3;
  bool coarsing=false,refining=false;
  if (ray.actualKS==0){step=previous2Daxi[1]/3;}
  //cout<<"on the nextAxiPoint routine"<<endl;

  //get previous t
  const CFreal t=(previous2Daxi[0]-ray.startPoint[2])/ray.direction[2]; // t=(z-z0)/dz
  //cout<<"direction 3D= "<<ray.direction[0]<<' '<<ray.direction[1]<<' '<<ray.direction[2]<<endl;
  //cout<<"initial t= "<<t<<endl;
  convert3Dto2DAxi(t,ray,fIni);
  convert3Dto2DAxi(t+step,ray,fNew);
  convert3Dto2DAxi(t+step/2,ray,fMed);
  while ( true ){
    //cout<<"current step: "<<step<<endl;
    for (CFuint i=0; i<2; ++i)
      fMedAprox[i]=( fIni[i]+fNew[i] )/2;
    curvatureRatio=MathFunctions::getDistance(fMed,fMedAprox);///MathFunctions::getDistance(fIni,fNew);
    //cout<<"error: "<<error<<endl;
    if (curvatureRatio>threshold){
      if (coarsing){
        //cout<<"New 2D AXi (fOld): "<<fOld[0]<<' '<<fOld[1]<<endl;
        //cout<<"New t= "<<(t+step/2)<<endl;
        return fOld;
      }
      step/=2;
      fNew=fMed;
      convert3Dto2DAxi(t+step/2,ray,fMed);
      refining=true;
    }
    else{
      if (refining || step>maxStep){
        //cout<<"New 2D AXi (fNew): "<<fNew[0]<<' '<<fNew[1]<<endl;
        //cout<<"New t= "<<(t+step)<<endl;
        return fNew;
      }
      step*=2;
      fOld=fNew;
      fMed=fNew;
      convert3Dto2DAxi(t+step,ray,fNew);
      coarsing=true;
    }
  }
}
*/
/////////////////////////////////////////////////////////////////////////////

CFreal RadiativeTransferMonteCarlo::get3Ddist(RealVector axi1, RealVector axi2, Ray ray){

  const CFreal t2=(axi2[0]-ray.startPoint[2])/ray.direction[2];
  const CFreal dt=abs(t2-(axi1[0]-ray.startPoint[2])/ray.direction[2]);
  //cout<<"dt= "<<dt<<endl;
  CFreal norm=0.;
  for(CFuint i=0;i<3;++i){
    norm+=pow(ray.direction[i]*dt,2.);
  }
  return sqrt(norm);
}
/////////////////////////////////////////////////////////////////////////////

CFuint RadiativeTransferMonteCarlo::rayTracing(Ray& beam,
                                               CFuint begin,
                                               EntityType& entity,
                                               const string& trs)
{  
  CFuint endEntityID = 0;
  //CFreal previousK = 0.;
  RealVector start2DAxi(2),end2DAxi(2);
  DataHandle<CFreal> normals = socket_normals.getDataHandle();

  //stringstream ss;
  //ss<<"intersections_"<<_myP<<".txt";
  //string ss2=ss.str();
  //ofstream debug(ss2.c_str(),ios::app);

//  debug<<endl<<"-1 "<<"-1 "<<endl;

  //cout<<endl<<"RAYTRACING "<<"ID: "<<beam.startEntityID<<endl;

  if(begin==0){ // ==>> the startPoint coincide with the barycentre of the cell with startEntityID as ID;
    //      used for rays that are emitted by a cell
    //cout<<"Start Point in the Baricentre: ";
    //cout<<"BEGUIN=0"<<endl;

    CFLog(DEBUG_MAX, "RadiativeTransferMonteCarlo::rayTracing() => begin==0\n");
    _particleTracking.setupStartPoint(beam.startEntityID);
    const RealVector& start = _particleTracking.getStartPoint();
    for (CFuint d = 0; d < _dim; ++d) {
      beam.startPoint[d] = start[d];
    }
    if (_Axi) {
      //Take the [z,r] coordinates, pick tetha and convert to a 3D cartesian point
      CFreal pi=  MathConsts::CFrealPi();
      CFreal tetha=m_rand.uniformRand(-pi,pi);
      //cout<<"tetha= "<<tetha<<endl;
      for (CFuint dd=0;dd<2;++dd){
        start2DAxi[dd]=beam.startPoint[dd];
      }
      //cout<<"StartPoint3D: "<<endl;
      RealVector cart3D(3);
      cart3D=convertCylindricalToCartesian(beam.startPoint[0],beam.startPoint[1],tetha);
      for (CFuint dd=0; dd<3;++dd){
        beam.startPoint[dd]=cart3D[dd];
      }
    }
    //previousK = getK(beam.startEntityID);
  }

  if(begin==1){ // ==>> the startPoint is into the cell with startEntityID as ID but does not coincide with the barycentre of this cell;
    //      used for rays that are reflected by a wall or that start from a partition face
    CFLog(DEBUG_MAX, "RadiativeTransferMonteCarlo::rayTracing() => begin==1\n");
    //cout<<"BEGUIN=1"<<endl;
    if(!_Axi){
      _particleTracking.setupStartPoint(beam.startEntityID, beam.startPoint);
    }

    //previousK = getK(beam.startEntityID);
    //for (CFuint d = 0; d < _dim; ++d) {
    //  debug<<beam.startPoint[d]<<' ';
    //  cout<<start[d]<<' ';
    //}

  }
  
  if(begin==2){ // ==>> the startPoint coincide with the barycentre of the face with startEntityID as ID that belonging to "trs";
    //      used for rays that are emitted by a wall face
    CFLog(DEBUG_MAX, "RadiativeTransferMonteCarlo::rayTracing() => begin==2 with TRS "<< trs <<"\n");
    _particleTracking.setupStartPoint(beam.startEntityID, trs);
    const RealVector& start = _particleTracking.getStartPoint();
    for (CFuint d = 0; d < _dim; ++d) {
      beam.startPoint[d] = start[d];
//      debug<<beam.startPoint[d]<<' ';
    }

    if (_Axi) {
      //Create the 3D Cylindric coordinates and convert to 3D cartesian

      CFreal pi=  MathConsts::CFrealPi();

      CFreal tetha=m_rand.uniformRand(-pi,pi);
      for (CFuint dd=0;dd<2;++dd){
        start2DAxi[dd]=beam.startPoint[dd];
      }
      RealVector cart3D=convertCylindricalToCartesian(beam.startPoint[0],beam.startPoint[1],tetha);
      for (CFuint dd=0; dd<3;++dd){
        beam.startPoint[dd]=cart3D[dd];
      }
      //cout<<"StartPoint3D: "<<beam.startPoint[0]<<' '<<beam.startPoint[1]<<' '<<beam.startPoint[2]<<endl;

      //Create the normal vector

      CFuint faceID=_particleTracking.getWallFaceID();
      CFuint startID=faceID*_dim;
      vector<CFreal> faceNormal(3);
      CFreal nz = -1 * normals[startID    ];     // Normals point outwards
      CFreal nr = -1 * normals[startID + 1];     // so we reverse them
      //const CFreal D=1/sqrt( nz*nz+nr*nr );
      //nr*=D;
      //nz*=D;

      //cout<<"Start Point: "<<start2DAxi[0]<<' '<<start2DAxi[1]<<endl;
      //cout<<"axi normals: "<<nz<<' '<<nr<<endl;
      CFreal tetha2=atan2(beam.startPoint[1],beam.startPoint[0]);
      faceNormal[0] =  nr * cos(tetha2);
      faceNormal[1] =  nr * sin(tetha2);
      faceNormal[2] =  nz;
      //*******************************************************
      // WARNING!!
      //
      // Normals geneated with Gmsh2CFmesh point inwards!
      //
      //*******************************************************
      CFreal invNorm=1/sqrt(pow(faceNormal[0],2) +pow(faceNormal[1],2) +pow(faceNormal[2],2) );
      for(CFuint i=0;i<3; ++i){
        faceNormal[i]*=invNorm;
      }
//      cout<<"normal face: "<<faceNormal[0]<<' '<<faceNormal[1]<<' '<<faceNormal[2]<<endl;

      //Generate Directions

      vector<CFreal> direction(3);
      m_rand.hemiDirections(3,faceNormal,direction);
      for(CFint i=0; i<3; ++i){
        beam.direction[i]=direction[i];
      }

      //cout<<"Direction: "<<beam.direction[0]<<' '<<beam.direction[1]<<' '<<beam.direction[2]<<endl;
      beam.tt=0.;


    }
    //previousK = getWallK(_particleTracking.getWallFaceID());

    //cout<<"actual cell ID= "<<actualCellID<<endl;

  }

//  debug<<endl;

  if(!_Axi){
  _particleTracking.setupPath(beam.direction, ParticleTracking::DIRECTION);
  }

  
  /* compute the path. Note that: 1) startCell is the first cell of the computed path
   *                              2) actualCell is the step by step computed cell along the path
   *                              3) endEntityID is the ID of the entity (cell, wallFace or partitionFace) where the ray is absorbed or goes through
   *                              4) exitFace is the face between the actualCell and the next cell which the ray goes through (=newCell)
   */
  CFuint actualCellID;
  if(!_Axi){
    actualCellID=_particleTracking.getStartCellID();
    _actualPoint = _particleTracking.getStartPoint();
    _intersectionPointOld = _particleTracking.getStartPoint();
  }
  else{
    if(beam.emittingEntityType>=0 && begin==2){ //if the ray starts from a face
      actualCellID=beam.emittingEntityType;
    }
    else{
      actualCellID=beam.startEntityID;
    }
    _actualPoint=get2DAxicoord(beam);
    _intersectionPointOld=get2DAxicoord(beam);

  }
  CFuint nbCrossedCells = 0;
  CFuint nbIter = 0;

  CFLog(DEBUG_MAX, "RadiativeTransferMonteCarlo::rayTracing() => actualCellID = " << actualCellID << "\n");

  bool foundEntity = false;
  //cout<<"Start K= "<<previousK<<endl;
  //cout<<"Beam KS= "<<beam.KS<<endl;
  while(nbIter <= _maxVisitedCells){
    nbIter++;
    CFint newCellID,exitFaceID;
    GeoEntOut out;
    if(_Axi){
      //cout<<"beam.tt before: "<<beam.tt<<endl;
      out=_particleTracking.myAxiRayTracing(beam,actualCellID);
      //cout<<"beam.tt after: "<<beam.tt<<endl;
      exitFaceID=out.exitFaceID;
      newCellID=out.exitCellID;

//      if (begin==2 ){
//        cout<<"exitFaceID= "<<exitFaceID<<" newCellID= "<<newCellID<<endl;
//      }
    }
    else{
      // global cellID
      newCellID = _particleTracking.tracking(actualCellID);
      exitFaceID = _particleTracking.getExitFaceID();
    }
    
    CFLog(DEBUG_MAX, "RadiativeTransferMonteCarlo::rayTracing() => newCellID = " << newCellID << "\n");
    // local (in this processor) face ID

    CFLog(DEBUG_MAX, "RadiativeTransferMonteCarlo::rayTracing() => exitFaceID = " << exitFaceID << "\n");
    /*
    if (exitFaceID<0 && _Axi){ //if the end point is inside the given cell in the axi case
       //cout<<"New EndPoint"<<endl;
       CFreal cellK=CELLK;
       CFreal dist=get3Ddist(_intersectionPointOld,end2DAxi,beam);
       beam.actualKS+=dist*cellK;
       //cout<<"new K: "<<beam.actualKS<<endl;
       start2DAxi=end2DAxi;
       CFreal arrEnd2DAxi[2];
       for(CFuint d=0; d<2;d++)
         arrEnd2DAxi[d]=end2DAxi[d];
       _particleTracking.setupStartPoint(actualCellID,arrEnd2DAxi);
       end2DAxi=nextAxiPoint(end2DAxi,beam);
       for(CFuint d=0; d<2;d++)
         arrEnd2DAxi[d]=end2DAxi[d];
       _particleTracking.setupPath(arrEnd2DAxi,ParticleTracking::ENDPOINT);
       _intersectionPoint=start2DAxi;
       _intersectionPointOld=_intersectionPoint;
       //debug<<_intersectionPoint[0]<<' '<<_intersectionPoint[1]<<endl;
       _internalPoint = _particleTracking.getInternalPoint();
       //cout<<"Internal Point: "<<_internalPoint[0]<<' '<<_internalPoint[1]<<' '<<endl;
    }
    */
    // proceed only if exitFaceID >= 0
    if(exitFaceID>=0){

      //cout<<"current K= "<<beam.actualKS<<endl;
      //cout<<"Arived to cell with ID= "<<newCellID<<" and Face ID= "<<exitFaceID<<endl;
      
      CFuint reflexCounter = 0;

      // compute optical lenght
      if(!_Axi){
      _internalPoint = _particleTracking.getInternalPoint();
      _intersectionPoint = _particleTracking.getIntersectionPoint();
      }
      else{
       _internalPoint=get2DAxicoord(beam);
       _intersectionPoint=get2DAxicoord(beam);
      }


      CFreal sK = MathFunctions::getDistance(_intersectionPoint, _intersectionPointOld);
      if (_Axi){
        sK=get3Ddist(_intersectionPoint, _intersectionPointOld, beam);
        /*
        const CFreal r = _intersectionPointOld[YY];
        const CFreal R = _intersectionPoint[YY];
        // SANNA CHANGE
        //const CFreal gamma = atan2(beam.direction[2], beam.direction[1]);
        const CFreal gamma = atan(beam.direction[2]/beam.direction[1]);
        const CFreal yz = sqrt(R*R - r*r*sin(gamma)*sin(gamma)) - r*cos(gamma);
        const CFreal x = abs(_intersectionPoint[XX] - _intersectionPointOld[XX]);
        sK = sqrt(yz*yz + x*x);
        */
      }
      CFreal cellK=CELLK;
      if (_testcaseID==0){
        cellK=getK(newCellID,beam);
      }

      //cout<<"distance: "<<sK<<endl;
      beam.actualKS += sK*cellK;
      //cout<<"updated K before checking absortion:"<<beam.actualKS<<endl;
      //      if(_Axi){
      //        CFreal pi=3.1415;
      //        beam.actualKS= cellK*sK*(abs(pi*_intersectionPointOld[1]+pi*_intersectionPoint[1])+1);
      //      }
      _actualPoint = _internalPoint;
      _intersectionPointOld = _intersectionPoint;
      //cout<<"New K= "<<beam.actualKS<<endl;
      if(beam.actualKS >= beam.KS){ // photon absorbed by a cell
        //if(begin==1){
        //  cout<<"ABSORVED BY CELL!!"<<endl;
        //}

//        for (CFuint i=0;i<_dim;++i){
//          debug<<_internalPoint[i]<<' ';
//        }
//        debug<<endl;

        //static CFuint nbAbsortions=0;
        //++nbAbsortons;

        endEntityID = actualCellID;
        entity = INTERNAL_CELL;

        //cout<<" ray tracing:  beam.actualKS = "<<beam.actualKS<<"\t  beam.KS = "<<beam.KS<<"\t  beam.startEntityID = "
       //     <<beam.startEntityID<<"\t  endEntityID = "<<endEntityID<<"\t  nbCrossedCells = "<<nbCrossedCells<<endl; //print informations to check the tracing
       // cout<<"absorved by entity with id "<<endEntityID<<endl;
        //cout<<"number of Absortions= "<<nbAbsortions;
        foundEntity = true;
        break;
      }


      // compute the destiny of the photon
      const bool isSym = binary_search(_symmetryIDs.begin(),_symmetryIDs.end(),exitFaceID);
      const bool isWall = binary_search(_wallIDs.begin(),_wallIDs.end(),exitFaceID);
      if (isWall || isSym) {
        //if(begin==1){
        //  cout<<"is wall or sym!"<<endl;
        //}
        //cout<<"beam.starPOINT: "<<beam.startPoint[0]<<' '<<beam.startPoint[1]<<' '<<beam.startPoint[2]<<endl;
        //cout<<"intersection point: ";
        //for (CFuint i=0;i<_dim;++i){
        //          cout<<_intersectionPoint[i]<<' ';
        //}
        //cout<<endl;

        if (isWall) {
          m_isWallFaceAbsorbing[exitFaceID] = true;
        }

        // the photon reaches the wall
        const CFreal wallK = getWallK(exitFaceID);
        //m_MTRand.seed(exitFaceID);
        //if (begin==2)
        //  cout<<_myP<<" wall!! "<<beam. emittingEntityType<<endl;
        // this should be one for Symmetry
        const CFreal reflectionProbability = (isWall) ? m_rand.uniformRand() : 1.; //M_rand.uniformRand();
        if(reflectionProbability >= wallK){ // the photon is reflected by the wall
          reflexCounter++;

          //cout<<"Reflected with \nIncoming direction: "
          //   <<beam.direction[0]<<' '<<beam.direction[1]<<' '<<beam.direction[2]<<endl;
          const CFuint startID = exitFaceID*_dim;
          CFreal nx,ny,nz;
          if(!_Axi){
            nx = - normals[startID    ];
            ny = - normals[startID + 1];
            nz = (_dim == 2) ? 0. : - normals[startID + 2];
          }
          else{
                   nz = - normals[startID    ];
            CFreal nr = - normals[startID + 1];
            //CFreal invNorm=1/sqrt(nz*nz+nr*nr);
            //nz*=invNorm;
            //nr*=invNorm;
            const RealVector pos=get3Dcoord(beam);
            CFreal tetha2=atan2(pos[1],pos[0]);
            nx = cos(tetha2) * nr;
            ny = sin(tetha2) * nr;
            //nz= nz;

            //cout<<"axi normals: "<<nz<<' '<<nr<<endl;

          }

          const CFreal invFaceLength = 1./sqrt(nx*nx + ny*ny + nz*nz);
          nx *= invFaceLength;
          ny *= invFaceLength;
          nz *= invFaceLength;
          CFreal normal[3]={nx,ny,nz};
          CFreal projection=0.0;
          CFuint dim2=(_Axi? 3 : _dim);
          for (CFuint i = 0; i < dim2; i++)
            projection += normal[i]*beam.direction[i];

           //const CFreal invReflDir = 1./_reflectedDirection.norm2();

          //beam.startEntityID = actualCellID; // take the cell in the fluid domain, not the wall one (newCell)
          //_particleTracking.setupStartPoint(beam.startEntityID, beam.startPoint);

          //cout<<"beam.startPoint= "<<beam.startPoint[0]<<' '<<beam.startPoint[1]<<' '<<beam.startPoint[2]<<endl;
          //cout<<"beam.direction= " <<beam.direction[0] <<' '<<beam.direction[1] <<' '<<beam.direction[2] <<endl;

          if(!_Axi){
            for(CFuint i = 0; i<_dim; i++ )
              beam.direction[i]=beam.direction[i]-2*projection*normal[i];

            for(CFuint s=0; s<_dim; ++s){
              beam.startPoint[s] = _intersectionPoint[s];
            }
            _particleTracking.setupStartPoint(actualCellID, beam.startPoint);
            _particleTracking.setupPath(beam.direction, ParticleTracking::DIRECTION);

            //_intersectionPoint=_particleTracking.getIntersectionPoint();
            //_internalPoint = _particleTracking.getInternalPoint();
            _intersectionPointOld=_intersectionPoint;
            _actualPoint=_intersectionPoint;

          }
          else{
            //cout<<"beam.tt= "<<beam.tt<<endl;
            RealVector startPoint=get3Dcoord(beam);
            for (CFuint i=0; i<3; ++i){
              beam.startPoint[i]=startPoint[i];
            }
            //_actualPoint=get2DAxicoord(beam);
            //_intersectionPoint=get2DAxicoord(beam);
            _intersectionPointOld=_intersectionPoint;
            // if (begin==1)
            //cout<<"intercetion point3D: "<<startPoint[0]<<' '<<startPoint[1]<<' '<<startPoint[2]<<endl;
            beam.tt=0.;

            for(CFuint i = 0; i<3; i++ )
              beam.direction[i]=beam.direction[i]-2*projection*normal[i];
            }
          //if (begin==1){
          //cout<<"Face normal: "<<nx<<' '<<ny<<' '<<nz<<endl;
          //cout<<"Outcoming direction: "<<beam.direction[0]<<' '<<beam.direction[1]<<' '<<beam.direction[2]<<endl;
          //cout<<"***********************************************"<<endl;
          //}
          //previousK = getK(beam.startEntityID);

          //if(begin==1){
          //  cout<<"END reflected!"<<endl;
          //}
          continue;
        }
        else{ // the photon is absorbed by the wall
          //cout<<"ABSORVED!"<<endl;
          endEntityID = exitFaceID;
          entity = WALL_FACE;
          foundEntity = true;
          break;
        }
      }

      if(binary_search(_boundaryIDs.begin(),_boundaryIDs.end(),exitFaceID)){ // photon disappears in a boundary
        //if (begin==1)
        //cout<<"Disappeared in boundary "<<exitFaceID<<endl;
//        for (CFuint i=0;i<_dim;++i){
//          debug<<_intersectionPoint[i]<<' ';
//        }
//        debug<<endl;
        endEntityID = 0;
        entity = DISAPPEARED;
        foundEntity = true;
        break;
      }

      if(binary_search(_partitionIDs.begin(),_partitionIDs.end(),exitFaceID)){ // photon goes in another partition
        CFLog(DEBUG_MAX, "RadiativeTransferMonteCarlo::rayTracing()  => PARTITION exitFaceID = " << exitFaceID << "\n");

        if(!_Axi){
          endEntityID = _particleTracking.getExitFaceGlobalGhostID();

          // SANNA CHANGE
          const RealVector& start = _intersectionPointOld;
          //const RealVector& start = _particleTracking.getInternalPoint();
          for (CFuint d = 0; d < _dim; ++d) {
            beam.startPoint[d] = start[d];
          }
          beam.startEntityID = endEntityID;
        }
        else{
          //const RealVector start = get3Dcoord(beam);
          for (CFuint d = 0; d < 3; ++d) {
          }
          beam.startEntityID = _particleTracking.getExitFaceGlobalGhostID(out);
        }


        entity = PARTITION_GHOSTCELL;
        foundEntity = true;
        break;

      }

      //previousK = getK(newCellID);

      actualCellID = newCellID;
      nbCrossedCells++;
    }
    else{
      //cout<<"enter negligle"<<endl;
      entity = NEGLIGIBLE;
      return 0;
    }
  }

  //cout<<" nbCrossedCells = "<<nbCrossedCells<<endl;
//  debug.close();
//  if (begin==2){
//    cout<<_myP<<" END raytracing ID"<<beam.emittingEntityType<<endl;
//  }
  return endEntityID;
}

/////////////////////////////////////////////////////////////////////////////

void RadiativeTransferMonteCarlo::build_MPI_strut_Ray(MPI_Datatype* MPI_Ray_ptr)
{ 

  // Ray instance from which we get the displacements
  Ray ray;
  
  // number of elements in each "block" of the new type
  int block_lengths[10];

  // displacement of each element from start of new type
  MPI_Aint displacements[10];

  // MPI types of the elements
  MPI_Datatype typelist[10];

  // use for calculating displacements
  MPI_Aint start_address;
  MPI_Aint address;

  // set block lengths
  block_lengths[0] = 3;
  block_lengths[1] = 3;
  block_lengths[2] = 1;
  block_lengths[3] = 1;
  block_lengths[4] = 1;
  block_lengths[5] = 1;
  block_lengths[6] = 1;
  block_lengths[7] = 1;
  block_lengths[8] = 1;
  block_lengths[9] = 1;
  
  // set typelist
  typelist[0] = MPI_DOUBLE;
  typelist[1] = MPI_DOUBLE;
  typelist[2] = MPI_DOUBLE;
  typelist[3] = MPI_DOUBLE;
  typelist[4] = MPI_UNSIGNED;
  typelist[5] = MPI_UNSIGNED;
  typelist[6] = MPI_UNSIGNED;
  typelist[7] = MPI_UNSIGNED;
  typelist[8] = MPI_DOUBLE;
  typelist[9] = MPI_DOUBLE;
  
  // first element is at displacement 0
  displacements[0] = 0;

  // calculate other displacements relative to the first element
  MPI_Address(&ray.direction[0], &start_address);
  
  // find address of Ray.startPonit[0] and displacement from Ray.directio[0]
  MPI_Address(&ray.startPoint[0], &address);
  displacements[1] = address - start_address;
  
  // find address of Ray.KS and displacement from Ray.directio[0]
  MPI_Address(&ray.KS, &address);
  displacements[2] = address - start_address;
  
  // find address of Ray.actualKS and displacement from Ray.directio[0]
  MPI_Address(&ray.actualKS, &address);
  displacements[3] = address - start_address;
  
  // find address of Ray.startEntityID and displacement from Ray.directio[0]
  MPI_Address(&ray.startEntityID, &address);
  displacements[4] = address - start_address;

  // find address of Ray.emittingEntityID and displacement from Ray.directio[0]
  MPI_Address(&ray.emittingEntityID, &address);
  displacements[5] = address - start_address;
  
  // find address of Ray.emittingEntityTypa and displacement from Ray.directio[0]
  MPI_Address(&ray.emittingEntityType, &address);
  displacements[6] = address - start_address;
  
  // find address of Ray.emittingProcessorRank and displacement from Ray.directio[0]
  MPI_Address(&ray.emittingProcessorRank, &address);
  displacements[7] = address - start_address;

  // find address of Ray.tt and displacement from Ray.directio[0]
  MPI_Address(&ray.tt, &address);
  displacements[8] = address - start_address;

  // find address of Ray.wavelength and displacement from Ray.directio[0]
  MPI_Address(&ray.wavelength, &address);
  displacements[9] = address - start_address;
  
  // build the derived datatype
  MPI_Type_struct(10, block_lengths, displacements, typelist, MPI_Ray_ptr);
  
  // commit new datatype
  MPIError::getInstance().check("MPI_Type_commit",
                                "RadiativeTransferMonteCarlo::build_MPI_strut_Ray()",
                                MPI_Type_commit(MPI_Ray_ptr));
}

/////////////////////////////////////////////////////////////////////////////

void RadiativeTransferMonteCarlo::computeHeatFlux()
{
  CFLog(VERBOSE, "RadiativeTransferMonteCarlo computeHeatFlux()\n");

//  cout<<endl<<"COMPUTE HEAT FLUX"<<endl;
//  cout<<"delta wave: "<<m_deltaWavs[m_spectralIdx]<<endl;
  
  DataHandle<RealVector> nstates = socket_nstates.getDataHandle();
  
  if (_testcaseID > 0) {cf_always_assert(m_spectralIdx == 0);}
  
  // create sendBuffer vector
  vector<vector<Qout> > sendBuffer3build(_nP);
  
  // compute leaving radiation heat flux for each cell
  CellTrsGeoBuilder::GeoData& cellData = m_cellBuilder.getDataGE();
  const CFuint nCells = cellData.trs->getLocalNbGeoEnts();
  DataHandle<CFreal> volumes = socket_volumes.getDataHandle();
  
  for(CFuint c=0; c<nCells; ++c){

    cellData.idx = c;
    GeometricEntity *const cell = m_cellBuilder.buildGE();
    CFuint localID = cell->getState(0)->getLocalID();
    CFuint nBP = m_overlapCellRanks[localID].size();
    CFreal NBP = double(nBP); //CHECK
    if(nBP == 0 ) NBP = 1.0;
    //cout<<"Cell with ID: "<<cell->getState(0)->getGlobalID()<<" has "<<nBP<<" overlaps"<<endl;
    const CFuint cellID = cell->getID();
    State *const innerState = cell->getState(0);
    CFreal Temperature = (*innerState)[m_temperatureID];
    
    if( (_testcaseID == 1) || (_testcaseID == 3) ){
      Temperature = 10000;
    }

    if(_testcaseID == 2){ // slab testcase
      const CFuint nNodes = cell->nbNodes();
      const vector<Node*>& nodes = *cell->getNodes();
      CFreal height = 0.0;
      for(CFuint f=0; f<nNodes; ++f){
        height += (*(nodes[f]))[YY];
      }
      height /= nNodes;
      Temperature = 10000.*height;
    }

    //CFreal Temperature = (*innerState)[_dim + 1];
    //CFreal r = _CellIDmap.find(cellID);
    //CFreal Temperature = (15000 - 5000*r);
    
    CFreal Volume =volumes[cellID];
    if(_Axi){
      //the "volume" extracted is the area of the cell
      const CFreal pi= MathConsts::CFrealPi();
      const RealVector centroid = cell->computeCentroid();
      Volume*= 2.*pi*centroid[YY];
      // store volumes
      m_axiVolumes[cellID] = Volume;
      //Volume=volumes[cellID]; //get area
    }
    
    const CFreal sigma = PhysicalConsts::StephanBolzmann();
    CFreal kCell = CELLK;//getEm(cell->getState(0)->getGlobalID()); // [W/m^3/str/m]

    // 4 comes from the integration in the solid angle
    // kCell is the linear absortion Coefficient [m^-1]
    CFreal Q_Out=0.;
      Q_Out = 4.*sigma*kCell*pow(Temperature,4.)*Volume; // [W]
      if (_testcaseID == 0) {
        cf_assert(m_spectralIdx < m_deltaWavs.size());
        CFreal I_out=getEm( cell->getState(0)->getGlobalID() ); // spectral radiative intensity
        Q_Out =4.*MathConsts::CFrealPi()*I_out*Volume;
        //cout<<"Iout: "<<I_out<<endl;
      }
//    cout<<"Cell number: "<<c<<" With Q_Out= "<<Q_Out<<endl;
    _RHSGas[c] += Q_Out;
    
    const CFreal e = (_NRAYS == 0) ? 0.0 : (Q_Out/(_NRAYS*NBP)); // energy carried by each photon
    // const CFreal e = Q_Out/(_NRAYS*NBP); // energy carried by each photon
    const CFuint startIdx = c*(_NRAYS + 1);
    
    for(CFuint i=startIdx; (_AbsorbingCellIDs[i] != -1); ++i){
      Qout qout;
      qout.ID = _AbsorbingCellIDs[i];
      qout.Type = _AbsorbingCellType[i];
      qout.Q = e;
      if(_AbsorbingCellRanks[i] < (CFint)_nP){
        sendBuffer3build[_AbsorbingCellRanks[i]].push_back(qout);
      }
    }

    m_cellBuilder.releaseGE();
  }
  
  // compute leaving radiation heat flux for each wall face
  FaceTrsGeoBuilder::GeoData& WallFacesData = m_wallFaceBuilder.getDataGE();
  const CFuint nbWallTrs = _wallTrsNames.size();
  DataHandle<CFreal> surfaces = socket_faceAreas.getDataHandle();
  CFuint ff = 0;
  for(CFuint j=0; j<nbWallTrs; ++j){
    SafePtr<TopologicalRegionSet> WallFaces = MeshDataStack::getActive()->getTrs(_wallTrsNames[j]);
    WallFacesData.trs = WallFaces;
    const CFuint nbFacesWall = WallFaces->getLocalNbGeoEnts();
    for(CFuint f=0; f<nbFacesWall; ++f){
      WallFacesData.idx = f;
      WallFacesData.isBFace = true;
      GeometricEntity *const face = m_wallFaceBuilder.buildGE();
      const CFuint faceID = face->getID();
      State *const innerState = face->getState(0);
      const CFuint localID = innerState->getLocalID();
      cf_assert(localID < m_overlapCellRanks.size());
      const CFuint nBP = m_overlapCellRanks[localID].size();
      CFreal NBP = double(nBP);
      if(nBP == 0) NBP = 1.0;
      
      // here you need Twall
      const CFuint nbNodesInFace = face->nbNodes();
      CFreal Temperature = 0.;// 0.5*((*innerState)[m_temperatureID] + (*ghostState)[m_temperatureID]);
      for(CFuint n=0; n<nbNodesInFace; ++n){
        Temperature += nstates[face->getNode(n)->getLocalID()][m_temperatureID];
      }
      Temperature *= static_cast<CFreal>(1./nbNodesInFace);
      //cout<<"temperature: "<<Temperature<<endl;
      
      CFreal Surface = surfaces[faceID];
      if(_Axi){
        CFuint nNodes = face->nbNodes();
        const vector<Node*>& nodes = *face->getNodes();
        CFreal Ybaricentre = 0.0;
        for(CFuint j=0; j<nNodes; ++j){
          Ybaricentre += (*(nodes[j]))[YY];
        }
        Ybaricentre /= nNodes;
        Surface *= Ybaricentre*2.0*MathConsts::CFrealPi();
      }
      cf_assert(ff < m_wallSurfaces.size());
      m_wallSurfaces[ff] = Surface;
      
      CFreal Q_Out = 0.;
      if (m_usePlanck) {
        //const CFreal lambda = getCurrentLambda()*1e-10;
        //CFreal intPlank = 0.;intPlankSimpson(Temperature,lambda,m_deltaWavs[m_spectralIdx]); //WRONG
        //const CFreal Em = getWallEm(faceID);          // this is the surface emissivity, not the spectral emissivity power
        Q_Out = 0.;//Em*intPlank*MathConsts::CFrealPi()*Surface;
        _RHSWall[ff] += Q_Out;
      }
      else {
        const CFreal sigma = PhysicalConsts::StephanBolzmann();
        const CFreal Em = getWallEm(faceID);          // this is the surface emissivity, not the spectral emissivity power
        Q_Out = sigma*Em*pow(Temperature,4.)*Surface; // this is dimensionally correct
        _RHSWall[ff] = Q_Out;
      }
      
      const CFreal e = (_NRAYS == 0) ? 0.0 : (Q_Out/(_NRAYS*NBP)); // energy carried by each photon
      //cout<<"photon energy: "<<e<<endl;
      //   const CFreal e = Q_Out/(_NRAYS*NBP); // energy carried by each photon
      const CFuint startIdx = ff*(_NRAYS + 1);
      for(CFuint i=startIdx; (_AbsorbingWallFaceIDs[i] != -1); ++i){
        Qout qout;
        qout.ID   = _AbsorbingWallFaceIDs [i];
        qout.Type = _AbsorbingWallFaceType[i];
        qout.Q    = e;
        sendBuffer3build[ _AbsorbingWallFaceRanks[i] ].push_back(qout);
      }
      m_wallFaceBuilder.releaseGE();
      ff++;
    }
  }
  
  // total exchange of leaving radiation heat flux
  
  // create MPI ray structure
  MPI_Datatype MPI_Qout_ptr;
  build_MPI_strut_Qout(&MPI_Qout_ptr);
  
  // build send buffers
  CFuint sizeSendBuffer3 = 0;
  for(CFuint i=0; i<_nP; ++i){
    const CFuint addedSize = sendBuffer3build[i].size();
    m_sendCounts[i] = addedSize;
    sizeSendBuffer3 += addedSize;
  }
  
  // build send buffers and counts vectors
  vector<Qout> sendBuffer3;
  vector<Qout> recvBuffer3;
  sendBuffer3.reserve(max<CFuint>(sizeSendBuffer3, 1));
  recvBuffer3.reserve(max<CFuint>(sizeSendBuffer3, 1));
  
  sendBuffer3.clear();
  for(CFuint i=0; i<_nP; ++i){
    const CFuint addedSize = sendBuffer3build[i].size();
    for(CFuint j=0; j<addedSize; ++j){
      sendBuffer3.push_back(sendBuffer3build[i][j]);
    }
  }
  
  // set and verify stop condition
  // if no one partition has to send beam to any other (it means also that all photons has been absorbed or dispersed)
  CFuint totalStop = 0;
  CFuint localStop = (sizeSendBuffer3 == 0) ? 0 : 1;
  MPI_Allreduce(&localStop, &totalStop, 1, MPI_UNSIGNED, MPI_MAX, _comm);

  if (totalStop != 0) {
    cf_assert(m_sendCounts.size() > 0);
    cf_assert(m_recvCounts.size() > 0);
    MPI_Alltoall(&m_sendCounts[0], 1, MPI_UNSIGNED, &m_recvCounts[0], 1, MPI_UNSIGNED, _comm);
    
    // build displacements vectors
    m_recvDispls[0] = 0;
    m_sendDispls[0] = 0;
    if(_nP>1){
      for(CFuint p=1; p<_nP; p++){
        m_recvDispls[p] = m_recvDispls[p-1] + m_recvCounts[p-1];
        m_sendDispls[p] = m_sendDispls[p-1] + m_sendCounts[p-1];
      }
    }
    
    const CFuint LastDisplacement3 = m_recvDispls[_nP-1] + m_recvCounts[_nP-1];
    
    // build receive buffers
    if (LastDisplacement3 > recvBuffer3.size()) {
      recvBuffer3.resize(LastDisplacement3);
    }
    
    // total exchange rays data
    MPI_Alltoallv(&sendBuffer3[0], &m_sendCounts[0], &m_sendDispls[0], MPI_Qout_ptr,
                  &recvBuffer3[0], &m_recvCounts[0], &m_recvDispls[0], MPI_Qout_ptr, _comm);
    
    // compute radiation heat source
    for(CFuint i=0; i<LastDisplacement3; ++i){
      Qout& qout = recvBuffer3[i];
      if(qout.Type == -1){ // cell
        if(_EmittingCellIDs.exists(qout.ID)){
          CFuint idx = _EmittingCellIDs.find(qout.ID);
          cf_assert(idx < _RHSGas.size());
          cf_assert(m_spectralIdx < m_deltaWavs.size());
          _RHSGas[idx] -= qout.Q; // emitted - absorbed
        }
      }
      
      if(qout.Type >= 0){ //face // Type stores the globalID of the adiacent cell
        if(_EmittingWallFaceIDs.exists(qout.Type)){
          CFuint idx = _EmittingWallFaceIDs.find(qout.Type);
          cf_assert(idx < _RHSWall.size());
          cf_assert(m_spectralIdx < m_deltaWavs.size());
          _RHSWall[idx] -= qout.Q; // emitted - absorbed
        }
      }
    }
    
    // add the absorbed radiation power to the copies of the overlapped cells in other partitions
    vector<vector<Qout> > sendBuffer4build(_nP);
    for(CFuint i=0; i<LastDisplacement3; ++i){
      cf_assert(i < recvBuffer3.size());
      Qout& qout = recvBuffer3[i];
      CFuint globalID = qout.ID;
      const CFint type = qout.Type;
      
      // SANNA
      if(type >= 0) {
        globalID = type;
      }
      
      if(type == -1 || (type >= 0 && _EmittingWallFaceIDs.exists(type))) {
        // CFLog(INFO, "globalID = " << globalID << "\n");
        const CFuint localID = _CellG2LIDmap.find(globalID);
        //CFLog(INFO, "globalID = " << globalID << " FOUND \n");

        cf_assert(localID < m_overlapCellRanks.size());
        const CFuint nBelongingPartitions = m_overlapCellRanks[localID].size();
        if(nBelongingPartitions > 0){
          for(CFuint pr=0; pr<nBelongingPartitions; ++pr){
            cf_assert(localID < m_overlapCellRanks.size());
            const CFuint proc = m_overlapCellRanks[localID][pr];
            if(proc != _myP){
              cf_assert(proc < sendBuffer4build.size());
              sendBuffer4build[proc].push_back(qout);
            }
          }
        }
      }
    }

    CFuint sizeSendBuffer4 = 0;
    std::vector<int> m_sendCounts4(_nP);
    std::vector<int> m_recvCounts4(_nP);
    vector<Qout> sendBuffer4;
    vector<Qout> recvBuffer4;
    for(CFuint i=0; i<_nP; ++i){
      const CFuint addedSize4 = sendBuffer4build[i].size();
      m_sendCounts4[i] = addedSize4;
      sizeSendBuffer4 += addedSize4;
      for(CFuint j=0; j<addedSize4; ++j){
        sendBuffer4.push_back(sendBuffer4build[i][j]);
      }
    }
    
    MPI_Alltoall(&m_sendCounts4[0], 1, MPI_UNSIGNED, &m_recvCounts4[0], 1, MPI_UNSIGNED, _comm);
    
    std::vector<int> m_sendDispls4(_nP,0);
    std::vector<int> m_recvDispls4(_nP,0);
    m_recvDispls4[0] = 0;
    m_sendDispls4[0] = 0;
    if(_nP>1){
      for(CFuint p=1; p<_nP; p++){
        m_recvDispls4[p] = m_recvDispls4[p-1] + m_recvCounts4[p-1];
        m_sendDispls4[p] = m_sendDispls4[p-1] + m_sendCounts4[p-1];
      }
    }
    const CFuint LastDisplacement4 = m_recvDispls4[_nP-1] + m_recvCounts4[_nP-1];
    recvBuffer4.resize(LastDisplacement4);
    MPI_Alltoallv(&sendBuffer4[0], &m_sendCounts4[0], &m_sendDispls4[0], MPI_Qout_ptr,
                  &recvBuffer4[0], &m_recvCounts4[0], &m_recvDispls4[0], MPI_Qout_ptr, _comm);

    for(CFuint i=0; i<LastDisplacement4; ++i){
      Qout& qout = recvBuffer4[i];
      if(qout.Type == -1){ // cell
        if(_EmittingCellIDs.exists(qout.ID)){
          CFuint idx = _EmittingCellIDs.find(qout.ID);
          cf_assert(idx < _RHSGas.size());
          cf_assert(m_spectralIdx < m_deltaWavs.size());
          _RHSGas[idx] -= qout.Q; // emitted - absorbed
        }
      }
      
      if(qout.Type >= 0){ // wall
        if(_EmittingWallFaceIDs.exists(qout.Type)){
          const CFuint idx = _EmittingWallFaceIDs.find(qout.Type);
          cf_assert(idx < _RHSWall.size());
          cf_assert(m_spectralIdx < m_deltaWavs.size());
          _RHSWall[idx] -= qout.Q; // emitted - absorbed
        }
      }
    }
  }
    CFreal totalSum, mySum=0.;
    for(CFuint ii=0;ii<_RHSWall.size();++ii){
      mySum+=_RHSWall[ii];
    }

    MPI_Reduce(&mySum, &totalSum, 1, MPI_DOUBLE, MPI_SUM, 0, _comm);

    if(_myP==0){
     ofstream sRHF("spectralRadiativeFlux.txt",ios::app);
     CFreal wav=(m_usePlanck)? m_wavReduced(0,m_spectralIdx) : 0;
     sRHF<< wav <<' '<<totalSum<<endl;
     sRHF.close();
    }

}

/////////////////////////////////////////////////////////////////////////////

void RadiativeTransferMonteCarlo::build_MPI_strut_Qout(MPI_Datatype* MPI_Qout_ptr)
{ 

  // Qout instance from which we get the displacements
  Qout qout;

  // number of elements in each "block" of the new type
  int block_lengths[3];

  // displacement of each element from start of new type
  MPI_Aint displacements[3];

  // MPI types of the elements
  MPI_Datatype typelist[3];

  // use for calculating displacements
  MPI_Aint start_address;
  MPI_Aint address;

  // set block lengths
  block_lengths[0] = 1;
  block_lengths[1] = 1;
  block_lengths[2] = 1;
  
  // set typelist
  typelist[0] = MPI_UNSIGNED;
  typelist[1] = MPI_UNSIGNED;
  typelist[2] = MPI_DOUBLE;
  
  // first element is at displacement 0
  displacements[0] = 0;

  // calculate other displacements relative to the first element
  MPI_Address(&qout.ID, &start_address);

  // find address of Qout.Type and displacement from Qout.ID
  MPI_Address(&qout.Type, &address);
  displacements[1] = address - start_address;

  // find address of Qout.Q and displacement from Qout.ID
  MPI_Address(&qout.Q, &address);
  displacements[2] = address - start_address;

  // build the derived datatype
  MPI_Type_struct(3, block_lengths, displacements, typelist, MPI_Qout_ptr);

  // commit new datatype
  MPI_Type_commit(MPI_Qout_ptr);

}

/////////////////////////////////////////////////////////////////////////////

void RadiativeTransferMonteCarlo::buildCellIDradiusmap()
{ 
  CellTrsGeoBuilder::GeoData& cellData = m_cellBuilder.getDataGE();
  const CFuint nCells = cellData.trs->getLocalNbGeoEnts();
  RealVector centrePoint(_dim);
  
  for(CFuint i=0; i<nCells; ++i){
    cellData.idx = i;
    GeometricEntity* const cell = m_cellBuilder.buildGE();
    const CFuint cellID = cell->getState(0)->getGlobalID();
    _particleTracking.computeAverage(*cell->getNodes(), cell->nbNodes(), centrePoint);
    
    CFreal radius = 0.;
    if(_dim == 2){
      radius = centrePoint[0];
    }
    if(_dim == 3){
      radius = sqrt(centrePoint[0]*centrePoint[0] + centrePoint[1]*centrePoint[1]);
    }
    _CellIDmap.insert(cellID,radius);
    m_cellBuilder.releaseGE();
  }
  _CellIDmap.sortKeys();
}

//////////////////////////////////////////////////////////////////////////////

void RadiativeTransferMonteCarlo::getRHFcircularSections()
{ 
  CFreal T = 10000;
  CFreal sigma = PhysicalConsts::StephanBolzmann();
  CFreal NonD = sigma*pow(T,4.);

  DataHandle<CFreal> volumes = socket_volumes.getDataHandle();

  CellTrsGeoBuilder::GeoData& cellData = m_cellBuilder.getDataGE();
  const CFuint nCells = cellData.trs->getLocalNbGeoEnts();
  vector<CFreal> Rbv;
  vector<CFreal> Qv;
  vector<CFreal> CIDv;

  vector<CFreal> globalRbv;
  vector<CFreal> globalQv;
  vector<CFreal> globalCIDv;

  for(CFuint i=0; i<nCells; ++i){
    cellData.idx = i;
    GeometricEntity* const cell = m_cellBuilder.buildGE();
    CFuint CellID = cell->getState(0)->getGlobalID();
    CFuint nNodes = cell->nbNodes();
    const vector<Node*>& nodes = *cell->getNodes();

    CFreal Xbaricentre = 0.0;
    CFreal Ybaricentre = 0.0;
    CFreal Zbaricentre = 0.0;
    for(CFuint j=0; j<nNodes; ++j){
      Xbaricentre += (*(nodes[j]))[XX];
      Ybaricentre += (*(nodes[j]))[YY];
      Zbaricentre += (*(nodes[j]))[ZZ];
    }
    Xbaricentre /= nNodes;
    Ybaricentre /= nNodes;
    Zbaricentre /= nNodes;
    CFreal Rbaricentre = pow((pow(Xbaricentre,2.) + pow(Ybaricentre,2.)),0.5);
    CFreal Q = _RHSGas[i]/volumes[i]/NonD;
    if (cell->getState(0)->isParUpdatable()){
      Rbv.push_back(Rbaricentre);
      Qv.push_back(Q);
      CIDv.push_back(CellID);
    }
    //else{
    // cout<<"gostcell"<<endl;
    //}
    m_cellBuilder.releaseGE();

  }

  int nLC = Rbv.size();
  vector<int> nGC(_nP);
  MPI_Gather(&nLC, 1, MPI_INT, &nGC[0], 1, MPI_INT, 0, _comm);
  vector<int> displ(_nP,0);
  for(CFuint d=1; d<_nP; ++d){
    displ[d] = displ[d-1]+nGC[d-1];
  }
  int lastDisplacement = displ[_nP-1] + nGC[_nP-1];
  globalRbv.resize(lastDisplacement);
  MPI_Gatherv(&Rbv[0], Rbv.size(), MPI_DOUBLE, &globalRbv[0], &nGC[0], &displ[0], MPI_DOUBLE, 0, _comm);
  globalQv.resize(lastDisplacement);
  MPI_Gatherv(&Qv[0], Qv.size(), MPI_DOUBLE, &globalQv[0], &nGC[0], &displ[0], MPI_DOUBLE, 0, _comm);
  globalCIDv.resize(lastDisplacement);
  MPI_Gatherv(&CIDv[0], CIDv.size(), MPI_DOUBLE, &globalCIDv[0], &nGC[0], &displ[0], MPI_DOUBLE, 0, _comm);

  if(_myP == 0 ){

    ofstream debug("results_3D.txt");
    ofstream debug2("DivergenceOfHeatFlux.plt");
    debug2<<"TITLE=\"Compiled Heat Divergence Values\""<<endl
          <<"VARIABLES=\"Adimentional Radius\",\"Divergence of heat flux\""<<endl
          <<"ZONE T=\""<<globalQv.size()<<" cels, "<<_NRAYS<<" photons per cell\" "<<endl;


    for(CFuint i=0;i<globalRbv.size();++i){
      debug <<globalRbv[i]<<' '<<globalQv[i]<<' '<<globalCIDv[i]<<endl;
    }
    debug.close();

    CFreal rMax=*max_element(globalRbv.begin(),globalRbv.end());
    CFreal rMin=*min_element(globalRbv.begin(),globalRbv.end());
    CFreal step=(rMax-rMin)/double(m_nbSteps);
    vector<CFreal> boundaries(m_nbSteps),qq(m_nbSteps),r2(m_nbSteps),size_r2(m_nbSteps);
    for (CFuint i=0; i<m_nbSteps;++i){
      boundaries[i]=rMin+step*i;
    }

    for (CFuint i=0;i<globalRbv.size();++i){
      for (CFuint j=1;j<m_nbSteps;++j){
        if (globalRbv[i]<=boundaries[j] && globalRbv[i]>=boundaries[j-1]){
          r2[j] +=globalRbv[i];
          qq[j]+=globalQv[i];
          ++size_r2[j];
          break;
        }
      }
    }
    for (CFuint i=0;i<m_nbSteps;++i){
      if (size_r2[i]>0){
        r2[i]/=size_r2[i];
        qq[i]/=size_r2[i];
        debug2<<r2[i]<<' '<<qq[i]<<endl;
      }
    }
    debug2.close();
  }
}

//////////////////////////////////////////////////////////////////////////////

void RadiativeTransferMonteCarlo::getRHFaxyCylinder()
{
  CFreal T = 10000;
  CFreal sigma = PhysicalConsts::StephanBolzmann();
  CFreal NonD = sigma*pow(T,4.);

  CellTrsGeoBuilder::GeoData& cellData = m_cellBuilder.getDataGE();
  const CFuint nCells = cellData.trs->getLocalNbGeoEnts();
  vector<CFreal> Rbv;
  vector<CFreal> Qv;
  //vector<CFreal> CIDv;

  vector<CFreal> globalRbv;
  vector<CFreal> globalQv;
  //vector<CFreal> globalCIDv;

  for(CFuint i=0; i<nCells; ++i){
    cellData.idx = i;
    GeometricEntity* const cell = m_cellBuilder.buildGE();
    //CFuint CellID = cell->getState(0)->getGlobalID();
    CFuint nNodes = cell->nbNodes();
    const vector<Node*>& nodes = *cell->getNodes();

    CFreal Q=_RHSGas[i]/m_axiVolumes[i]/NonD;

    CFreal Zbaricentre = 0.0;
    CFreal Rbaricentre = 0.0;
    for(CFuint j=0; j<nNodes; ++j){
      Zbaricentre += (*(nodes[j]))[XX];
      Rbaricentre += (*(nodes[j]))[YY];
    }
    Zbaricentre /= nNodes;
    Rbaricentre /= nNodes;
    Rbv.push_back(Rbaricentre);
    Qv.push_back(Q);
    m_cellBuilder.releaseGE();
  }

  int nLC = Rbv.size();
  vector<int> nGC(_nP);
  MPI_Gather(&nLC, 1, MPI_INT, &nGC[0], 1, MPI_INT, 0, _comm);
  vector<int> displ(_nP,0);
  for(CFuint d=1; d<_nP; ++d){
    displ[d] = displ[d-1]+nGC[d-1];
  }
  int lastDisplacement = displ[_nP-1] + nGC[_nP-1];
  globalRbv.resize(lastDisplacement);
  MPI_Gatherv(&Rbv[0], Rbv.size(), MPI_DOUBLE, &globalRbv[0], &nGC[0], &displ[0], MPI_DOUBLE, 0, _comm);
  globalQv.resize(lastDisplacement);
  MPI_Gatherv(&Qv[0], Qv.size(), MPI_DOUBLE, &globalQv[0], &nGC[0], &displ[0], MPI_DOUBLE, 0, _comm);
  //globalCIDv.resize(lastDisplacement);
  //MPI_Gatherv(&CIDv[0], CIDv.size(), MPI_DOUBLE, &globalCIDv[0], &nGC[0], &displ[0], MPI_DOUBLE, 0, _comm);

  if(_myP == 0 ){

    ofstream debug("results_axi.txt");
    ofstream debug2("DivergenceOfHeatFlux.plt");
    debug2<<"TITLE=\"Compiled Heat Divergence Values\""<<endl
          <<"VARIABLES=\"Adimentional Radius\",\"Divergence of heat flux\""<<endl
          <<"ZONE T=\""<<globalQv.size()<<" cels, "<<_NRAYS<<" photons per cell\" "<<endl;

    for(CFuint i=0;i<globalRbv.size();++i){
      debug <<globalRbv[i]<<' '<<globalQv[i]<<endl;
    }
    debug.close();


    CFreal rMax=*max_element(globalRbv.begin(),globalRbv.end());
    CFreal rMin=*min_element(globalRbv.begin(),globalRbv.end());
    CFreal step=(rMax-rMin)/double(m_nbSteps);
    vector<CFreal> boundaries(m_nbSteps),qq(m_nbSteps),r2(m_nbSteps),size_r2(m_nbSteps);
    for (CFuint i=0; i<m_nbSteps;++i){
      boundaries[i]=rMin+step*i;
    }

    for (CFuint i=0;i<globalRbv.size();++i){
      for (CFuint j=1;j<m_nbSteps;++j){
        if (globalRbv[i]<=boundaries[j] && globalRbv[i]>=boundaries[j-1]){
          r2[j] +=globalRbv[i];
          qq[j]+=globalQv[i];
          ++size_r2[j];
          break;
        }
      }
    }
    for (CFuint i=0;i<m_nbSteps;++i){
      if(size_r2[i]>0){
        r2[i]/=size_r2[i];
        qq[i]/=size_r2[i];
        debug2<<r2[i]<<' '<<qq[i]<<endl;
      }
    }
    debug2.close();
  }

}

//////////////////////////////////////////////////////////////////////////////

void RadiativeTransferMonteCarlo::getRHFslab()
{   
  DataHandle<CFreal> volumes = socket_volumes.getDataHandle();
  CellTrsGeoBuilder::GeoData& cellData = m_cellBuilder.getDataGE();
  const CFuint nCells = cellData.trs->getLocalNbGeoEnts();
  vector<pair<CFreal, pair<CFreal, CFreal> > > PairsVector;
  PairsVector.reserve(nCells);
  for(CFuint i=0; i<nCells; ++i){
    cellData.idx = i;
    GeometricEntity* const cell = m_cellBuilder.buildGE();
    CFuint CellID = cell->getID();
    CFreal Volume = volumes[CellID];
    CFuint nNodes = cell->nbNodes();
    const vector<Node*>& nodes = *cell->getNodes();
    CFreal Xbaricentre = 0.0;
    CFreal Ybaricentre = 0.0;
    CFreal Zbaricentre = 0.0;
    for(CFuint j=0; j<nNodes; ++j){
      Xbaricentre += (*(nodes[j]))[XX];
      Ybaricentre += (*(nodes[j]))[YY];
      Zbaricentre += (*(nodes[j]))[ZZ];
    }
    Xbaricentre /= nNodes;
    Ybaricentre /= nNodes;
    Zbaricentre /= nNodes;
    //if((Xbaricentre>0.4) && (Xbaricentre<0.5)){
    CFreal Q = _RHSGas[i]/Volume;
    pair<CFreal, CFreal> Pair1(Ybaricentre,Q);
    pair<CFreal, pair<CFreal, CFreal> > Pair2(Zbaricentre, Pair1);
    PairsVector.push_back(Pair2);
    //}
    m_cellBuilder.releaseGE();
  }
  
  sort(PairsVector.begin(),PairsVector.end(), PairSort<CFreal, pair<CFreal, CFreal> >());
  
  CFuint counterZ = 1;
  CFuint counterY = 1;
  for(CFuint g=1; g<PairsVector.size(); ++g){
    if(MathTools::MathChecks::isEqualWithError(PairsVector[g].first, PairsVector[g-1].first, 1e-2)){
      counterY++;
    }
    else{
      counterZ++;
      counterY = 1;
    }
  }
  
  vector<pair<CFreal, CFreal> > Middle;
  Middle.reserve(counterY);
  CFuint middleIdx = counterY*ceil(counterZ/2);
  for(CFuint l=0; l<counterY; l++){
    Middle.push_back(PairsVector[middleIdx+l].second);
  }
  sort(Middle.begin(),Middle.end());

  CFreal h = 1.0; //0.01;
  //  cout<<" height of a slab; "<<endl;
  //cin>>h;
  
  vector<CFreal> height(counterY,0.0);
  vector<CFreal> RH(counterY,0.0);
  for(CFuint l=0; l<counterY; l++){
    height[l] = Middle[l].first;
    RH[l] = Middle[l].second;
  }
  vector<CFreal> RHF(counterY,0.0);
  RHF[0] = RH[0]*h;
  for(CFuint l=1; l<counterY; l++){
    RHF[l] = RH[l]*h + RHF[l-1];
  }

  vector<CFreal> RH2(counterY,0.0);
  for(CFuint l=1; l<counterY; l++){
    RH2[l] = (RHF[l] - RHF[l-1])/h;
  }
  RH2[0] = RH2[1];
  
  const string fileName = "results.txt" + StringOps::to_str(_myP);
  ofstream fout(fileName.c_str());
  
  if (_nP > 1) {fout << counterY << endl;}
  for(CFuint i=0; i<counterY; i++){
    cout<<"  height = "<<height[i]<<" \t "<<" Heat = "<<RH[i]<<endl;
    fout<<" "<<height[i]<<" \t "<<RH[i]<<endl;
  }
  fout.close();

  const string fileName2 = "results2.txt" + StringOps::to_str(_myP);
  ofstream fout2(fileName2.c_str());
  if (_nP > 1) {fout2 << counterY << endl;}
  for(CFuint i=0; i<counterY; i++){
    cout<<"  height = "<<height[i]<<" \t "<<" heat flux = "<<RHF[i]<<endl;
    fout2<<" "<<height[i]<<" \t "<<RHF[i]<<endl;
  }
  fout2.close();

  const string fileName3 = "results3.txt" + StringOps::to_str(_myP);
  ofstream fout3(fileName3.c_str());
  if (_nP > 1) {fout3 << counterY << endl;}
  for(CFuint i=0; i<counterY; i++){
    cout<<"  height = "<<height[i]<<" \t "<<" heat = "<<RH2[i]<<endl;
    fout3<<" "<<height[i]<<" \t "<<RH2[i]<<endl;
  }
  fout3.close();
  
  MPI_Barrier(_comm);
  
  // assembly of different output files
  if (_nP > 1 && _myP == 0) {
    const CFreal dy = 0.01; //std::abs(RHF[1] - RHF[0]);
    const CFreal nbIntervals = h/dy;
    vector<CFreal> output(static_cast<int>(nbIntervals+1), 0.);
    
    CFreal y = 0.;
    CFreal heat = 0;
    for (CFuint i = 0; i < _nP; ++i) {
      const string file = "results.txt" + StringOps::to_str(i);
      ifstream fin(file.c_str());
      CFuint nbEntries = 0;
      fin >> nbEntries;
      for (CFuint e = 0; e < nbEntries; ++e) {
        fin >> y >> heat;
        cf_assert(static_cast<int>(y/dy) < (int)output.size());
        output[static_cast<int>(y/dy)] += heat;
      }
    }
    
    ofstream outfile("results_final.txt");
    cout << "outfile.size() " << output.size() << endl;
    for (CFuint i = 0; i < output.size(); ++i) {
      outfile << dy*i + 0.005 << " " << output[i] << endl;
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void RadiativeTransferMonteCarlo::setFaceIDs(const vector<string>& names,
                                             vector<CFuint>& faceIDs,
                                             bool isWall)
{  
  // AL: the ordering of the wall TRS's must be consistent with the one in AeroForces
  FaceTrsGeoBuilder::GeoData& facesData = m_wallFaceBuilder.getDataGE();
  const CFuint nbWallTrs = names.size();
  for(CFuint j=0; j<nbWallTrs; ++j){
    SafePtr<TopologicalRegionSet> WallFaces = MeshDataStack::getActive()->getTrs(names[j]);
    facesData.trs = WallFaces;
    const CFuint nbFacesWall = WallFaces->getLocalNbGeoEnts();
    faceIDs.reserve(faceIDs.size() + nbFacesWall);
    for(CFuint i=0; i<nbFacesWall; ++i){
      facesData.idx = i;
      GeometricEntity *const face = m_wallFaceBuilder.buildGE();
      faceIDs.push_back(face->getID());
      
      if (isWall) {
        m_mapWallFaceIDtoCellID.insert(face->getID(), face->getState(0)->getGlobalID());
      }
      
      m_wallFaceBuilder.releaseGE();
    }
  }
  sort(faceIDs.begin(), faceIDs.end());
  
  if (isWall) {
    m_mapWallFaceIDtoCellID.sortKeys();
  }
}


/////////////////////////////////////////////////////////////////////////////

void RadiativeTransferMonteCarlo::myReduceSpectra(const RealMatrix& indata,
                                                RealMatrix& outdata)
{
  CFLog(INFO, "RadiativeTransferMonteCarlo::myReduceSpectra() => start\n");

  cf_always_assert(indata.size() > 0);
  cf_always_assert(outdata.size() > 0);
  cf_always_assert(outdata.nbRows() == indata.nbRows());

  // here implement the spectra reduction
  const CFreal minWav     = m_radLibrary->getMinWavelength();
  const CFreal maxWav     = m_radLibrary->getMaxWavelength();
  const CFreal nInStride  = m_radLibrary->getWavelengthStride();
  const CFreal nOutStride = m_reducedSpectraSize;

  m_deltaWav = 1e-10;

  vector<CFreal> inWav(nInStride);
  vector<CFreal> inEm(nInStride);
  vector<CFreal> inAm(nInStride);

  const CFreal dWav= (maxWav-minWav ) / nOutStride;

  for (CFuint s = 0; s < outdata.nbRows(); ++s) {

     //extract the vectors from the indata Matrix
     for(CFuint i=0; i<nInStride; ++i){
        inWav[i] = indata(s,i*3  );
        inEm [i] = indata(s,i*3+1);
        inAm [i] = indata(s,i*3+2);
     }

     vector<CFreal> tempOutWav,tempOutAm,tempOutEm;

     for( CFuint outIdx= 0, inIdx=0; outIdx<nOutStride; ++outIdx){
        CFreal outWavMin = minWav    + dWav*outIdx;
        CFreal outWavMax = outWavMin + dWav;

       if(inIdx>0){
          //add first point at outWavMin
          CFreal startWavEm=linearInterpol(inWav[inIdx],inEm[inIdx],inWav[inIdx+1],inEm[inIdx+1],outWavMin);
          CFreal startWavAm=linearInterpol(inWav[inIdx],inAm[inIdx],inWav[inIdx+1],inAm[inIdx+1],outWavMin);

          tempOutWav.push_back( outWavMin );
          tempOutEm.push_back( startWavEm );
          tempOutAm.push_back( startWavAm );
          if(inWav[inIdx+1]<outWavMax){
            ++inIdx;
          }
        }
        while(inWav[inIdx+1]<outWavMax && inIdx<inWav.size()){
          //add middle points
          tempOutWav.push_back( inWav[inIdx] );
          tempOutAm.push_back(  inAm[inIdx]  );
          tempOutEm.push_back(  inEm[inIdx]  );
          ++inIdx;
       }
       if(inIdx<inWav.size() ){
         //add last point at outWavMax
         CFreal endWavEm=linearInterpol(inWav[inIdx],inEm[inIdx],inWav[inIdx+1],inEm[inIdx+1],outWavMax);
         CFreal endWavAm=linearInterpol(inWav[inIdx],inAm[inIdx],inWav[inIdx+1],inAm[inIdx+1],outWavMax);

         tempOutWav.push_back( outWavMax );
         tempOutEm.push_back( endWavEm );
         tempOutAm.push_back( endWavAm );
       }

       //calculate the integral and add it in the reduced coeff vectors
       m_wavReduced(s,outIdx) = ( tempOutWav.back()+tempOutWav.front() )/2.;
       CFreal sumEm=0.,sumAm=0.;
       for(CFuint i=0; i<tempOutWav.size()-1; ++i){
          sumEm+= ( tempOutEm[i+1]+tempOutEm[i] )/2.*( tempOutWav[i+1]-tempOutWav[i] );
          sumAm+= ( tempOutAm[i+1]+tempOutAm[i] )/2.*( tempOutWav[i+1]-tempOutWav[i] );
       }

       m_emReduced(s,outIdx) = sumEm*m_deltaWav;
       m_amReduced(s,outIdx) = sumAm/( tempOutWav.back()-tempOutWav.front() );

       //clear the temporary vectors for new calculation
       tempOutWav.clear();
       tempOutEm.clear();
       tempOutAm.clear();
     }

     if(_myP==0 && s==0){
       CFreal sum=0;
       for(CFuint y=0;y<m_emReduced.nbCols();++y){
         //cout<<m_emReduced(s,y)<<' ';
         sum+=m_emReduced(s,y);
       }
       cout<<endl<<"IntSum="<<sum<<endl;
     }

  }

  CFLog(INFO, "RadiativeTransferMonteCarlo::myReduceSpectra() => END\n");
}

/////////////////////////////////////////////////////////////////////////////

void RadiativeTransferMonteCarlo::reduceSpectra(const RealMatrix& indata, 
                                                RealMatrix& outdata)
{
  //outdata=indata;
  const CFreal minWav    = CFreal( m_radLibrary->getMinWavelength()    );
  const CFreal maxWav    = CFreal( m_radLibrary->getMaxWavelength()    );
  const CFreal inStride  = CFreal( m_radLibrary->getWavelengthStride() );
  const CFreal dInStride      = (maxWav-minWav)/inStride*1e-10;
  //return;
  CFLog(INFO, "RadiativeTransferMonteCarlo::reduceSpectra() => start\n");
  
  cf_always_assert(indata.size() > 0);
  cf_always_assert(outdata.size() > 0);
  cf_always_assert(outdata.nbRows() == indata.nbRows());
  
  // here implement the spectra reduction
  const CFuint inStride2  = m_radLibrary->getWavelengthStride();
  const CFuint outStride2 = m_reducedSpectraSize;
  const CFuint dWav = inStride2/(outStride2-2); // oustride has additional first and last point
  cf_assert(inStride2%(outStride2-2) == 0); // inStride must be a multiple of outStride
  
  if (_myP == 0) {
    ofstream fout1("indata.txt");
    fout1 << "TITLE = Original spectrum of radiative properties\n";
    fout1 << "VARIABLES = Wavl EmCoef AbCoef \n";
    for (CFuint i = 0; i < indata.nbCols()/3; ++i) {
      fout1 << indata(0,i*3) << " " << indata(0,i*3+1) << " " << indata(0,i*3+2) << endl;
    }
  }
  
  CFLog(VERBOSE, "RadiativeTransferMonteCarlo::reduceSpectra() => original nb wavelengths = " << outdata.nbCols()/3 << "\n");
  CFLog(VERBOSE, "RadiativeTransferMonteCarlo::reduceSpectra() => reduced nb wavelengths  = " << indata.nbCols()/3 << "\n");

  const CFreal ovdWav = 1./CFreal(dWav);
  for (CFuint s = 0; s < outdata.nbRows(); ++s) {
    for (CFuint i = 0; i < m_reducedSpectraSize; ++i) {
      if (i > 0 && i < m_reducedSpectraSize-1) {
        CFuint start = (i-1)*dWav*3;
        const CFuint index = i*3;
        //CFLog(DEBUD_MIN, "RadiativeTransferMonteCarlo::reduceSpectra() => " << index << " => [" << start << ", " << end << "]\n");

        outdata(s,index  ) = 0.;
        outdata(s,index+1) = 0.;
        outdata(s,index+2) = 0.;

        for (CFuint j = 0; j < dWav; ++j, start+=3) {
          cf_assert( index+2 < outdata.nbCols() );
          cf_assert( start+2 <  indata.nbCols() );

          outdata(s,index )  +=  indata(s,start);
          outdata(s,index+1) +=  indata(s,start+1);
          outdata(s,index+2) +=  indata(s,start+2);
        }

        //Average Wavelength:
        outdata(s,index)   *= ovdWav;
        //Average Emission coefficient
        //outdata(s,index+1) *= 1/(tempWav[dWav]-tempWav[0]);
        //Average Absorption coefficient
        //outdata(s,index+2) /= (tempWav[dWav]-tempWav[0]);
      }
      else if (i == 0)	{
        // first entry
        outdata(s,0) = indata(s,0);
        outdata(s,1) = indata(s,1);
        outdata(s,2) = indata(s,2);
      }
      else if (i == m_reducedSpectraSize-1) {
        // last entry
        const CFuint ostart = outdata.nbCols()-3;
        const CFuint istart = indata.nbCols()-3;
        cf_assert(ostart+2 < outdata.nbCols());
        cf_assert(istart+2 < indata.nbCols());

        outdata(s,ostart)   = indata(s,istart);
        outdata(s,ostart+1) = indata(s,istart+1);
        outdata(s,ostart+2) = indata(s,istart+2);
      }
    }
  }
  
  if (_myP == 0) {
    ofstream fout2("outdata.txt");
    fout2 << "TITLE = Reduced spectrum of radiative properties\n";
    fout2 << "VARIABLES = Wavl EmCoef AbCoef \n";
    for (CFuint i = 0; i < outdata.nbCols()/3; ++i) {
      fout2 << outdata(0,i*3) << " " << outdata(0,i*3+1) << " " << outdata(0,i*3+2) << endl;
    }
  }
  
  CFLog(INFO, "RadiativeTransferMonteCarlo::reduceSpectra() => Total Delta = ["
        << indata(0,0) << " - " << indata(0, indata.nbCols()-3) << "]\n");
 /*
  CFuint final = 0;
  for (CFuint i = 1; i < m_reducedSpectraSize-1; ++i) {
    CFuint start = (i-1)*dWav*3;
    CFLog(VERBOSE, "Delta [" << i-1 << "] = [" << indata(0,start) << " - ");
    m_deltaWavs[i-1] = 0.;
    cf_assert(i-1 < m_deltaWavs.size());
    for (CFuint j = 0; j < dWav; ++j, start+=3) {
      // average wavelength
      if (start+3 <= (indata.nbCols()-3)) {
        final = start+3;
      }
      else {
        break;
      }
      
      CFLog(DEBUG_MIN,  "(" << start << "-" << final << ")");
      m_deltaWavs[i-1] += indata(0,final) - indata(0,start);
    }
    CFLog(VERBOSE, indata(0,final) << "]\n");
  }
  for(CFuint ii=1;ii<)
  outdata(s,index+1)
  */
  m_deltaWavs = dInStride;
  
  CFLog(INFO, "RadiativeTransferMonteCarlo::reduceSpectra() => END\n");
}
//////////////////////////////////////////////////////////////////////////////

void RadiativeTransferMonteCarlo::flagOverlapCells()
{  
  CellTrsGeoBuilder::GeoData& cellData = m_cellBuilder.getDataGE();
  SafePtr<TopologicalRegionSet> cells = MeshDataStack::getActive()->getTrs("InnerCells");
  cellData.trs = cells;
  const CFuint nCells = cellData.trs->getLocalNbGeoEnts();
  
  // flag overlap elements
  m_isOverlapCell.resize(nCells, false);
  
  // suppose that just one space method is available
  SafePtr<SpaceMethod> spaceMethod = getMethodData().getCollaborator<SpaceMethod>();
  SafePtr<CellCenterFVM> fvmcc = spaceMethod.d_castTo<CellCenterFVM>();
  cf_assert(fvmcc.isNotNull());
  SafePtr<NodalStatesExtrapolator<CellCenterFVMData> > nodalExtrap = 
    fvmcc->getData()->getNodalStatesExtrapolator();
  
  // first layer of overlap
  for(CFuint c=0; c<nCells; ++c){
    cellData.idx = c;
    GeometricEntity *const cell = m_cellBuilder.buildGE();
    const State* const state = cell->getState(0);
    
    if (!state->isParUpdatable()) {
      m_isOverlapCell[state->getLocalID()] = true;
      
      // consider all vertex neighbors of the current ghost state
      const CFuint nbNodes = cell->nbNodes();
      for (CFuint iNode = 0; iNode < nbNodes; ++iNode) {
        const CFuint nodeID = cell->getNode(iNode)->getLocalID();
        const vector<State*>& stencil = nodalExtrap->getNodalStateNeighbors(nodeID);
        for (CFuint in = 0; in < stencil.size(); ++in) {
          // neglect ghost states which do not correspond to an actual cell
          if (!stencil[in]->isGhost()) {
            m_isOverlapCell[stencil[in]->getLocalID()] = true;
          }
        }
      }
    }
    
    m_cellBuilder.releaseGE();
  }
  
  // second layer of overlap (if second order reconstruction is used)
  DataHandle<State*, GLOBAL> states = socket_states.getDataHandle();
  if (MeshDataStack::getActive()->getNbOverlapLayers() == 2) {
    for(CFuint c=0; c<nCells; ++c) {
      if (m_isOverlapCell[c] && states[c]->isParUpdatable()) {
        cellData.idx = c;
        GeometricEntity *const cell = m_cellBuilder.buildGE();

        // consider all vertex neighbors of the current ghost state
        const CFuint nbNodes = cell->nbNodes();
        for (CFuint iNode = 0; iNode < nbNodes; ++iNode) {
          const CFuint nodeID = cell->getNode(iNode)->getLocalID();
          const vector<State*>& stencil = nodalExtrap->getNodalStateNeighbors(nodeID);
          for (CFuint in = 0; in < stencil.size(); ++in) {
            // neglect ghost states which do not correspond to an actual cell
            if (!stencil[in]->isGhost()) {
              m_isOverlapCell[stencil[in]->getLocalID()] = true;
            }
          }
        }

        m_cellBuilder.releaseGE();
      }
    }
  }
  
  
  
  // here every processor will broadcast its overlap global cell IDs to all
  // other processors so that each processor can keep accountancy of the
  // ranks that share the same IDs
  vector<CFuint> overlapCellIDs;
  // overlapCellIDs.reserve(nCells); // oversized memory
  
  m_overlapCellRanks.resize(nCells);
  
  CFLog(DEBUG_MIN, "OVERLAP => \n");
  CFuint overlapSize = 0;
  for (CFuint i = 0; i < m_isOverlapCell.size(); ++i) {
    if (m_isOverlapCell[i]) {
      CFLog(DEBUG_MIN, states[i]->getGlobalID() << " ");
      overlapCellIDs.push_back(states[i]->getGlobalID());
      overlapSize++;
      // consider the current rank
      m_overlapCellRanks[i].reserve(_nP); // preallocation of memory to be fast later
      m_overlapCellRanks[i].push_back(_myP);
    }
  }
  
  CFLog(DEBUG_MIN, "\n");
  sort(overlapCellIDs.begin(), overlapCellIDs.end());
  
  CFuint maxNbCells = 0;
  CFuint localNbCells = nCells;
  MPI_Allreduce(&localNbCells, &maxNbCells, 1, MPI_UNSIGNED, MPI_MAX, _comm);
  cf_assert(maxNbCells >= nCells);
  
  vector<CFuint> bCastOverlapCellIDs(maxNbCells+1); // oversized buffer
  
  for (CFuint root = 0; root < _nP; ++root) {
    if (root == _myP) {
      // first entry is the size of the overlap to be broadcast
      bCastOverlapCellIDs[0] = overlapSize;
      // copy local data to buffer
      for (CFuint k = 0;  k < overlapSize; ++k) {
        cf_assert(k < overlapCellIDs.size());
        cf_assert(k+1 < bCastOverlapCellIDs.size());
        bCastOverlapCellIDs[1+k] = overlapCellIDs[k];
      }
    }
    
    MPI_Bcast(&bCastOverlapCellIDs[0], bCastOverlapCellIDs.size(), MPI_UNSIGNED, root, _comm);
    
    if (root != _myP) {
      const CFuint nbC = bCastOverlapCellIDs[0]+1;
      for (CFuint c = 1; c < nbC; ++c) {
        cf_assert(c < bCastOverlapCellIDs.size());
        const CFuint globalID = bCastOverlapCellIDs[c];
        if (binary_search(overlapCellIDs.begin(), overlapCellIDs.end(), globalID)) {
          // consider the current rank
          const CFuint localCellID = _particleTracking.getLocalCellID(globalID);
          cf_assert(localCellID < m_overlapCellRanks.size());
          m_overlapCellRanks[localCellID].push_back(root);
        }
      }
    }
  }
  
  for (CFuint i = 0; i < m_overlapCellRanks.size(); ++i) {
    if (m_overlapCellRanks[i].size() > 0) {
      CFLog(DEBUG_MIN, "Overlap cellID [" << states[i]->getGlobalID() << "] in proc: ");
      for (CFuint p = 0; p < m_overlapCellRanks[i].size(); ++p) {
        CFLog(DEBUG_MIN, m_overlapCellRanks[i][p] << " ");
      }
      CFLog(DEBUG_MIN, "\n");
    }
  }
}

////////////////////////////////////////////////////////////////////////////// 

void RadiativeTransferMonteCarlo::buildCellIDmap()
{
  CellTrsGeoBuilder::GeoData& cellData = m_cellBuilder.getDataGE();
  const CFuint nCells = cellData.trs->getLocalNbGeoEnts();
  for(CFuint i=0; i<nCells; ++i){
    cellData.idx = i;
    GeometricEntity *const cell = m_cellBuilder.buildGE();
    _CellG2LIDmap.insert(cell->getState(0)->getGlobalID(), cell->getState(0)->getLocalID());
    m_cellBuilder.releaseGE();
  }
  _CellG2LIDmap.sortKeys();
}

//////////////////////////////////////////////////////////////////////////////

/*void RadiativeTransferMonteCarlo::fillWallHeatFluxes()
{
  // loop over the boundary faces and accumulate flux on the local wall face (check on isUpdatable()?)
  
  // 1) if you hit the inlet boundary you stop and assign to the RHSWall
  
  // 2) if you hit a partition
  //      you communicate to all processes sharing the overlap:
  //      - the L\R global stateIDs of the partition face
  //      - the current accumulated heat flux
  //      you identify the local face ID within the cell and you keep on marching in the opposite direction
  //      then back to 1) or 2) till completion
  }*/

//////////////////////////////////////////////////////////////////////////////

} // namespace RadiativeTransfer

} // namespace COOLFluiD






















