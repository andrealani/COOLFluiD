#ifndef COOLFluiD_RadiativeTransfer_RadiativeTransferMonteCarloHSNBHSNB_hh
#define COOLFluiD_RadiativeTransfer_RadiativeTransferMonteCarloHSNBHSNB_hh

//////////////////////////////////////////////////////////////////////////////

#include <numeric>
#include <boost/random.hpp>

#include "Framework/DataProcessingData.hh"
#include "Framework/FaceTrsGeoBuilder.hh"
#include "Config/ConfigObject.hh"
#include "Framework/ConvectiveVarSet.hh"
#include "Framework/DofDataHandleIterator.hh"
#include "RandomNumberGenerator.hh"
#include "LagrangianSolver/LagrangianSolver.hh"
#include "RadiativeTransfer/PostProcess/PostProcess.hh"
#include "RadiativeTransfer/RadiationLibrary/RadiationPhysicsHandler.hh"
#include "RadiativeTransfer/RadiationLibrary/Models/HSNB/HSNBRadiator.hh"
#include "LagrangianSolver/ParallelVector/ParallelVector.hh"
#include "RadiativeTransfer/Solvers/MonteCarlo/HSNBPhotonData.hh"
#include "RadiativeTransfer/Solvers/MonteCarlo/HSNBLocalParameterSet.hh"
#include "RadiativeTransfer/Solvers/MonteCarlo/HSNBParamSynchronizer.hh"

#include "RadiativeTransfer/Solvers/MonteCarlo/PhotonData.hh"


#include <utility>
#include <time.h>
#include <cmath>
#include <algorithm>
#include <boost/progress.hpp>
#include <boost/random.hpp>

#include "MathTools/MathChecks.hh"
#include "MathTools/MathConsts.hh"
#include "RadiativeTransfer/RadiativeTransfer.hh"
#include "Framework/DataProcessingData.hh"
#include "Framework/ConvectiveVarSet.hh"
#include "Framework/DataSocketSource.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/GeometricEntity.hh"
#include "Framework/GeometricEntityPool.hh"
#include "Framework/FaceTrsGeoBuilder.hh"
#include "Framework/CellTrsGeoBuilder.hh"
#include "Framework/PhysicalConsts.hh"
#include "Common/CFPrintContainer.hh"
#include "Common/MPI/MPIError.hh"
#include "FiniteVolume/CellCenterFVM.hh"
#include "MathTools/MathFunctions.hh"
#include "Common/MPI/MPIStructDef.hh"
#include "Framework/SocketBundleSetter.hh"
#include "LagrangianSolver/ParallelVector/ParallelVector.hh"
#include "Common/Stopwatch.hh"


//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {
    namespace FiniteVolume {
      class CellCenterFVMData;
    }
  }

 namespace RadiativeTransfer {

//////////////////////////////////////////////////////////////////////////////


typedef LagrangianSolver::Particle<HSNBPhotonData> Photon;

//////////////////////////////////////////////////////////////////////////////

template<class PARTICLE_TRACKING>
class RadiativeTransferMonteCarloHSNB : public Framework::DataProcessingCom
{
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
   static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor
   */
   RadiativeTransferMonteCarloHSNB(const std::string& name);

  /**
   * Default destructor
   */
   ~RadiativeTransferMonteCarloHSNB();

  /**
   * Configures the command.
   */
   void configure(Config::ConfigArgs& args);

  /**
   * PreProcess phase
   */
   void setup();

  /**
   * executeOnTrs
   */
   void executeOnTrs();

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
   std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

  /**
   * Returns the DataSocket's that this command provides as sources
   * @return a vector of SafePtr with the DataSockets
   */
   std::vector<Common::SafePtr<Framework::BaseDataSocketSource> > providesSockets();

private:

  /**
   * MonteCarlo
   */
   void MonteCarlo();

  /**
   * Traces a single photon until it is either fully absorbed or has left the computational domain of this
   * process. After every raytracing step the HSNBParameters are added to a HSNBPhotonTrace by calling the
   * corresponding function of the HSNBRadiator (addStateParams(...)) and partial absorption is computed in the
   * currentCell. If the photon passes over to the domain of another process its Trace is committed to the
   * HSNBParamSynchronizer.
   */
  CFuint rayTracing(Photon& photon);

  CFuint rayTracingFullAbsorption(Photon& beam);
  /**
   * build vector of radiative heat source along a single radius in the middle of the cilinder
   */
   void getRHF();

  /**
   * build vector of radiative heat source along a radius computed as average along a circular sections in the middle of the cilinder
   */
   void getRHFslab();

  /**
   * compute Heat Flux
   */
   void computeHeatFlux();

  /// Main raytracing loop is contained here:  We generate photons until all photons to be generated
  /// have been traced (toGenerateCellPhotons+toGenerateFacePhotons=0). If one of the conditions to force
  /// synchronization is met we synchronize / pass photons and HSNBPhotonTraces between processes by calling
  /// m_paramSynchronizer.synchronize() and m_lagrangianSolver.sincronize()
  void computePhotons();

  /// Check whether one of the different conditions to enforce sync of photons between comp. domains is met
  /// No more photons will be generated if one of the conditions is met and the rayTracing cycle will jump
  /// to the synchronization ONCE ALL PHOTONS RECEIVED FROM OTHER PROCESSES HAVE BEEN TRACED
  void checkSyncForceConditions();

  /// Generates a photon and associates it with a wavenumber, an energy fraction and a random direction
  /// The photon parameters are generated by the HSNBRadiator following a probability distribution
  /// for emission over all mechanisms.
  bool getCellPhotonData(Photon& ray);




private:

  LagrangianSolver::LagrangianSolver<HSNBPhotonData, PARTICLE_TRACKING> m_lagrangianSolver;

  /// Socket for the Gas Radiative Heat Source
  Framework::DataSocketSource < CFreal > socket_qrad; //GasRadiativeHeatSource

  /// the socket to the radiative heat flux at the wall faces
  Framework::DataSocketSource < CFreal > socket_qradFluxWall;

  /// handle to the face normals
  Framework::DataSocketSink< CFreal> socket_normals;

  /// storage of the nodes
  Framework::DataSocketSink < Framework::Node* , Framework::GLOBAL > socket_nodes;

  /// storage of the states
  Framework::DataSocketSink < Framework::State* , Framework::GLOBAL > socket_states;

  /// storage of the ghost states
  Framework::DataSocketSink<Framework::State*> socket_gstates;

  /// storage of the nodal state's
  Framework::DataSocketSink < RealVector > socket_nstates;

  /// socket for the time-averaged Joule heat source storage
  Framework::DataSocketSink<CFreal> socket_volumes;

  /// storage of the face areas
  Framework::DataSocketSink<CFreal> socket_faceAreas;

  /// rank of the processor containing the right state for partition faces
  Framework::DataSocketSink<CFuint> socket_rankPartitionFaces;

  /// IDs corresponding to the cell for which the normal point outward
  Framework::DataSocketSink<CFint> socket_isOutward;

  /// storage of face centroids
  Framework::DataSocketSink<CFreal> socket_faceCenters;

  /// pointer to the RadiationPhysicsHandler
  Common::SharedPtr<RadiationPhysicsHandler> m_radiation;

  /// temporary array for direction
  RealVector m_direction;

  /// temporary array for entry direction
  RealVector m_entryDirection;

  /// temporary array for exit direction
  RealVector m_exitDirection;

  /// temporary array for position
  RealVector m_position;

  /// temporary array for normal
  RealVector m_normal;

  /// temporary array for face normal in 3D
  RealVector m_faceNormal3;

  /// temporary array for cartesian position in 3D
  RealVector m_cartPosition3;

  /// temporary array for parametric coordinate in 3D
  RealVector m_sOut3;

  /// Current datatype with the correct array sizes for synchro
  MPIStruct HSNBUserdatatype;

  /// number of dimension
  CFuint m_dim;

  ///my process
  CFuint m_myProcessRank;

  ///number processes
  CFuint m_nbProcesses;

  CFreal m_tolerance;

  bool m_dynamicRayDistribution;

  SharedPtr<HSNBRadiator> m_HSNBRadiator;

  SafePtr<Radiator> m_baseRadiatorPtr;

  HSNBParamSynchronizer m_paramSynchronizer;

  ///MPI_COMM_WORLD
  MPI_Comm m_comm;

  /// cell builder
  Framework::GeometricEntityPool<Framework::CellTrsGeoBuilder> m_cellBuilder;

  /// wall face builder
  Framework::GeometricEntityPool<Framework::FaceTrsGeoBuilder> m_wallFaceBuilder;

  MathTools::CFMat<CFuint> m_wallTypes;

  Common::SelfRegistPtr< PostProcess > m_postProcess;

  /// inner TRS boundary names
  std::string m_radLibraryName;

  std::string m_postProcessName;

  /// number of rays that each element emits, ONLY USED IF m_dynamicRayDistribution=false!
  CFuint m_nbRaysElem;

  /// Maximum number of rays to be traced between synchronization of photon traces
  CFuint m_nbRaysCycle;

  /// total number of rays to be distributed among all cells. Cells with higher emissive power will receive a larger fraction
  CFuint m_nbRaysTotal;

  /// Flag is set to true if one of the conditions to jump to the synchronization is met. Is used to prevent excessive memory usage
  bool m_enforceSync;

  /// total number of cells crossed by all traces in this raytracing cycle
  CFuint m_nbCrossedCellsCycle;

  /// maximum number of visited cells for a single trace
  CFuint m_maxVisitedCellsTrace;

  /// maximum total number of visited cells for one raytracing cycle. One of the conditions to enforce sync.
  CFuint m_maxVisitedCellsCycle;

  /// maximum time in seconds allowed beforce a sync between procs is enforced
  CFreal m_maxSecondsBetweenSyncs;

  /// True if it is an axisymmetric simulation
  bool m_isAxi;

  RandomNumberGenerator m_rand;

  CFuint m_sendBufferSize;

  /// If this margin is exceeded no more traces will be commited and sync is enforced
  CFuint m_maxTraceBufferByteSize;

  /// Max number of traces crossing the domain boundary to be committed to the send buffer before sync is forced.
  CFuint m_maxTraceBufferTracesCommitted;

  /// Counts the number of traces committed in the current cycle
  CFuint m_nbTracesCommittedCycle;

  /// The total size allowed for all buffers used in trace parallelization in Byte
  ///
  /// The default value = 0 allows for unlimited buffer memory usage (if no other limitations
  /// such as maximum sendbuffer size etc. are imposed). Note that obeying this
  /// memory limitation will potentially force processes to idle until enough memory is free
  CFuint m_maxGlobalBufferByteSize;

  RealVector m_ghostStateInRadPowers;

  CFuint m_dim2;

  void getTotalEnergy();

  void printPhoton(const Photon &photon);

  void printPerformanceMetrics();

  CFreal m_totalRadPower;
  CFreal m_totalSubdomainRadPower;
  CFreal m_totalGlobalRadPower;

  /// Percentage of the total NB of rays held by all processes
  CFreal m_loadPercentages;
  /// Maximum NB of rays to be computed by a single process
  CFreal m_maxRaysSingleProcess;

  CFuint m_nbCrossedCellsTotal;

  CFreal m_crossedCellsAllPercentages;





  vector< CFuint> plotProcessIds;
  //TRACKING END

  LagrangianSolver::ParallelVector<CFreal> m_stateRadPower;
  LagrangianSolver::ParallelVector<CFreal> m_stateNetRadPower;
  LagrangianSolver::ParallelVector<CFreal> m_stateInRadPowers;

  CFuint m_nbDiatomics;
  CFuint m_nbAtomics;
  CFuint m_nbContinua;
  CFuint m_nbNonThickDiatomics;
  CFuint m_nbThickDiatomics;

  //Performance analysis metrics:
  CFreal m_avgTimePerPhoton;
  CFreal m_maxAvgTimePerPhoton;
  CFreal m_avgProcessEfficiency;
  CFreal m_avgProcessEfficiencyGlobal;
  CFuint m_avgCellsPerPhotonCycle;
  CFuint m_avgCellsPerPhotonTotal;
  CFreal m_avgTimePerPhotonGlobal;
  CFuint m_avgCellsPerPhotonGlobal;
  CFuint m_avgNewPhotonsTracedCycle;
  CFuint m_avgNewPhotonsTracedCycleGlobal;
  CFreal m_avgIdleTime;
  CFreal m_avgIdleTimeGlobal;
  CFuint m_avgPhotonsCycle;
  CFint m_nbPhotonsTracedCycle;
  CFreal m_avgCycleTime;
  bool m_trackTelemetry;
  CFuint m_nbCrossingForty;

  CFuint m_debugTraceUID;

  Stopwatch<WallTime> m_stpCycle;
  Stopwatch<WallTime> m_stpTracing;
  CFint m_nbIterations;
  CFreal m_emittedPhotonEnergyState;
  CFreal m_totalWallAbsorption;



  HSNBDataContainer m_curTraceSet;

  vector<CFuint> m_nbPhotonsState;

  vector<CFreal> m_ghostStateRadPower;
  vector<CFuint> m_nbPhotonsGhostState;

  vector<CFint> m_randomTrajectories;

  //FIXME: persistent memory for cell and face iterators
  CFuint  m_iphoton_cell_fix, m_istate_cell_fix;
  CFuint  m_iphoton_face_fix, m_igState_face_fix;
  CFuint m_istate_nextCell_fix;

  CFreal m_energyCheck;

  CFreal m_relaxationFactor;

  CFuint m_debugPhotonsComitted;

  CFint m_previousCellID;


  bool getFacePhotonData(Photon &ray);
}; // end of class RadiativeTransferMonteCarloHSNB

//////////////////////////////////////////////////////////////////////////////

template<class PARTICLE_TRACKING>
void RadiativeTransferMonteCarloHSNB<PARTICLE_TRACKING>::defineConfigOptions
(Config::OptionList& options)
{
  using namespace std;
  using namespace COOLFluiD::Common;

  
  options.addConfigOption< CFuint >("numberOfRays","number of rays sent by each element.");
  options.addConfigOption< CFuint >("MaxNbVisitedCellsTrace","Maximum number of visited cells for a single trace.");
  options.addConfigOption< bool >("Axi","True if it is an axisymmetric simulation.");

  options.addConfigOption< bool >("DynamicRayDistribution","Assign photons to cells proportional to the emissive power of the cell instead of using a fixed number of photons /cell");
  options.addConfigOption< bool >("TrackTelemetry","Performs additional synchronization to compute telemetry data (average cells per photon e.g...). Takes extra time.");

  options.addConfigOption< string >("PostProcessName","Name of the post process routine");
  options.addConfigOption< CFuint >("sendBufferSize","Size of the buffer for communication");
  options.addConfigOption< CFuint >("nbRaysCycle","Number of rays to emit before communication step");
  options.addConfigOption< CFreal >("relaxationFactor","Relaxation Factor");
  options.addConfigOption< CFreal >("MaxSecondsBetweenSyncs","The maximum allowable time between synchronization steps");

  options.addConfigOption< CFuint >("MaxGlobalBufferByteSize","The maximum allowable memory usage by all buffers used for trace parallelization (in total, in byte)");

  options.addConfigOption< CFuint >("totalNbRays","total number of rays to be distributed over all cells in accordance with the cells emissive power");
  options.addConfigOption< CFuint >("traceBufferMaxByteSize","Margin to enforce sync for photon traces to avoid excessive memory usage.");
  options.addConfigOption< CFuint >("MaxNbVisitedCellsCycle","If the total number of cells crossed in a raytracing cycle exceeds this margin sync is enforced after the current trace is finished.");

  //  options.addConfigOption< bool >("plotTrajectories","Photon trajectories will be exported to a tecplot geometry plot");
//  options.addConfigOption< CFuint >("MaxNbTrajectories","Maximum number of trajectories to be plotted");
//  options.addConfigOption< trajectoryExportType >("exportType", "Determines selection criterion for trajectory ids (\"Random\" or \"ConstantSpacing\"");

}

//////////////////////////////////////////////////////////////////////////////

template<class PARTICLE_TRACKING>
RadiativeTransferMonteCarloHSNB<PARTICLE_TRACKING>::RadiativeTransferMonteCarloHSNB
(const std::string& name):
  DataProcessingCom(name),
  m_lagrangianSolver(name),
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
  socket_faceCenters("faceCenters"),
  m_radiation(new RadiationPhysicsHandler("RadiationPhysicsHandler")),
  m_direction(),
  m_entryDirection(),
  m_exitDirection(),
  m_position(),
  m_normal(),
  m_faceNormal3(),
  m_cartPosition3(),
  m_sOut3()
{
  using namespace std;
  using namespace COOLFluiD::Framework;
  using namespace COOLFluiD::Common;

  addConfigOptionsTo(this);

  m_tolerance=0.001;





  m_postProcessName = "PostProcessNull";
  setParameter("PostProcessName",&m_postProcessName);

  m_nbRaysElem=300;
  setParameter("numberOfRays",&m_nbRaysElem);

  m_nbRaysTotal=40e6;
  setParameter("totalNbRays",&m_nbRaysTotal);

  m_maxVisitedCellsTrace = 10000;
  setParameter("MaxNbVisitedCellsTrace",&m_maxVisitedCellsTrace);

  m_maxVisitedCellsCycle = 100000;
  //m_maxVisitedCellsCycle = 100e6;
  setParameter("MaxNbVisitedCellsCycle",&m_maxVisitedCellsCycle);

  m_maxSecondsBetweenSyncs=60.0;
  setParameter("MaxSecondsBetweenSyncs",&m_maxSecondsBetweenSyncs);

  //Default 0 means no limitation
  m_maxGlobalBufferByteSize = 0;
  setParameter("MaxGlobalBufferByteSize",&m_maxGlobalBufferByteSize);

  m_dynamicRayDistribution=true;
  setParameter("DynamicRayDistribution",&m_dynamicRayDistribution);

  m_trackTelemetry=false;
  setParameter("TrackTelemetry",&m_trackTelemetry);

  m_isAxi = false;
  setParameter("Axi",&m_isAxi);

  //m_sendBufferSize=10000;
  m_sendBufferSize=100000;
  setParameter("sendBufferSize", &m_sendBufferSize);

  //100 MB
  m_maxTraceBufferByteSize=1e8;


  //2MB
  //m_traceBufferMaxByteSize=2e8;

  setParameter("traceBufferMaxByteSize", &m_maxTraceBufferByteSize);

  m_nbRaysCycle = m_sendBufferSize / 2;
  //m_nbRaysCycle = 50000;
  setParameter("nbRaysCycle", &m_nbRaysCycle);

  m_relaxationFactor = 1.;
  setParameter("relaxationFactor", &m_relaxationFactor);




}

/////////////////////////////////////////////////////////////////////////////

template<class PARTICLE_TRACKING>
RadiativeTransferMonteCarloHSNB<PARTICLE_TRACKING>::~RadiativeTransferMonteCarloHSNB()
{
}

/////////////////////////////////////////////////////////////////////////////

template<class PARTICLE_TRACKING>
void RadiativeTransferMonteCarloHSNB<PARTICLE_TRACKING>::configure(Config::ConfigArgs& args)
{
  CFAUTOTRACE;

  using namespace COOLFluiD::Framework;
  using namespace COOLFluiD::Common;

  DataProcessingCom::configure(args);
  cf_assert(m_radiation.isNotNull());
  configureNested ( m_radiation.getPtr(), args );

  m_postProcess = Environment::Factory< PostProcess >::getInstance().
    getProvider(m_postProcessName)->create(m_postProcessName);
  cf_assert(m_postProcess.isNotNull());
  configureNested ( m_postProcess.getPtr(), args );
}

//////////////////////////////////////////////////////////////////////////////

template<class PARTICLE_TRACKING>
vector<SafePtr<BaseDataSocketSink> > RadiativeTransferMonteCarloHSNB<PARTICLE_TRACKING>::
needsSockets()
{
  std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > result;

  result.push_back(&socket_nodes);
  result.push_back(&socket_states);
  result.push_back(&socket_nstates);
  result.push_back(&socket_gstates);
  result.push_back(&socket_normals);
  result.push_back(&socket_volumes);
  result.push_back(&socket_faceAreas);
  result.push_back(&socket_rankPartitionFaces);
  result.push_back(&socket_isOutward);
  result.push_back(&socket_faceCenters);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

template<class PARTICLE_TRACKING>
vector<SafePtr<BaseDataSocketSource> > RadiativeTransferMonteCarloHSNB<PARTICLE_TRACKING>::
providesSockets()
{
  std::vector<Common::SafePtr<Framework::BaseDataSocketSource> > result;
  result.push_back(&socket_qrad);
  result.push_back(&socket_qradFluxWall);
  return result;
}

/////////////////////////////////////////////////////////////////////////////

template<class PARTICLE_TRACKING>
void RadiativeTransferMonteCarloHSNB<PARTICLE_TRACKING>::setup()
{
  using namespace std;
  using namespace COOLFluiD::Framework;
  using namespace COOLFluiD::MathTools;
  using namespace COOLFluiD::Common;
  using namespace COOLFluiD::LagrangianSolver;

  const std::string nsp = getMethodData().getNamespace();

  // MPI parameters
  m_myProcessRank = PE::GetPE().GetRank(nsp);
  m_nbProcesses = PE::GetPE().GetProcessorCount(nsp);
  m_comm = PE::GetPE().GetCommunicator(nsp);
  MPIError::getInstance().init(m_comm, m_myProcessRank);

  // Limit the max traces per cycle to prevent excessive load balance
  m_maxTraceBufferTracesCommitted=m_nbRaysCycle / m_nbProcesses;

  // set dimensions
  m_dim = PhysicalModelStack::getActive()->getDim();
  m_dim2=(m_isAxi? 3 : m_dim);

  Framework::SocketBundle sockets;
  sockets.states      = socket_states;
  sockets.gstates     = socket_gstates;
  sockets.nodes       = socket_nodes;
  sockets.nstates     = socket_nstates;
  sockets.normals     = socket_normals;
  sockets.isOutward   = socket_isOutward;
  sockets.volumes     = socket_volumes;
  sockets.faceCenters = socket_faceCenters;
  sockets.faceAreas   = socket_faceAreas;

  m_direction.resize(m_dim2);
  m_entryDirection.resize(m_dim2);
  m_exitDirection.resize(m_dim2);
  m_position.resize(m_dim2);
  m_normal.resize(m_dim2);
  m_faceNormal3.resize(3);
  m_cartPosition3.resize(3);
  m_sOut3.resize(3);

  m_stateRadPower.setDataSockets(sockets);
  m_stateInRadPowers.setDataSockets(sockets);

  // copy the content of the "initialNodalStates" (if available) into "nstates"
  const string dhName = nsp + "_initialNodalStates";
  int localSize  = 0;
  int globalSize = 0;
  vector<CFreal> bufferVec;
  CFreal* buffer = CFNULL;

  if (MeshDataStack::getActive()->getDataStorage()->checkData(dhName)) {
   cf_assert(PE::GetPE().GetRank(nsp) == 0);
   DataHandle<CFreal> initialNodalStates =
      MeshDataStack::getActive()->getDataStorage()->getData<CFreal>(dhName);
    buffer = &initialNodalStates[0];
    localSize = (int) initialNodalStates.size();
  }

  MPIError::getInstance().check
    ("MPI_Allreduce", "RadiativeTransferMonteCarloHSNB::setup()",
      MPI_Allreduce(&localSize, &globalSize, 1, MPIStructDef::getMPIType(&localSize), MPI_MAX, m_comm));

  if (globalSize > 0) {
    if (buffer == CFNULL) {
     bufferVec.resize(globalSize);
     buffer = &bufferVec[0];
    }

    MPIError::getInstance().check
      ("MPI_Bcast", "RadiativeTransferMonteCarloHSNB::setup()",
        MPI_Bcast(buffer, globalSize, MPIStructDef::getMPIType(buffer), 0, m_comm));

    /// storage of nstates (states in nodes)
    DataHandle<RealVector> nstates = sockets.nstates.getDataHandle();
    const CFuint nbNodalVars = nstates[0].size();
    // initialNodalStates includes all nodal values for the whole mesh (NOT SCALABLE!!!)
    // we need to copy into the nstates only the nodal values corresponding to local nodes
    cf_assert(nstates.size()*nbNodalVars <= globalSize);
    DataHandle<Node*, GLOBAL> nodes = socket_nodes.getDataHandle();
    for (CFuint i = 0; i < nstates.size(); ++i) {
      const CFuint globalNodeID = nodes[i]->getGlobalID();
      const CFuint startn = globalNodeID*nbNodalVars;
      for (CFuint j = 0; j < nbNodalVars; ++j) {
        const CFuint idx = startn+j;
        cf_assert(idx < globalSize);
        nstates[i][j] = buffer[idx];
     }
    }
  }

  if (MeshDataStack::getActive()->getDataStorage()->checkData(dhName)) {
   MeshDataStack::getActive()->getDataStorage()->deleteData<CFreal>(dhName);
  }

  //initialize ParticleTracking
  m_lagrangianSolver.setDataSockets(sockets);
  m_lagrangianSolver.setupSendBufferSize(m_sendBufferSize);

  //initialize PostProcessign
  m_postProcess->setDataSockets(sockets);

  m_radiation->setupDataSockets(sockets);
  m_radiation->setup();
  // setting up the cell builder
  m_cellBuilder.setup();
  m_cellBuilder.getGeoBuilder()->setDataSockets(socket_states, socket_gstates, socket_nodes);


  CellTrsGeoBuilder::GeoData& cellData = m_cellBuilder.getDataGE();

//  How to get cell volumes:
//  cellData.idx = cellIdx;
//  GeometricEntity *const cell = m_cellBuilder.buildGE();
//  cell->computeVolume();
//  m_cellBuilder.releaseGE();

  SafePtr<TopologicalRegionSet> cells = MeshDataStack::getActive()->getTrs("InnerCells");
  cellData.trs = cells;

  // setting up the wall face builder
  m_wallFaceBuilder.setup();
  m_wallFaceBuilder.getGeoBuilder()->setDataSockets(socket_states, socket_gstates, socket_nodes);
  m_wallFaceBuilder.getDataGE().isBFace = true;

  const CFuint nCells = cellData.trs->getLocalNbGeoEnts();

  m_radiation->configureTRS();
  m_radiation->setupAxiFlag(m_isAxi);


  //Setup the mechanism counts to allocate the sync buffer correctly
  m_baseRadiatorPtr = m_radiation->getCellDistPtr(0)->getRadiatorPtr();
  m_HSNBRadiator=SharedPtr<HSNBRadiator>(static_cast<HSNBRadiator*>(&(*m_baseRadiatorPtr)));
  m_HSNBRadiator->setAbsorptionTolerance(m_tolerance);
  m_nbContinua=m_HSNBRadiator->getNbContinua();
  m_nbAtomics=m_HSNBRadiator->getNbAtoms();
  m_nbNonThickDiatomics=m_HSNBRadiator->getNbNonThickDiatomics();
  m_nbThickDiatomics=m_HSNBRadiator->getNbThickDiatomics();
  m_nbDiatomics=m_nbNonThickDiatomics+m_nbThickDiatomics;



  CFuint nbSpecies=PhysicalModelStack::getActive()->getImplementor()->
        getPhysicalPropertyLibrary<PhysicalChemicalLibrary>()->getNbSpecies();

  std::cout << "RadiativeTransferMonteCarloHSNB::setup => " << " m_nbContinua=" << m_nbContinua
            << " m_nbAtomics=" << m_nbAtomics << " m_nbNonThickDiatomics=" << m_nbNonThickDiatomics <<
            " m_nbThickDiatomics=" << m_nbThickDiatomics << " m_nbThickDiatomics=" << m_nbThickDiatomics << std::endl;

               m_paramSynchronizer.setup(m_myProcessRank,m_nbProcesses,m_nbThickDiatomics,m_nbNonThickDiatomics, m_nbContinua,m_nbAtomics,nbSpecies, m_HSNBRadiator->co2Exists());

  //std::cout << "NB DIATOMICS " << m_nbDiatomics<< " OF WHICH " << m_nbNonThickDiatomics << " ARE NONTHICK, NB ATOMICS " << m_nbAtomics << ", NB CONTINUA " << m_nbContinua  << std::endl;

  //Use the HSNBPhotonData to send the information needed to compose the variable length
  //MPI Datatype

  HSNBPhotonData hsnbPhotonData;

  //In doubt about this
  //CFint mechType= (int)hsnbPhotonData.mechType;

  int counts2[14] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1};
  MPIStructDef::buildMPIStruct<CFint,CFint,CFint,CFuint,CFuint,CFint,
           CFuint,CFuint,CFreal,CFreal,CFreal,CFreal,CFreal,CFreal>
          (&hsnbPhotonData.globalTrajectoryId, &hsnbPhotonData.indexWithinTrajectory,
           &hsnbPhotonData.fatherProcessId, &hsnbPhotonData.targetProcessId,
           &hsnbPhotonData.prevProcessId, &hsnbPhotonData.mechType,
           &hsnbPhotonData.mechanismIndex, &hsnbPhotonData.nbCrossedCells,
           &hsnbPhotonData.energyFraction, &hsnbPhotonData.wavelength,
           &hsnbPhotonData.KS, &hsnbPhotonData.transm, &hsnbPhotonData.energyResiduum,
           &hsnbPhotonData.wallTransmissivity,
           counts2 , HSNBUserdatatype);



  m_lagrangianSolver.setupParticleDatatype(HSNBUserdatatype.type);



  // resize the storage of qrad in each cell
  socket_qrad.getDataHandle().resize(nCells);
  m_stateInRadPowers.resize(nCells,0.);
  m_ghostStateInRadPowers.resize(socket_gstates.getDataHandle().size() , 0.);

  vector<string> wallTrsNames, boundaryTrsNames;
  m_radiation->getBoundaryTRSnames(boundaryTrsNames);
  m_radiation->getWallTRSnames(wallTrsNames);

  m_lagrangianSolver.setFaceTypes(wallTrsNames, boundaryTrsNames );

  // preallocation of memory for qradFluxWall
  CFuint nbFaces = 0;
  FaceTrsGeoBuilder::GeoData& WallFacesData = m_wallFaceBuilder.getDataGE();
  for(CFuint j=0; j<wallTrsNames.size(); ++j){
    SafePtr<TopologicalRegionSet> WallFaces = MeshDataStack::getActive()->getTrs(wallTrsNames[j]);
    WallFacesData.trs = WallFaces;
    nbFaces += WallFaces->getLocalNbGeoEnts();
  }
  socket_qradFluxWall.getDataHandle().resize(nbFaces);
}

/////////////////////////////////////////////////////////////////////////////

template<class PARTICLE_TRACKING>
void RadiativeTransferMonteCarloHSNB<PARTICLE_TRACKING>::executeOnTrs()
{
  using namespace std;
  using namespace COOLFluiD::Framework;
  using namespace COOLFluiD::MathTools;
  using namespace COOLFluiD::Common;
  using namespace COOLFluiD::LagrangianSolver;

  CFLog(VERBOSE, "RadiativeTransferMonteCarloHSNB::executeOnTrs() START\n");

  // for testcases == 0, the following must be recomputed everytime
  Stopwatch<WallTime> s;
  s.restart();

  // re-initialize the heat fluxes to 0
  //m_stateInRadPowers=0.;
  //m_RHSGas.assign(m_RHSGas.size(), 0.);



  MonteCarlo();
  computeHeatFlux();



  CFLog(INFO, "MonteCarlo() took " << s << "s\n");

  CFLog(VERBOSE, "RadiativeTransferMonteCarloHSNB::executeOnTrs() END\n");
}

/////////////////////////////////////////////////////////////////////////////

template<class PARTICLE_TRACKING>
void RadiativeTransferMonteCarloHSNB<PARTICLE_TRACKING>::getTotalEnergy()
{



  using namespace std;
  using namespace COOLFluiD::Framework;
  using namespace COOLFluiD::MathTools;
  using namespace COOLFluiD::Common;
  using namespace COOLFluiD::Numerics::FiniteVolume;
  using namespace COOLFluiD::LagrangianSolver;

  CFuint nbStates = m_radiation->getNbStates();
  CFuint nbGhostStates = m_radiation->getNbGhostStates();
  CFreal stateRadPower =  0.0;
  CFreal gStateRadPower =  0.0;

  m_stateRadPower.resize(nbStates);
  m_stateNetRadPower.resize(nbStates);
  m_nbPhotonsState.resize(nbStates);

  m_ghostStateRadPower.resize(nbGhostStates);
  m_nbPhotonsGhostState.resize(nbGhostStates);

  m_HSNBRadiator=SharedPtr<HSNBRadiator>(static_cast<HSNBRadiator*>(&(*m_baseRadiatorPtr)));


  //CFreal nbPhotons = nbStates ;//* m_nbRaysElem / m_radiation->getNumberLoops();
  //cout<<"nbPhotons: "<< nbPhotons <<endl;

  CFLog(VERBOSE, "RadiativeTransferMonteCarloHSNB::getTotalEnergy() => m_nbRaysTotal=" << m_nbRaysTotal << "\n");
  CFLog(VERBOSE, "RadiativeTransferMonteCarloHSNB::getTotalEnergy() => m_nbRaysCycle=" << m_nbRaysCycle << "\n");
  CFLog(VERBOSE, "RadiativeTransferMonteCarloHSNB::getTotalEnergy() => m_nbRaysCell=" << m_nbRaysElem << "\n");

  CFLog(INFO, "RadiativeTransferMonteCarloHSNB::getTotalEnergy() => Starting computation of emission for " << nbStates << " states. \n");
  boost::progress_display* progressBar = NULL;
  if (m_myProcessRank == 0) progressBar = new boost::progress_display(nbStates+nbGhostStates);

  if (m_dynamicRayDistribution) {

      //Get the total radiative power of the Cells

      if (m_HSNBRadiator->readEmissionData(m_stateRadPower, m_ghostStateRadPower,stateRadPower, gStateRadPower)==false) {
          //If the emission data had not been precomputed and loaded -> do the work now:
          for(CFuint state=0; state<nbStates; ++state){
              if ( !m_radiation->isStateNull(state) ){
                  //cout<<"is not ghost!"<<endl;

                  m_stateRadPower[state]= m_radiation->getCellDistPtr(state)
                          ->getRadiatorPtr()->getSpectraLoopPower();

                  m_stateNetRadPower=0.0;

                  stateRadPower += m_stateRadPower[state];


                  //cf_assert(isnormal(m_stateRadPower[state])==true);

                  CFLog(DEBUG_MED, "RadiativeTransferMonteCarlo::getTotalEnergy =>  m_stateRadPower[" << state << "]= " << m_stateRadPower[state] << ", currentCellVolume=" << m_radiation->getCellDistPtr(state)->getRadiatorPtr()->getCurrentCellVolume() << "\n");
                  //      std::cout << "RadiativeTransferMonteCarlo::getTotalEnergy => m_stateRadPower["<< state << "]= " << m_stateRadPower[state] << ", currentCellVolume=" << m_radiation->getCellDistPtr(state)->getRadiatorPtr()->getCurrentCellVolume() << "\n";
                  //m_baseRadiatorPtr = m_radiation->getCellDistPtr(state)->getRadiatorPtr();
                  //m_HSNBRadiator=SharedPtr<HSNBRadiator>(dynamic_cast<HSNBRadiator*>(&(*m_baseRadiatorPtr)));

              }
              else{
                  //cout<<"is ghost!"<<endl;
                  m_stateRadPower[state] = .0;
              }

              if (m_myProcessRank == 0)  ++*(progressBar);
          }


          //  for(CFuint state=0; state<nbStates; ++state){
          //    m_HSNBRadiator=SharedPtr<HSNBRadiator>(dynamic_cast<HSNBRadiator*>(&(*m_baseRadiatorPtr)));
          //    m_stateRadPower[state]/=m_HSNBRadiator->getTotalEmissivePower();
          //  }

          //Get the total radiative power of the Wall Faces
          for(CFuint gstate=0; gstate<nbGhostStates; ++gstate){
              if ( !m_radiation->isGhostStateNull(gstate) ){
                  m_ghostStateRadPower[gstate]= m_radiation->getWallDistPtr(gstate)
                          ->getRadiatorPtr()->getSpectraLoopPower();

                  //      cout<<"m_ghostStateRadPower[" << gstate << "]=" << m_ghostStateRadPower[gstate]<<endl;

                  gStateRadPower += m_ghostStateRadPower[gstate];
              }
              else{
                  //cout<<"is ghost!"<<endl;
                  m_ghostStateRadPower[gstate] = 0.;
              }

              if (m_myProcessRank == 0)  ++*(progressBar);
          }
          //cout<<" walls rad power : "<< gStateRadPower<<endl;
      }


      m_totalSubdomainRadPower = stateRadPower + gStateRadPower;



      //std::cout << "Total emitted power = " << m_totalSubdomainRadPower << std::endl;

      MPI_Allreduce(&m_totalSubdomainRadPower, &m_totalGlobalRadPower, 1,
                    Common::MPIStructDef::getMPIType(&m_totalSubdomainRadPower), MPI_SUM, m_comm);

      //Determine raynumber to be distributed among this domain according to the domain's share of the total emissive power
      m_nbRaysTotal*=m_totalSubdomainRadPower/m_totalGlobalRadPower;


      CFLog(INFO, "RadiativeTransferMonteCarlo::getTotalEnergy =>  total rad. power is " << m_totalGlobalRadPower << "\n");

      std::cout << "RadiativeTransferMonteCarloHSNB::getTotalEnergy => P " << m_myProcessRank << " emits " << m_totalSubdomainRadPower << " / " << m_totalGlobalRadPower << " which corresponds to " <<
                   100*m_totalSubdomainRadPower/m_totalGlobalRadPower << "% or "<< m_nbRaysTotal << " rays for "<< nbStates+nbGhostStates
                << " cells / faces. The ratio rays/cells is "<< m_nbRaysTotal / (nbStates+nbGhostStates) << ". \n";

      CFuint curPhotonCount=0;

//      std::cout<<"statePhotons"<<endl;
      for(CFuint i= 0; i< m_stateRadPower.size(); ++i){

          m_nbPhotonsState[i] =( m_stateRadPower[i] > 0. )? (m_nbRaysTotal*m_stateRadPower[i]/m_totalSubdomainRadPower) : 0 ;
          curPhotonCount+=m_nbPhotonsState[i];

          // cout<<"nbPhotons: m_nbPhotonsState[i]="<<m_nbPhotonsState[i]<< ' '<<m_stateRadPower[i] <<endl;
      }

//      std::cout<<"ghostPhotons"<<endl;
      for(CFuint i= 0; i< m_ghostStateRadPower.size(); ++i){

          m_nbPhotonsGhostState[i] =( m_ghostStateRadPower[i] > 0.)? (m_nbRaysTotal*m_ghostStateRadPower[i]/m_totalSubdomainRadPower) : 0 ;
          curPhotonCount+=m_nbPhotonsGhostState[i];

          //cout<<"nbPhotons: m_nbPhotonsGhostState[i]="<<m_nbPhotonsGhostState[i]<< ' '<<m_ghostStateRadPower[i] <<endl;
      }


  }
  else {

      //Associate a fixed number of rays with every cell in the grid which is potentially wasteful but might improve resolution at boundaries:
      //Get the total radiative power of the Cells
      m_nbRaysTotal=0;

      //Try to read a precomputed table first (if it exists)
      if (m_HSNBRadiator->readEmissionData(m_stateRadPower, m_ghostStateRadPower,stateRadPower, gStateRadPower)==false) {
          //If the emission data could not be read from a precomputed table -> do the work now!
          for(CFuint state=0; state<nbStates; ++state){
              if ( !m_radiation->isStateNull(state) ){
                  //cout<<"is not ghost!"<<endl;
                  m_stateRadPower[state]= m_radiation->getCellDistPtr(state)
                          ->getRadiatorPtr()->getSpectraLoopPower();
                  stateRadPower += m_stateRadPower[state];
                  //cout<<"state radPower: "<<m_stateRadPower[state]<<endl;
                  //cout<<"THE OTER: m_axi Volume: "<<m_axiVolumes[i]<<endl;
                  //m_stateRadPower[i] = cellK*sigma*pow(T,4.)*m_axiVolumes[i];
                  CFLog(DEBUG_MED, "RadiativeTransferMonteCarlo::getTotalEnergy =>  m_stateRadPower[" << state << "]= " << m_stateRadPower[state] << ", currentCellVolume=" << m_radiation->getCellDistPtr(state)->getRadiatorPtr()->getCurrentCellVolume() << "\n");
              }
              else{
                  //cout<<"is ghost!"<<endl;
                  m_stateRadPower[state] = .0;
              }

              if (m_myProcessRank == 0)  ++*(progressBar);
          }

          //Get the total radiative power of the Wall Faces
          for(CFuint gstate=0; gstate<nbGhostStates; ++gstate){
              if ( !m_radiation->isGhostStateNull(gstate) ){
                  m_ghostStateRadPower[gstate]= m_radiation->getWallDistPtr(gstate)
                          ->getRadiatorPtr()->getSpectraLoopPower();
                  //cout<<"gstate radPower: "<<m_ghostStateRadPower[gstate]<<endl;

                  gStateRadPower += m_ghostStateRadPower[gstate];
              }
              else{
                  //cout<<"is ghost!"<<endl;
                  m_ghostStateRadPower[gstate] = 0.;
              }

              if (m_myProcessRank == 0)  ++*(progressBar);
          }
          //cout<<" walls rad power : "<< gStateRadPower<<endl;
      }

      m_totalRadPower = stateRadPower + gStateRadPower;
      std::cout << "P " << m_myProcessRank << " m_totalRadPower=" << m_totalRadPower << "\n";

      //cout<<"statePhotons"<<endl;
      CFuint nbRaysState=0;
      for(CFuint i= 0; i< m_stateRadPower.size(); ++i){
          //CFuint test = CFuint(m_stateRadPower[i]/m_totalRadPower*nbPhotons);
          nbRaysState=( m_stateRadPower[i] > 0. )? m_nbRaysElem : 0 ;
          m_nbPhotonsState[i] =nbRaysState;
          m_nbRaysTotal+=nbRaysState;
          //cout<<"nbPhotons: "<<m_nbPhotonsState[i]<< ' '<<m_stateRadPower[i] <<endl;
      }

      //cout<<"ghostPhotons"<<endl;
      for(CFuint i= 0; i< m_ghostStateRadPower.size(); ++i) {
          //CFuint test = CFuint(m_stateRadPower[i]/m_totalRadPower*nbPhotons);
          nbRaysState=( m_ghostStateRadPower[i] > 0.)? m_nbRaysElem : 0 ;
          m_nbPhotonsGhostState[i] = nbRaysState;
          m_nbRaysTotal+=nbRaysState;
          //cout<<"nbPhotons: "<<m_nbPhotonsGhostState[i]<< ' '<<m_ghostStateRadPower[i] <<endl;
      }

  }


//  MPI_Reduce(&m_nbRaysTotal, &m_maxRaysSingleProcess, 1,
//        Common::MPIStructDef::getMPIType(&m_nbRaysTotal), MPI_MAX, 0, m_comm);


  //If no precomputed table exists yet -> write it to a file so that it can be reused
  m_HSNBRadiator->exportEmissionData(m_stateRadPower, m_ghostStateRadPower,stateRadPower, gStateRadPower);


delete progressBar;

//  std::cout << "P "  << m_myProcessRank << " finished getTotalEnergy() \n";
}

/////////////////////////////////////////////////////////////////////////////

template<class PARTICLE_TRACKING>
bool RadiativeTransferMonteCarloHSNB<PARTICLE_TRACKING>::
getCellPhotonData(Photon &ray)
{
  using namespace std;
  using namespace COOLFluiD::Framework;
  using namespace COOLFluiD::MathTools;
  using namespace COOLFluiD::Common;
  using namespace COOLFluiD::Numerics::FiniteVolume;
  using namespace COOLFluiD::LagrangianSolver;

  //static CellTrsGeoBuilder::GeoData& cellData = m_cellBuilder.getDataGE();
  const CFuint nbStates = m_nbPhotonsState.size();
  CFuint nbRaysCell;

//  std::cout << "Enter loop with m_istate_cell_fix=" << m_istate_cell_fix << "\n";
//  CFLog(VERBOSE, "ENTER THE LOOP " << m_istate_cell_fix << "\n");

  //m_energyCheck=0.0;
  for(;m_istate_cell_fix<nbStates; ++m_istate_cell_fix) {
    CFLog(DEBUG_MAX, "RadiativeTransferMonteCarloHSNB::getCellPhotonData => Current state is "  << m_istate_cell_fix << "\n");
    for(;m_iphoton_cell_fix<m_nbPhotonsState[ m_istate_cell_fix ];) {
//              cout<<"P " << m_myProcessRank<< ": m_nbPhotonsState[state]: "<<m_nbPhotonsState[m_istate_cell_fix]<<" m_iphoton_cell_fix="
//                 << m_iphoton_cell_fix << "\n";

      CFLog(DEBUG_MIN, "RadiativeTransferMonteCarloHSNB::getCellPhotonData => m_istate_cell_fix=" << m_istate_cell_fix << ", m_iphoton_cell_fix=" << m_iphoton_cell_fix << ", m_nbPhotonsState[=" << m_istate_cell_fix << "]="
            << m_nbPhotonsState[ m_istate_cell_fix ]  << "\n");

      //cellData.idx =  state;


      //BT: NOTE: I rewrote getRandomEmission in the HSNBRadiator to only generate the direction
      // wavelength has to be set in HSNBRadiator::generatePhotonData
      m_radiation->getCellDistPtr( m_istate_cell_fix )->
    getRadiatorPtr()->getRandomEmission(ray.userData.wavelength, m_direction );

      //cout<<"ray directions: ";
//      CFLog(INFO, "RadiativeTransferMonteCarloHSNB::getCellPhotonData => Emission direction: [ " );
      for(CFuint ii=0; ii<m_dim2; ++ii){
          ray.commonData.direction[ii]=m_direction[ii];
//          CFLog(INFO, ray.commonData.direction[ii] << " ");
      }
//      CFLog(INFO, "] \n" );



      //Get the beam max optical path Ks
      ray.userData.KS = - std::log( m_rand.uniformRand() );

      //Get cell center
      static Framework::DataHandle<Framework::State*, Framework::GLOBAL> states
    = socket_states.getDataHandle();

      //NOTE: This might be replaced by a computation for a random location in the cell volume
      //which might speed up convergence as demonstrated in tests with the 1D mesh
      Node& baricenter = (*states[ m_istate_cell_fix ]).getCoordinates();


      CFLog(DEBUG_MED, "RadiativeTransferMonteCarloHSNB::getCellPhotonData => Emission position: [ " );
      for(CFuint i=0;i<m_dim;++i){
          ray.commonData.currentPoint[i]=baricenter[i];
          CFLog(DEBUG_MED, ray.commonData.currentPoint[i] << " ");
      }
      CFLog(DEBUG_MED, "] \n" );

      for(CFuint i=m_dim;i<m_dim2;++i){
          //  cout<<baricenter[i]<<' ';
          ray.commonData.currentPoint[i] = 0.;
      }



      //CFLog(INFO, "RadiativeTransferMonteCarloHSNB::getCellPhotonData BEFORE INIT=>m_istate_nextCell_fix=" << m_istate_nextCell_fix
      //      << ", ray.userData.energyFraction="
      //      << ray.userData.energyFraction << ", ray.userData.mechType=" << ray.userData.mechType
      //      << ", ray.userData.wavelength=" << ray.userData.wavelength << "\n");

      //Assign an energy, emitting mechanism and a wavelength / wavenumber to the new photon using the HSNBRadiator
      m_HSNBRadiator->generateCellPhotonData(m_istate_nextCell_fix, ray.userData.energyFraction,ray.userData.mechanismIndex, ray.userData.mechType,ray.userData.wavelength);

//       DEBUG ONLY
//            CFLog(INFO, "RadiativeTransferMonteCarloHSNB::getCellPhotonData => Emission direction: [ " );
//            for(CFuint ii=0; ii<m_dim2; ++ii){
//                ray.commonData.direction[ii]=ii;
//                CFLog(INFO, ray.commonData.direction[ii] << " ");
//            }
//            CFLog(INFO, "] \n" );
//            ray.userData.energyFraction=61.4424;
//            ray.userData.mechanismIndex=1;
//            ray.userData.mechType=0;
//            ray.userData.wavelength=4.0448e-05;
//            ray.commonData.direction[0]=  0.128173;
//            ray.commonData.direction[1]= -0.297606;
//            ray.commonData.direction[2]=  0.946046;



//      CFLog(INFO, "RadiativeTransferMonteCarloHSNB::getCellPhotonData AFTER INIT=>m_istate_nextCell_fix=" << m_istate_nextCell_fix
//            << ", ray.userData.energyFraction="
//            << ray.userData.energyFraction << ", ray.userData.mechType=" << ray.userData.mechType
//            << ", ray.userData.wavelength=" << ray.userData.wavelength << ", ray.userData.mechanismIndex=" << ray.userData.mechanismIndex<< "\n");

      ray.userData.energyResiduum=ray.userData.energyFraction;

      m_emittedPhotonEnergyState+=ray.userData.energyFraction;
      m_stateNetRadPower[m_istate_cell_fix]+=ray.userData.energyFraction;

      //m_energyCheck+=ray.userData.energyFraction;

      //Get the remaining information
      //const CFuint cellID = cell->getID();
      CFuint cellID = m_radiation->getCurrentCellStateID();
      ray.commonData.cellID = cellID;

      ray.userData.prevProcessId=m_myProcessRank;
      ray.userData.nbCrossedCells=0;
      ray.userData.transm=1.0;
      ray.userData.wallTransmissivity=1.0;
      ray.userData.fatherProcessId=m_myProcessRank;

      cf_assert(!isnan(ray.userData.energyFraction));

      CFLog(DEBUG_MAX, "RadiativeTransferMonteCarloHSNB::getCellPhotonData() => Generate photon. MechID: " << ray.userData.mechanismIndex << ", mechType: "<< ray.userData.mechType<<", lambda: " <<
            ray.userData.wavelength << "\n");

      //ray.userData.energyFraction= m_stateRadPower[ m_istate_cell_fix ]/CFreal(m_nbPhotonsState[ m_istate_cell_fix ]);

//      CFLog(INFO, "RadiativeTransferMonteCarloHSNB::getCellPhotonData => " << std::abs(m_emittedPhotonEnergyState-m_stateRadPower[m_istate_cell_fix]) << "\n");

//      cf_assert(std::abs(m_emittedPhotonEnergyState-m_stateRadPower[m_istate_cell_fix])<m_tol);
//      std::cout << "ASSERT ENERGY CONSERVATION, CURRENT PHOTON ENERGY="<< ray.userData.energyFraction <<  ", m_emittedPhotonEnergyState=" << m_emittedPhotonEnergyState << " / " << m_stateRadPower[m_istate_cell_fix] << "\n";

//      cout << ray.userData.indexWithinTrajectory << " INDEX OF NEW RAY" << endl;
      //m_cellBuilder.releaseGE();
      ++m_iphoton_cell_fix;
      //CFLog(DEBUG_MED, "EXIT THE INNER LOOP m_istate_cell_fix:" << m_istate_cell_fix << ", m_iphoton_cell_fix: "<< m_iphoton_cell_fix << "\n");
      return true;
    }
    //CFLog(INFO, "RadiativeTransferMonteCarloHSNB::getCellPhotonData => m_energyCheck=" << m_energyCheck << " / m_stateRadPower[m_istate_cell_fix]=" << m_stateRadPower[m_istate_cell_fix] << "\n");
    //m_energyCheck=0.0;
    m_istate_nextCell_fix=m_istate_cell_fix+1;
    m_emittedPhotonEnergyState=0.0;

    //If we are not done with all cells: Get ray distribution for the next state
    if (m_istate_nextCell_fix!=nbStates) {
        CFLog(DEBUG_MED, "RadiativeTransferMonteCarloHSNB::getCellPhotonData => Set state "<< m_istate_nextCell_fix << "\n");
        //BT: Getting the CellDistPtr(state) sets the current state in the radiator
        m_baseRadiatorPtr = m_radiation->getCellDistPtr(m_istate_nextCell_fix)->getRadiatorPtr();
        m_HSNBRadiator=SharedPtr<HSNBRadiator>(static_cast<HSNBRadiator*>(&(*m_baseRadiatorPtr)));

        nbRaysCell=m_nbPhotonsState[m_istate_nextCell_fix];
        CFLog(DEBUG_MAX, "RadiativeTransferMonteCarloHSNB::getCellPhotonData => m_nbPhotonsState[" << m_istate_nextCell_fix
              << "]=" << nbRaysCell << "\n");

        //Setup the next state: The total number of rays for this cell is distributed among all emitting mechanisms
        //according to the fraction of the total emission associated with each mechanism
        //HAS TO BE DONE BEFORE GENERATING PHOTONS IN A NEW CELL
        m_HSNBRadiator->setState(nbRaysCell);

        //CFLog(INFO, "RadiativeTransferMonteCarloHSNB::getCellPhotonData => EXPECT=" << m_stateRadPower[m_istate_nextCell_fix] << " / ACTUAL="  << m_HSNBRadiator->getCurrentCellEmissivePower() << "\n");
        //CFLog(INFO, "RadiativeTransferMonteCarloHSNB::getCellPhotonData => m_nbPhotonsState[" << m_istate_nextCell_fix << "]="<< nbRaysCell << "\n");
    }


    m_iphoton_cell_fix=0;
  }
  return false;
}

/////////////////////////////////////////////////////////////////////////////

template<class PARTICLE_TRACKING>
bool RadiativeTransferMonteCarloHSNB<PARTICLE_TRACKING>::getFacePhotonData(Photon &ray)
{
  using namespace std;
  using namespace COOLFluiD::Framework;
  using namespace COOLFluiD::MathTools;
  using namespace COOLFluiD::Common;
  using namespace COOLFluiD::Numerics::FiniteVolume;
  using namespace COOLFluiD::LagrangianSolver;

  //static CellTrsGeoBuilder::GeoData& cellData = m_cellBuilder.getDataGE();
  //m_ photon=0, gState=0;
  const CFuint nbGstates = m_nbPhotonsGhostState.size();

  for( ; m_igState_face_fix < nbGstates ; ++ m_igState_face_fix)
  {
    for( ; m_iphoton_face_fix < m_nbPhotonsGhostState[ m_igState_face_fix  ]; )
    {
      //Get directions
      m_radiation->getWallDistPtr( m_igState_face_fix )->
    getRadiatorPtr()->getRandomEmission(ray.userData.wavelength, m_direction );



      const CFuint faceGeoID = m_radiation->getCurrentWallGeoID();
      const CFuint cellID = m_lagrangianSolver.getWallStateId( faceGeoID );

      static Framework::DataHandle<Framework::State*, Framework::GLOBAL> states
    = socket_states.getDataHandle();

      Node& cellCenter = (*states[cellID]).getCoordinates();

      //cout<<"direction = [";
      for(CFuint ii=0; ii<m_dim; ++ii){
    //    cout<< -directions[ii] << ' ';
        ray.commonData.direction[ii]= m_direction[ii];
      }
      //cout<<" ]; "<<endl;

      //Get the beam max optical path Ks
      ray.userData.KS = - std::log( m_rand.uniformRand() );
      // ray.actualKS = 0;

      ray.userData.mechType=WALLEMISSION;
      ray.userData.fatherProcessId=m_myProcessRank;
      ray.userData.prevProcessId=m_myProcessRank;
      ray.userData.nbCrossedCells=0;


      //Get the face center
      DataHandle<CFreal> faceCenters = socket_faceCenters.getDataHandle();

      //cout<<"baricenter= [ ";

      //cout<<" ]; "<<endl;

      m_cartPosition3[0] = ray.commonData.currentPoint[0];
      m_cartPosition3[1] = ray.commonData.currentPoint[1];
      m_cartPosition3[2] = ray.commonData.currentPoint[2];

      m_lagrangianSolver.getNormals(faceGeoID, m_cartPosition3, m_faceNormal3);

      // small correction to make sure the initial point is inside the cell
      // move the initial point 1% closer to the cell center
      CFreal tCenter = 0.;
      for(CFuint i=0;i<m_dim;++i){
          tCenter += m_faceNormal3[i] * (faceCenters[i] - cellCenter[i]);
      }

      for(CFuint i=0;i<m_dim;++i){
          //cout<<faceCenters[m_dim * faceGeoID + i]<<' ';
          ray.commonData.currentPoint[i]=faceCenters[m_dim * faceGeoID + i] +
                  .01 * m_faceNormal3[i]*tCenter;
      }

      //cout<<"normal= [" << m_faceNormal3 <<" ]; "<<endl;

      if (m_isAxi) {
          //rotate the position and vector to a random theta
          const CFreal theta = m_rand.uniformRand(-3.141516, 3.141516);
          const CFreal x = ray.commonData.currentPoint[0];
          const CFreal y = ray.commonData.currentPoint[1];

          ray.commonData.currentPoint[0] = x;
          ray.commonData.currentPoint[1] = y*std::cos(theta);
          ray.commonData.currentPoint[2] = y*std::sin(theta);

          m_cartPosition3[0] = ray.commonData.currentPoint[0];
          m_cartPosition3[1] = ray.commonData.currentPoint[1];
          m_cartPosition3[2] = ray.commonData.currentPoint[2];

          m_lagrangianSolver.getNormals(faceGeoID, m_cartPosition3, m_faceNormal3);

          m_rand.hemiDirections(3, m_faceNormal3, m_sOut3);

          ray.commonData.direction[0] = m_sOut3[0];
          ray.commonData.direction[1] = m_sOut3[1];
          ray.commonData.direction[2] = m_sOut3[2];
      }

      //Get the remaining information
      ray.commonData.cellID = cellID;

      ray.userData.energyFraction= m_ghostStateRadPower[m_igState_face_fix]/
    CFreal(m_nbPhotonsGhostState[m_igState_face_fix]);

      ray.userData.energyResiduum=ray.userData.energyFraction;
      ray.userData.transm=1.0;
      ray.userData.nbCrossedCells=0;
      ray.userData.wallTransmissivity=1.0;
      ray.userData.prevProcessId=m_myProcessRank;
      ray.userData.fatherProcessId=m_myProcessRank;

      if (isnormal(ray.userData.energyFraction)==false) {
          std::cout << "RadiativeTransferMonteCarloHSNB::getFacePhotonData => Singularity: ray.userData.energyFraction=" << ray.userData.energyFraction << "\n";
      }



      ++m_iphoton_face_fix;
      return true;
    }
    m_iphoton_face_fix=0;
  }
  return false;
}
/////////////////////////////////////////////////////////////////////////////

template<class PARTICLE_TRACKING>
void RadiativeTransferMonteCarloHSNB<PARTICLE_TRACKING>::checkSyncForceConditions()
{
    if (m_enforceSync==false) {
//        if (m_nbCrossedCellsCycle>=m_maxVisitedCellsCycle) {
//            m_enforceSync=true;
//        }
        if (m_paramSynchronizer.totalSendBufferByteSize()>=m_maxTraceBufferByteSize) {
//            std::cout << "RadiativeTransferMonteCarloHSNB::checkSyncForceConditions => m_paramSynchronizer.totalSendBufferByteSize=" << m_paramSynchronizer.totalSendBufferByteSize() << "\n";
            m_enforceSync=true;
        }
        else if (m_nbTracesCommittedCycle>m_maxTraceBufferTracesCommitted){
            m_enforceSync=true;
        }
        else if (m_stpTracing.read()>m_maxSecondsBetweenSyncs) {
            m_enforceSync=true;
        }

    }
}



/////////////////////////////////////////////////////////////////////////////

template<class PARTICLE_TRACKING>
void RadiativeTransferMonteCarloHSNB<PARTICLE_TRACKING>::computePhotons()
{
  using namespace std;
  using namespace COOLFluiD::Framework;
  using namespace COOLFluiD::MathTools;
  using namespace COOLFluiD::Common;
  using namespace COOLFluiD::Numerics::FiniteVolume;
  using namespace COOLFluiD::LagrangianSolver;

  CFLog(DEBUG_MAX, "RadiativeTransferMonteCarloHSNB::computeCellRays()\n");

  m_rand.seed(time(NULL)*(m_myProcessRank+1));

  CFLog(INFO, "RadiativeTransferMonteCarloHSNB::computeCellRays => seed=" << time(NULL)*(m_myProcessRank+1) << "\n");

  // CFuint totalnbPhotons =  (m_nbRaysElem )* m_radiation->getNbStates();

  CFuint toGenerateCellPhotons = 0;
  for(CFuint i=0;i<m_nbPhotonsState.size();++i){
    toGenerateCellPhotons+=m_nbPhotonsState[i];
  }



  CFuint toGenerateWallPhotons = 0;
  for(CFuint i=0;i<m_nbPhotonsGhostState.size();++i){
    toGenerateWallPhotons+=m_nbPhotonsGhostState[i];
  }


  //cout<<nbPhotons<<' '<<totalnbPhotons<<endl;

  //  cout<<"number of photons to emmit: "<<m_nbRaysCycle<< " "<<m_sendBufferSize<<endl;
  CFuint recvSize = 0;
  bool done = false;

  Photon photon;
  m_debugPhotonsComitted=0;

  //std::cout << " COMPUTE PHOTONS \n";

  m_avgCellsPerPhotonCycle=0;


  CFuint totalnbPhotons = toGenerateWallPhotons + toGenerateCellPhotons;

  CFLog(INFO, "RadiativeTransferMonteCarloHSNB::getCellPhotonData => totalnbPhotons=" << totalnbPhotons << ", toGenerateWallPhotons=" << toGenerateWallPhotons << ", toGenerateCellPhotons=" << toGenerateCellPhotons << "\n");

  //cf_assert(totalnbPhotons==m_nbRaysTotal);


  CFuint totalNbPhotonsLeftGlobal=totalnbPhotons;
  if (m_trackTelemetry) {
      //Get actual total number of photons:
      MPI_Reduce(&totalnbPhotons, &totalNbPhotonsLeftGlobal, 1,
                                                   Common::MPIStructDef::getMPIType(&totalnbPhotons), MPI_SUM, 0, m_comm);
      CFLog(INFO, "RadiativeTransferMonteCarloHSNB::getCellPhotonData => The net total photon count is " << totalNbPhotonsLeftGlobal << "\n");
  }


  boost::progress_display* progressBar = NULL;
  if (m_myProcessRank == 0) progressBar = new boost::progress_display(totalnbPhotons);

  //CFuint toGeneratePhotons = totalnbPhotons;



  Stopwatch<WallTime> stpIdle;
  Stopwatch<WallTime> s;
  s.restart();

  m_nbIterations=0;

  CFuint nbCellsBeforeRaytracing=0;
  CFuint totalNbPhotonsLeftGlobalLastCycle=totalNbPhotonsLeftGlobal;

  m_avgCellsPerPhotonTotal=0;
  m_avgPhotonsCycle=0;
  m_avgIdleTime=0;
  m_avgTimePerPhoton=0;
  m_avgCycleTime=0.0;
  m_avgCellsPerPhotonGlobal=0;
  m_avgTimePerPhotonGlobal=0;
  m_avgProcessEfficiency=0.0;
  m_avgNewPhotonsTracedCycle=0;
  m_debugTraceUID=0;
  m_totalWallAbsorption=0.0;
  m_nbCrossingForty=0;
  m_emittedPhotonEnergyState=0.0;
  CFuint genPhotonsCycle=0;

  CFuint maxNbPhotonsLeftGlobal=0;
  CFuint totalNbPhotonsLeftLocal=0;
  CFuint approxSecondsLeft=0.0;
  vector< Photon > photonStack;
  photonStack.reserve(m_sendBufferSize);

  //SET FIRST STATE
  m_baseRadiatorPtr = m_radiation->getCellDistPtr(m_istate_cell_fix)->getRadiatorPtr();
  m_HSNBRadiator=SharedPtr<HSNBRadiator>(static_cast<HSNBRadiator*>(&(*m_baseRadiatorPtr)));

  //Setup the first state: The total number of rays for this cell is distributed among all emitting mechanisms
  //according to the fraction of the total emission associated with each mechanism
  //HAS TO BE DONE BEFORE GENERATING PHOTONS IN A NEW CELL
  m_HSNBRadiator->setState(m_nbPhotonsState[ m_istate_cell_fix ]);

  bool stopPhotonGenOnBufferExceed  = (m_maxGlobalBufferByteSize != 0);

  if (stopPhotonGenOnBufferExceed) {
      CFLog(VERBOSE, "RadiativeTransferMonteCarloHSNB::getCellPhotonData => Maximum allowable total buffer usage is restricted to " << m_maxGlobalBufferByteSize << " bytes. \n");
  }

  while( !done ){

    genPhotonsCycle=0;
    m_nbTracesCommittedCycle=0;
    m_nbPhotonsTracedCycle=0;
    m_avgCellsPerPhotonCycle=0;
    m_nbIterations++;
    m_enforceSync=false;
    m_nbCrossedCellsCycle=0;
    m_stpCycle.restart();
    m_stpTracing.restart();

    recvSize = photonStack.size();
    //generate and raytrace the inner photons




    if (stopPhotonGenOnBufferExceed) {
        m_enforceSync=m_paramSynchronizer.totalMemoryUsageExceeds(m_maxGlobalBufferByteSize);
    }

    CFuint nbCellPhotons =
        std::min(std::max(CFint(m_nbRaysCycle) - CFint(recvSize),(CFint)0), CFint(toGenerateCellPhotons) );

    CFuint nbWallPhotons =
        std::min(std::max(CFint(m_nbRaysCycle) - CFint(recvSize) - CFint(nbCellPhotons),(CFint)0), CFint(toGenerateWallPhotons ));

    CFLog(DEBUG_MED, "RadiativeTransferMonteCarloHSNB::getCellPhotonData => EXPECT=" << m_stateRadPower[m_istate_cell_fix] << " / ACTUAL="  << m_HSNBRadiator->getCurrentCellEmissivePower() << "\n");
    //CFLog(INFO, "RadiativeTransferMonteCarloHSNB::computePhotons => m_nbPhotonsState[" << m_istate_cell_fix << "]="<< m_nbPhotonsState[ m_istate_cell_fix ] << "\n");

    for(CFuint i = 0; i< photonStack.size(); ++i ){

        //        printPhoton(photonStack[i]);
        m_curTraceSet=m_paramSynchronizer.getTraceStackAt(i);

        cf_assert(m_curTraceSet.trace->m_nbCellsCrossed==photonStack[i].userData.nbCrossedCells);

        // DEBUG RECEIVE PARAMETERS
        // m_curTraceSet.trace->print();
        // m_curTraceSet.print();

        nbCellsBeforeRaytracing=photonStack[i].userData.nbCrossedCells;
        rayTracing( photonStack[i] );
        //        rayTracingFullAbsorption(photonStack[i]);
        m_avgCellsPerPhotonCycle+=int(((photonStack[i].userData.nbCrossedCells-nbCellsBeforeRaytracing)-m_avgCellsPerPhotonCycle))/m_nbPhotonsTracedCycle;
        m_curTraceSet.reset();
    }

    checkSyncForceConditions();


    for(CFuint i=0; i < nbCellPhotons ; ++i ){
        //CFLog(INFO, "RadiativeTransferMonteCarloHSNB::computePhotons => i=" << i << ", nbCellPhotons=" << nbCellPhotons << ", toGenerateCellPhotons=" << toGenerateCellPhotons << "\n");

        //Stop tracing and skip to sync step.
        //The sync conditions are checked in the CheckSyncConditions().
        if (m_enforceSync) {
            break;
        }

        if(getCellPhotonData( photon )){
            //CFLog(INFO,"PHOTON: " << photon.userData.energyFraction<<' '<<photon.userData.KS<<'\n' );
//            std::cout << "P " << m_myProcessRank <<  " PHOTON: " << photon.userData.energyFraction<<' '<<photon.userData.KS<<'\n';
            //printPhoton(photon);

            //Use m_curTraceSet as second interface to photon Data, we can access the photon header data and
            //the HSNB trace from it
            m_curTraceSet.reset();
            m_curTraceSet.headerData=&photon.userData;
            genPhotonsCycle++;

            rayTracing(photon);

//          rayTracingFullAbsorption(photon);




            m_avgCellsPerPhotonCycle+=int((photon.userData.nbCrossedCells-m_avgCellsPerPhotonCycle))/m_nbPhotonsTracedCycle;

        }


        --toGenerateCellPhotons;
        if (m_myProcessRank == 0)  ++*(progressBar);
    }


    for(CFuint i=0; i < nbWallPhotons ; ++i ){

        //Stop tracing and skip to sync step.
        //The sync conditions are checked in the CheckSyncConditions().
        if (m_enforceSync) {
            break;
        }

        if(getFacePhotonData( photon )){
            //printPhoton(photon);
            m_curTraceSet.reset();
            m_curTraceSet.headerData=&photon.userData;
            genPhotonsCycle++;
            rayTracing(photon);
//          rayTracingFullAbsorption(photon);

            m_avgCellsPerPhotonCycle+=int((photon.userData.nbCrossedCells-m_avgCellsPerPhotonCycle))/m_nbPhotonsTracedCycle;
        }
        -- toGenerateWallPhotons;
        if (m_myProcessRank == 0)  ++*(progressBar);
    }

//    CFLog(INFO, "raytrace the outer photons \n");


    stpIdle.restart();
    m_stpTracing.stop();
    CFLog(DEBUG_MED, "RadiativeTransferMonteCarloHSNB::computePhotons => Raytracing took " << m_stpTracing.read() << " seconds. Start waiting for other processes to sync... \n");



    //m_paramSynchronizer.printSendTraces();



    //sincronize
    //      CFLog(INFO, "sincronizing\n");
    bool isLastPhoton = (toGenerateCellPhotons + toGenerateWallPhotons == 0);
    done = m_lagrangianSolver.sincronizeParticles(photonStack, isLastPhoton);
    stpIdle.stop();
//    cout << "P " << m_myProcessRank << " waited approx. " << stpIdle.read() << " seconds for other processes to call sync. NbCellsCrossed=" << m_nbCrossedCellsCycle <<  "\n";

//    Debug the synchronization of m_paramSynchronizer serially
//    int active=0;
//    if (m_myProcessRank == 0) {
//        active = 1;

//        m_paramSynchronizer.printSendTraces(true);

//        MPI_Send(&active, 1, MPIStructDef::getMPIType(&active), m_myProcessRank+1, 0, m_comm);
//    }
//    else {
//        MPI_Recv(&active, 1, MPIStructDef::getMPIType(&active), m_myProcessRank-1, 0, m_comm, MPI_STATUS_IGNORE);

//        m_paramSynchronizer.printSendTraces(true);

//        if (m_myProcessRank!=m_nbProcesses) {
//            MPI_Send(&active, 1, MPIStructDef::getMPIType(&active), m_myProcessRank+1, 0, m_comm);
//        }
//    }

//    CFLog(INFO, "RadiativeTransferMonteCarloHSNB::computePhotons =>  m_maxVisitedCellsCycle / avg =" << m_maxVisitedCellsCycle << "\n");

    //Send tons of HSNB Parameters
    if (done==false) {
        m_paramSynchronizer.synchronize(photonStack);
    }

    //CFLog(INFO, "PRINT RECV TRACES \n");

//    Debug the synchronization of m_paramSynchronizer serially
//    std::cout << "\n";
//    if (m_myProcessRank == 0) {
//        active = 1;

//        m_paramSynchronizer.printRecvTraces(true);

//        MPI_Send(&active, 1, MPIStructDef::getMPIType(&active), m_myProcessRank+1, 0, m_comm);
//    }
//    else {
//        MPI_Recv(&active, 1, MPIStructDef::getMPIType(&active), m_myProcessRank-1, 0, m_comm, MPI_STATUS_IGNORE);

//        m_paramSynchronizer.printRecvTraces(true);

//        if (m_myProcessRank!=m_nbProcesses) {
//            MPI_Send(&active, 1, MPIStructDef::getMPIType(&active), m_myProcessRank+1, 0, m_comm);
//        }
//    }

    m_stpCycle.stop();

    //Mesure performance
    ///////////////////////////////////

    if (m_trackTelemetry) {
        m_avgCycleTime+=(m_stpCycle.read()-m_avgCycleTime)/m_nbIterations;
        m_avgIdleTime+=(stpIdle.read()-m_avgIdleTime)/m_nbIterations;
        m_avgNewPhotonsTracedCycle+=int((genPhotonsCycle-m_avgNewPhotonsTracedCycle))/m_nbIterations;

        if (m_nbPhotonsTracedCycle==0) {
            m_avgTimePerPhoton=0.0;
        }
        else {
            m_avgTimePerPhoton=(m_stpTracing.read()/m_nbPhotonsTracedCycle);
        }

        m_avgCellsPerPhotonTotal+=int((m_avgCellsPerPhotonCycle-m_avgCellsPerPhotonTotal))/m_nbIterations;
        //cout << "P " << m_myProcessRank << " m_avgCellsPerPhotonTotal=" << m_avgCellsPerPhotonTotal << "\n";

        m_avgProcessEfficiency+=(m_stpTracing.read()/m_stpCycle.read()-m_avgProcessEfficiency)/m_nbIterations;

        totalNbPhotonsLeftGlobalLastCycle=totalNbPhotonsLeftGlobal;
        totalNbPhotonsLeftLocal=toGenerateCellPhotons+toGenerateWallPhotons+photonStack.size();

        MPI_Reduce(&m_avgCellsPerPhotonTotal, &m_avgCellsPerPhotonGlobal, 1,
                   Common::MPIStructDef::getMPIType(&m_avgCellsPerPhotonTotal), MPI_SUM, 0, m_comm);
        MPI_Reduce(&m_avgIdleTime, &m_avgIdleTimeGlobal, 1,
                   Common::MPIStructDef::getMPIType(&m_avgIdleTime), MPI_SUM, 0, m_comm);
        MPI_Reduce(&m_avgTimePerPhoton, &m_avgTimePerPhotonGlobal, 1,
                   Common::MPIStructDef::getMPIType(&m_avgTimePerPhoton), MPI_SUM, 0, m_comm);
        MPI_Reduce(&m_avgProcessEfficiency, &m_avgProcessEfficiencyGlobal, 1,
                   Common::MPIStructDef::getMPIType(&m_avgProcessEfficiency), MPI_SUM, 0, m_comm);
        MPI_Reduce(&m_avgNewPhotonsTracedCycle, &m_avgNewPhotonsTracedCycleGlobal, 1,
                   Common::MPIStructDef::getMPIType(&m_avgNewPhotonsTracedCycle), MPI_SUM, 0, m_comm);
        MPI_Reduce(&totalNbPhotonsLeftLocal, &maxNbPhotonsLeftGlobal, 1,
                   Common::MPIStructDef::getMPIType(&totalNbPhotonsLeftLocal), MPI_MAX, 0, m_comm);
        MPI_Reduce(&totalNbPhotonsLeftLocal, &totalNbPhotonsLeftGlobal, 1,
                                                     Common::MPIStructDef::getMPIType(&totalNbPhotonsLeftLocal), MPI_SUM, 0, m_comm);

        //    MPI_Reduce(&m_avgTimePerPhoton, &m_maxAvgTimePerPhoton, 1,
        //          Common::MPIStructDef::getMPIType(&m_avgTimePerPhoton), MPI_MAX, 0, m_comm);

        m_avgPhotonsCycle+=int((totalNbPhotonsLeftGlobalLastCycle-totalNbPhotonsLeftGlobal)-m_avgPhotonsCycle)/m_nbIterations;

        if (m_myProcessRank==0) {
            m_avgCellsPerPhotonGlobal/=m_nbProcesses;
            m_avgIdleTimeGlobal/=m_nbProcesses;
            m_avgTimePerPhotonGlobal/=m_nbProcesses;
            m_avgProcessEfficiencyGlobal/=m_nbProcesses;
        }


//        approxSecondsLeft=int(CFreal(maxNbPhotonsLeftGlobal)*m_avgTimePerPhotonGlobal);
        approxSecondsLeft=int((m_avgCycleTime/CFreal(m_avgPhotonsCycle))*CFreal(totalNbPhotonsLeftGlobal));


        //    printPerformanceMetrics();


        CFLog(VERBOSE,"RadiativeTransferMonteCarloHSNB::computePhotons => Received "<<photonStack.size()<< " photons and has generated "<< genPhotonsCycle <<" photons\n");
        CFLog(VERBOSE,"RadiativeTransferMonteCarloHSNB::computePhotons => Number of photons left for this process: "<< totalNbPhotonsLeftLocal <<"\n");
        CFLog(VERBOSE,"RadiativeTransferMonteCarloHSNB::computePhotons => Total number of photons left: "<< totalNbPhotonsLeftGlobal <<"\n");

        CFLog(VERBOSE,"RadiativeTransferMonteCarloHSNB::computePhotons => Max number of photons left / single process: "<< maxNbPhotonsLeftGlobal <<"\n");
        CFLog(VERBOSE,"RadiativeTransferMonteCarloHSNB::computePhotons => Approx. extrapolated net ray-tracing time left: "  << int(approxSecondsLeft/3600)<< " h "
              << int((approxSecondsLeft/60)%60) << " m " << int(approxSecondsLeft%60)<<  " s \n");
        CFLog(VERBOSE,"RadiativeTransferMonteCarloHSNB::computePhotons => Photon generation cycle took "<< m_stpCycle.read() <<" seconds. \n");
    }

    /////////////////////////////////

    //Clear all buffers in the synchronizer
    m_paramSynchronizer.reset();

    //std::cout << "P" << m_myProcessRank << " MEM: sizeof(m_paramSynchronizer)=" << sizeof(m_paramSynchronizer) << "\n";
  }
  delete progressBar;

  printPerformanceMetrics();
  m_paramSynchronizer.printPerformanceMetrics();

  CFLog(INFO,"\n RadiativeTransferMonteCarloHSNB::computePhotons() => Raytracing took "<<s.readTimeHMS().str()<<'\n');



}

/////////////////////////////////////////////////////////////////////////////

template<class PARTICLE_TRACKING>
void RadiativeTransferMonteCarloHSNB<PARTICLE_TRACKING>::printPhoton(const Photon& photon)
{
  using namespace std;
  using namespace COOLFluiD::Framework;
  using namespace COOLFluiD::MathTools;
  using namespace COOLFluiD::Common;
  using namespace COOLFluiD::Numerics::FiniteVolume;
  using namespace COOLFluiD::LagrangianSolver;

  CFLog(INFO,
    "Photon data: "      <<'\n'
    <<"  cellID: "          <<photon.commonData.cellID<<'\n'
    <<"  currentPoint: "    <<photon.commonData.currentPoint[0]<<' '
    <<photon.commonData.currentPoint[1]<<' '
    <<photon.commonData.currentPoint[2]<<'\n'
    <<"  direction: "       <<photon.commonData.direction[0]<<' '
    <<photon.commonData.direction[1]<<' '
    <<photon.commonData.direction[2]<<'\n'
    <<"  wavelength: "      <<photon.userData.wavelength<<'\n'
    <<"  KS: "              <<photon.userData.KS<<'\n'
    <<"  Energy Fraction: " <<photon.userData.energyFraction<<'\n'
    <<"  Energy Residuum: " <<photon.userData.energyResiduum<<'\n'
    <<"  Previous process: " << photon.userData.prevProcessId<<'\n'
      );
}

/////////////////////////////////////////////////////////////////////////////

template<class PARTICLE_TRACKING>
void RadiativeTransferMonteCarloHSNB<PARTICLE_TRACKING>::printPerformanceMetrics()
{

    CFLog(INFO, "RadiativeTransferMonteCarloHSNB::printPerformanceMetrics => Average computation time / cycle=" << m_avgCycleTime << " seconds. \n");
    CFLog(INFO, "RadiativeTransferMonteCarloHSNB::printPerformanceMetrics => Average number of photons completely traced / cycle="<< m_avgPhotonsCycle << "\n");
    CFLog(INFO, "RadiativeTransferMonteCarloHSNB::printPerformanceMetrics => Local average idle time=" << m_avgIdleTime << " seconds. \n");
    CFLog(INFO, "RadiativeTransferMonteCarloHSNB::printPerformanceMetrics => Global average idle time=" << m_avgIdleTimeGlobal << " seconds. \n");
    CFLog(INFO, "RadiativeTransferMonteCarloHSNB::printPerformanceMetrics => Local average time / photon=" << m_avgTimePerPhoton << " seconds. \n");
    CFLog(INFO, "RadiativeTransferMonteCarloHSNB::printPerformanceMetrics => Global average time / photon=" << m_avgTimePerPhotonGlobal << " seconds. \n");
    CFLog(INFO, "RadiativeTransferMonteCarloHSNB::printPerformanceMetrics => Average cells per photon / this cycle=" << m_avgCellsPerPhotonCycle << " cells. \n");
    CFLog(INFO, "RadiativeTransferMonteCarloHSNB::printPerformanceMetrics => Local average cells per photon / total=" << m_avgCellsPerPhotonTotal << " cells. \n");
    CFLog(INFO, "RadiativeTransferMonteCarloHSNB::printPerformanceMetrics => Global average cells per photon / total=" << m_avgCellsPerPhotonGlobal << " cells. \n");
    CFLog(INFO, "RadiativeTransferMonteCarloHSNB::printPerformanceMetrics => Avg global efficiency=" << m_avgProcessEfficiencyGlobal*100 << "%.  \n");


}

/////////////////////////////////////////////////////////////////////////////

template<class PARTICLE_TRACKING>
void RadiativeTransferMonteCarloHSNB<PARTICLE_TRACKING>::MonteCarlo()
{
  using namespace std;
  using namespace COOLFluiD::Framework;
  using namespace COOLFluiD::MathTools;
  using namespace COOLFluiD::Common;
  using namespace COOLFluiD::Numerics::FiniteVolume;
  using namespace COOLFluiD::LagrangianSolver;

  const CFuint nbLoops = m_radiation->getNumberLoops();

  //FIXME: reset persistent cell and face iterators
  m_iphoton_cell_fix=0, m_istate_cell_fix=0;
  m_iphoton_face_fix=0, m_igState_face_fix=0;

  for(CFuint i=0; i< nbLoops; ++i){
    m_radiation->setupWavStride(i);
    getTotalEnergy();


    computePhotons();

  }
}

/////////////////////////////////////////////////////////////////////////////

template<class PARTICLE_TRACKING>
CFuint RadiativeTransferMonteCarloHSNB<PARTICLE_TRACKING>::rayTracingFullAbsorption(Photon& beam)
{
    using namespace std;
    using namespace COOLFluiD::Framework;
    using namespace COOLFluiD::MathTools;
    using namespace COOLFluiD::Common;
    using namespace COOLFluiD::Numerics::FiniteVolume;
    using namespace COOLFluiD::LagrangianSolver;


    m_nbPhotonsTracedCycle++;
    CFuint nbCrossedCells = 0;
    CFuint nbIter = 0;
    CFint exitCellID = 0;
    CFint exitFaceID = 0;
    CFint currentCellID = 0;

    CFLog(DEBUG_MED, "RadiativeTransferMonteCarlo::rayTracing() => START\n");

    m_lagrangianSolver.newParticle(beam);

    CFLog(DEBUG_MED, "RadiativeTransferMonteCarlo::rayTracing() => particle ID: "<<beam.commonData.cellID<< "\n");


    HSNBPhotonData &beamData = m_lagrangianSolver.getUserDataPtr();
    m_curTraceSet.headerData=&beamData;
    exitCellID=m_lagrangianSolver.getExitCellID();


//    CFLog(INFO, "NEW TRACE \n");


    while(nbIter <= m_maxVisitedCellsTrace){
        //CFLog(INFO, "New step!\n");
        nbIter++;
        m_debugTraceUID++;

        currentCellID = exitCellID;

        m_lagrangianSolver.trackingStep();
        exitFaceID=m_lagrangianSolver.getExitFaceID();
        exitCellID=m_lagrangianSolver.getExitCellID();
        m_nbCrossedCellsCycle++;
        //TODO is done in addStateParams
        //m_curTraceSet.headerData->nbCrossedCells++;
        beam.userData.nbCrossedCells++;


//        CFLog(INFO, "m_nbCrossedCellsCycle= " << m_nbCrossedCellsCycle<< ", m_maxVisitedCellsCycle="<< m_maxVisitedCellsCycle << "\n");

//        if ((m_nbCrossedCellsCycle>=m_maxVisitedCellsCycle) && (m_enforceSync==false)) {
//            m_enforceSync=true;
//        }

        checkSyncForceConditions();

        if(exitFaceID>=0){

            const CFreal stepDistance=m_lagrangianSolver.getStepDistance();
            RealVector null;

            m_baseRadiatorPtr = m_radiation->getCellDistPtr(exitCellID)->getRadiatorPtr();
            m_HSNBRadiator=SharedPtr<HSNBRadiator>(static_cast<HSNBRadiator*>(&(*m_baseRadiatorPtr)));

            //All HSNB parameters needed for the computation of absorption along the current photon's trace are added
            m_HSNBRadiator->addStateParams(m_curTraceSet,currentCellID,stepDistance);

            //CFLog(INFO, "Absorption IN!\n");
            const CFreal cellK= m_HSNBRadiator->getAbsorption(m_curTraceSet);
                        //CFLog(INFO, "Absorption OUT!\n");

//            beamData.KS -= stepDistance*cellK;
//            CFLog(INFO, "KS=" << beamData.KS << ", optThick=" << cellK << "\n" );
//            beamData.KS -= cellK;
//            CFLog(INFO, "KS NEW=" << beamData.KS <<  "\n" );

            if((beamData.KS-cellK) <= 0.){ // photon absorbed by a cell
//                CFLog(INFO, "Full absorption \n" );
                const CFuint gEndId = currentCellID;
                const CFreal energyFraction = beamData.energyFraction;

                cf_assert(gEndId < m_stateInRadPowers.size());
                m_stateInRadPowers[gEndId]+=energyFraction;


                return currentCellID;
            }

            const CFuint faceType = m_lagrangianSolver.getFaceType(exitFaceID);

            if ( faceType == ParticleTracking::WALL_FACE){
                //CFLog(INFO,"HERE WALL !!\n");
                CommonData beam2;
                m_lagrangianSolver.getCommonData(beam2);
                for(CFuint i=0; i < m_dim2; ++i){
                    m_entryDirection[i]= beam2.direction[i];
                }

                m_lagrangianSolver.getExitPoint(m_position);

                //CFuint stateID = m_lagrangianSolver.getWallGhotsStateId(exitFaceID);
                const CFuint ghostStateID = m_lagrangianSolver.getWallGhotsStateId(exitFaceID);

                const CFreal wallK = m_radiation->getWallDistPtr(ghostStateID)
                        ->getRadiatorPtr()->getAbsorption( beamData.wavelength, m_entryDirection );

                m_lagrangianSolver.getNormals(exitFaceID, m_position, m_normal);

                const CFreal reflectionProbability =  m_rand.uniformRand();
                CFLog(DEBUG_MIN, "reflectionProbability[" << reflectionProbability << "] <= wallK["
                      << wallK << "]\n");

                if (reflectionProbability <= wallK){ // the photon is absorbed by the wall

                    const CFuint ghostStateID = m_lagrangianSolver.getWallGhotsStateId(exitFaceID);
                    m_ghostStateInRadPowers[ghostStateID] += beamData.energyFraction;
                    CFLog(DEBUG_MIN, "Rad power in ghostStateID[" << ghostStateID << "] = " <<
                          m_ghostStateInRadPowers[ghostStateID] << "\n");
                    //foundEntity = true;
                    m_totalWallAbsorption+=beamData.energyFraction;

                    return exitFaceID;
                }
                else {
                    m_radiation->getWallDistPtr(ghostStateID)->getReflectorPtr()->getRandomDirection
                            (beamData.wavelength, m_exitDirection, m_entryDirection, m_normal);
                    m_lagrangianSolver.newDirection( m_exitDirection );


                    CFLog(DEBUG_MED, "Particle reflected with Entry Direction[" << m_entryDirection
                          << "], Normal[ " << m_normal << "], Exit direction [" << m_exitDirection << "]\n");
                }
            }


            if (faceType == ParticleTracking::BOUNDARY_FACE ){
                //CFLog(INFO,"HERE BOUNDARY !!\n");
                //entity = DISAPPEARED;
                //foundEntity = true;
                return 0;
            }

            if(faceType == ParticleTracking::COMP_DOMAIN_FACE){
                beamData.nbCrossedCells=m_curTraceSet.trace->m_nbCellsCrossed;
                beamData.targetProcessId=m_lagrangianSolver.getTargetProcess(exitFaceID);

                m_lagrangianSolver.bufferCommitParticle(exitFaceID);

                m_debugPhotonsComitted++;
                //m_curTraceSet.headerData->targetProcessId=m_lagrangianSolver.getTargetProcess(exitFaceID);

                m_curTraceSet.headerData->targetProcessId=m_lagrangianSolver.getTargetProcess(exitFaceID);
                m_paramSynchronizer.bufferCommitParticle(m_curTraceSet);


//                if ((m_paramSynchronizer.totalSendBufferByteSize()>=m_traceBufferMaxByteSize) && (m_enforceSync==false)) {
//                    m_enforceSync=true;
 //               }

                checkSyncForceConditions();


                return 0;
            }

            nbCrossedCells++;
        }
        else{
            CFLog(DEBUG_MED, "RadiativeTransferMonteCarlo::rayTracing() => enter negligible\n");
            //entity = NEGLIGIBLE;


            return 0;
        }
    }
    CFLog(DEBUG_MED, "RadiativeTransferMonteCarlo::rayTracing() => Max number of steps reached! \n");
    return 0;
}


/////////////////////////////////////////////////////////////////////////////

template<class PARTICLE_TRACKING>
CFuint RadiativeTransferMonteCarloHSNB<PARTICLE_TRACKING>::rayTracing(Photon& beam)
{
    using namespace std;
    using namespace COOLFluiD::Framework;
    using namespace COOLFluiD::MathTools;
    using namespace COOLFluiD::Common;
    using namespace COOLFluiD::Numerics::FiniteVolume;
    using namespace COOLFluiD::LagrangianSolver;

    CFuint faceType;
    CFreal stepDistance;
    m_nbPhotonsTracedCycle++;
    CFuint nbCrossedCells = 0;
    CFuint nbIter = beam.userData.nbCrossedCells;
    CFint exitCellID = 0;
    CFint exitFaceID = 0;
    CFint currentCellID = 0;

    CFLog(DEBUG_MED, "RadiativeTransferMonteCarloHSNB::rayTracing() => START\n");

    m_lagrangianSolver.newParticle(beam);

    CFLog(DEBUG_MED, "RadiativeTransferMonteCarloHSNB::rayTracing() => particle ID: "<< beam.commonData.cellID<< "\n");

    HSNBPhotonData &beamData = m_lagrangianSolver.getUserDataPtr();
    //Associate the current photon with m_curTraceSet interface
    m_curTraceSet.headerData=&beamData;

    exitCellID=m_lagrangianSolver.getExitCellID();

    cf_assert(!isnan(beam.userData.energyResiduum));

    while(nbIter <= m_maxVisitedCellsTrace){
        //CFLog(INFO, "New step!\n");
        nbIter++;
        CFLog(DEBUG_MED, "RadiativeTransferMonteCarloHSNB::rayTracing => UID=" <<  m_debugTraceUID << "\n");
        m_debugTraceUID++;

        currentCellID = exitCellID;

        m_lagrangianSolver.trackingStep();
        exitFaceID=m_lagrangianSolver.getExitFaceID();
        exitCellID=m_lagrangianSolver.getExitCellID();
        m_nbCrossedCellsCycle++;
        //TODO is done in addStateParams
        //m_curTraceSet.headerData->nbCrossedCells++;
        beam.userData.nbCrossedCells++;

        checkSyncForceConditions();

        if(exitFaceID>=0){

            //Careful: The raytracer computes distances in [m] but the HSNB radiator computes absorption in [cm]
            //Thus we have to convert stepDistance by multiplying with a factor of 100 once it's added to the trace.
            stepDistance=m_lagrangianSolver.getStepDistance();

            m_baseRadiatorPtr = m_radiation->getCellDistPtr(exitCellID)->getRadiatorPtr();
            m_HSNBRadiator=SharedPtr<HSNBRadiator>(static_cast<HSNBRadiator*>(&(*m_baseRadiatorPtr)));

            //std::cout <<"P" <<m_myProcessRank << ": "<< m_curTraceSet.trace->m_thinDiatomics.size() << " thin Diatomics size \n";


            cf_assert(currentCellID < m_stateInRadPowers.size());
            CFLog(DEBUG_MED, "RadiativeTransferMonteCarloHSNB::rayTracing => distance=" <<  std::setprecision(16) << stepDistance*100 << ", state="<< currentCellID << "\n");

            //Add tracking step state parameters:
            m_HSNBRadiator->addStateParams(m_curTraceSet,currentCellID,stepDistance);

            CFreal tempRes=m_HSNBRadiator->computeAbsorbedEnergy(m_curTraceSet,currentCellID);

            CFLog(DEBUG_MED, "RadiativeTransferMonteCarloHSNB::rayTracing => absorbed energy "<< tempRes << ", m_stateInRadPowers["<< currentCellID << "]="<< m_stateInRadPowers[currentCellID] << "\n");
//            CFLog(INFO, "sizeof(m_HSNBRadiator)=" << sizeof(*m_HSNBRadiator) << "\n");
            m_stateInRadPowers[currentCellID]+=tempRes;
            CFLog(DEBUG_MED, "RadiativeTransferMonteCarloHSNB::rayTracing => m_stateInRadPowers["<< currentCellID << "]="<< m_stateInRadPowers[currentCellID] << "\n");


            cf_assert(!isnan(m_stateInRadPowers[currentCellID]));


            CFLog(DEBUG_MAX, "RadiativeTransferMonteCarloHSNB::rayTracing => computeAbsorbedEnergy=" << m_stateInRadPowers[currentCellID] << ", photon.transm="
                  << m_curTraceSet.headerData->transm << ", photon.energyResiduum" << m_curTraceSet.headerData->energyResiduum << "\n" );
            CFLog(DEBUG_MAX, "RadiativeTransferMonteCarloHSNB::rayTracing => m_stateInRadPowers[" << currentCellID << "]=" << m_stateInRadPowers[currentCellID]
                  << ", energyResiduum=" << m_curTraceSet.headerData->energyResiduum <<  ", rayEnergy=" << m_curTraceSet.headerData->energyFraction << " \n");


            //If no more energy is sensed we're done, the photon is fully absorbed and we stop the tracing
            if ((m_curTraceSet.headerData->transm<m_tolerance) || (m_curTraceSet.headerData->energyResiduum==0.0)) {
                CFLog(DEBUG_MED, "RadiativeTransferMonteCarloHSNB::rayTracing => Absorbed with residuum="<< m_curTraceSet.headerData->energyResiduum <<" \n");

                //Free memory if possible
                m_curTraceSet.reset();
                return currentCellID;
            }


            faceType = m_lagrangianSolver.getFaceType(exitFaceID);


            if ( faceType == ParticleTracking::WALL_FACE){
                //CFLog(INFO,"HERE WALL !!\n");
                CommonData beam2;
                m_lagrangianSolver.getCommonData(beam2);

                for(CFuint i=0; i < m_dim2; ++i){
                    m_entryDirection[i]= beam2.direction[i];
                }

                m_lagrangianSolver.getExitPoint(m_position);

                //CFuint stateID = m_lagrangianSolver.getWallGhotsStateId(exitFaceID);
                const CFuint ghostStateID = m_lagrangianSolver.getWallGhotsStateId(exitFaceID);

                // Use absorptivity / emissivity?
                const CFreal wallAbs=m_radiation->getWallDistPtr(ghostStateID)
                        ->getRadiatorPtr()->getAbsorption(beamData.wavelength, m_entryDirection);

//              const CFreal wallEmis=m_radiation->getWallDistPtr(ghostStateID)
//                                        ->getRadiatorPtr()->getEmissivity();

                //CFLog(DEBUG_MAX, "Before reflection: " <<m_curTraceSet.headerData->transm << "\n");

                const CFreal transp=m_curTraceSet.headerData->transm*(1-wallAbs);

                CFreal m_tempAbsorbedEnergy=m_curTraceSet.headerData->energyFraction*(m_curTraceSet.headerData->transm-transp);
                m_ghostStateInRadPowers[ghostStateID] +=  m_tempAbsorbedEnergy;

                cf_assert(!isnan(m_ghostStateInRadPowers[ghostStateID]));

                //Update photon energy fraction after wall absorption:
                m_curTraceSet.headerData->energyResiduum-= m_tempAbsorbedEnergy;
                m_curTraceSet.headerData->transm=transp;
                m_curTraceSet.headerData->wallTransmissivity*=(1-wallAbs);

                m_lagrangianSolver.getNormals(exitFaceID, m_position, m_normal);

                CFLog(DEBUG_MIN, "Rad power in ghostStateID[" << ghostStateID << "] = " <<
                      m_ghostStateInRadPowers[ghostStateID] << "\n");
                //foundEntity = true;

                //Check whether photon is fully absorbed by the wall, if not:  reflect it and continue trace
                if (m_curTraceSet.headerData->transm>=m_tolerance) {

                    //If the photon has not been fully absorbed: Reflect it
                    m_radiation->getWallDistPtr(ghostStateID)->getReflectorPtr()->getRandomDirection
                            (beamData.wavelength, m_exitDirection, m_entryDirection, m_normal);
                    m_lagrangianSolver.newDirection( m_exitDirection );

                    CFLog(DEBUG_MED, "Particle reflected with Entry Direction[" << m_entryDirection
                          << "], Normal[ " << m_normal << "], Exit direction [" << m_exitDirection << "]\n");
                }
                else {
                    CFLog(DEBUG_MED, "Full wall absorption: wallEmis=" << wallAbs << ", transp=" << transp <<  ", m_tempAbsorbedEnergy=" << m_tempAbsorbedEnergy
                          <<", residuum= " << m_curTraceSet.headerData->energyResiduum << "\n");
                    //Free memory if possible
                    m_curTraceSet.reset();
                    return currentCellID;
                }
            }


            if (faceType == ParticleTracking::BOUNDARY_FACE ){
                CFLog(DEBUG_MAX,"HERE BOUNDARY !!\n");
                m_totalWallAbsorption+=beamData.energyResiduum;

                //entity = DISAPPEARED;
                //foundEntity = true;
                return 0;
            }

            if(faceType == ParticleTracking::COMP_DOMAIN_FACE){
                CFLog(DEBUG_MAX,"HERE DOMAIN FACE!!\n");
                m_nbTracesCommittedCycle++;

                beamData.nbCrossedCells=m_curTraceSet.trace->m_nbCellsCrossed;
                beamData.targetProcessId=m_lagrangianSolver.getTargetProcess(exitFaceID);

                m_lagrangianSolver.bufferCommitParticle(exitFaceID);


                m_debugPhotonsComitted++;
                //m_curTraceSet.headerData->targetProcessId=m_lagrangianSolver.getTargetProcess(exitFaceID);

                m_curTraceSet.headerData->targetProcessId=m_lagrangianSolver.getTargetProcess(exitFaceID);
                m_paramSynchronizer.bufferCommitParticle(m_curTraceSet);

                checkSyncForceConditions();


                return 0;
            }

            nbCrossedCells++;
        }
        else{
            CFLog(DEBUG_MIN, "RadiativeTransferMonteCarloHSNB::rayTracing() => enter negligible\n");
            //entity = NEGLIGIBLE;


            return 0;
        }
    }
    CFLog(DEBUG_MIN, "RadiativeTransferMonteCarloHSNB::rayTracing() => Max number of steps reached! \n");
    return 0;
}

/////////////////////////////////////////////////////////////////////////////

template<class PARTICLE_TRACKING>
void RadiativeTransferMonteCarloHSNB<PARTICLE_TRACKING>::computeHeatFlux()
{
  using namespace std;
  using namespace COOLFluiD::Framework;
  using namespace COOLFluiD::MathTools;
  using namespace COOLFluiD::Common;
  using namespace COOLFluiD::Numerics::FiniteVolume;
  using namespace COOLFluiD::LagrangianSolver;

  CFLog(VERBOSE, "RadiativeTransferMonteCarloHSNB computeHeatFlux() START\n");


  DataHandle<State*, GLOBAL> states = socket_states.getDataHandle();
  DataHandle<CFreal> volumes = socket_volumes.getDataHandle();
  DataHandle<CFreal> gasRadiativeHeatSource = socket_qrad.getDataHandle();

  m_stateInRadPowers.sincronizeAssign();
  m_stateRadPower.sincronizeAssign();
  //fill qrad socket (loop states)
//  CFreal totalNetAb=0.0;
//  for(CFuint istate = 0; istate < states.size(); ++istate){
//      totalNetAb+=m_stateInRadPowers[istate];
//      cout << "m_stateInRadPowers[" << istate << "]=" << m_stateInRadPowers[istate] << "\n";
//  }
//  for(CFuint istate = 0; istate < states.size(); ++istate){
//      //cout << "m_stateRadPowers[" << istate << "]=" << m_stateRadPower[istate] << "\n";
//  }
//  CFreal totalNetEm=0.0;
//  for(CFuint istate = 0; istate < states.size(); ++istate){
//      totalNetEm+=m_stateNetRadPower[istate];
//      cout << "m_stateNetRadPowers[" << istate << "]=" << m_stateNetRadPower[istate] << "\n";
//  }

  int active=0;
//  if (m_myProcessRank == 0) {
//      active = 0;

//        for(CFuint istate = 0; istate < states.size(); ++istate){
//           cout << "m_stateInRadPowers[" << active << "]=" << m_stateInRadPowers[istate] << "\n";
//           active++;
//        }


//      MPI_Send(&active, 1, MPIStructDef::getMPIType(&active), m_myProcessRank+1, 0, m_comm);
//  }
//  else {
//      MPI_Recv(&active, 1, MPIStructDef::getMPIType(&active), m_myProcessRank-1, 0, m_comm, MPI_STATUS_IGNORE);

//      for(CFuint istate = 0; istate < states.size(); ++istate){
//         cout << "m_stateInRadPowers[" << active << "]=" << m_stateInRadPowers[istate] << "\n";
//         active++;
//      }

//      if (m_myProcessRank!=m_nbProcesses) {
//          MPI_Send(&active, 1, MPIStructDef::getMPIType(&active), m_myProcessRank+1, 0, m_comm);
//      }
//  }
//  std::cout << "Emitted power: " << totalNetEm << ", cell absorption="<< totalNetAb << ", wall absorption" << m_totalWallAbsorption  << ", crossing fourty="<< m_nbCrossingForty<< "\n";


  for(CFuint istate = 0; istate < states.size(); ++istate){
      //if (states[istate]->isParUpdatable() ){ //is it necessary ?
      CFreal volume = volumes[istate];
      volume *= (m_isAxi ? 6.283185307179586*(*states[istate]).getCoordinates()[YY] : 1.);
      //gasRadiativeHeatSource[istate] = (- m_stateInRadPowers[istate] + m_stateRadPower[istate]) / volume ;
//      gasRadiativeHeatSource[istate] = (1.-m_relaxationFactor)*gasRadiativeHeatSource[istate] +
//                                       m_relaxationFactor*(m_stateInRadPowers[istate]-m_stateRadPower[istate])/volume;

      gasRadiativeHeatSource[istate] = (1.-m_relaxationFactor)*gasRadiativeHeatSource[istate] +
                                       m_relaxationFactor*(m_stateInRadPowers[istate]-m_stateNetRadPower[istate])/volume;
//      cout<< m_stateInRadPowers[istate] <<" hsnbstate("<<istate<<")="
//              << gasRadiativeHeatSource[istate]<< " m_stateNetRadPower[istate]=" << m_stateNetRadPower[istate] << endl;

      //}
  }

  //Fill qradFluxWall socket (loop Trs idxs for each trs)
  //Needs to be consistent with AeroForces definition
  DataHandle<CFreal> faceAreas   = socket_faceAreas.getDataHandle();
  DataHandle<CFreal> faceCenters = socket_faceCenters.getDataHandle();
  DataHandle<CFreal> wallRadiativeHeatSource = socket_qradFluxWall.getDataHandle();

  vector<string> wallTrsNames;
  m_radiation->getWallTRSnames(wallTrsNames);
  const CFuint nbWallTrs = wallTrsNames.size();
  CFLog(DEBUG_MED, "RadiativeTransferMonteCarloHSNB::computeHeatFlux() => nbWallTrs = " << nbWallTrs << "\n");

  CFuint ff = 0;
  for(CFuint i=0; i<nbWallTrs; ++i){
    Framework::FaceTrsGeoBuilder::GeoData& facesData = m_wallFaceBuilder.getDataGE();
    SafePtr<TopologicalRegionSet> wallFaces =
      MeshDataStack::getActive()->getTrs(wallTrsNames[i]);

    facesData.trs = wallFaces;
    const CFuint nbFacesWall = wallFaces->getLocalNbGeoEnts();
    //cout<<"CPU "<< Common::PE::GetPE().GetRank()<<": trying TRS " << wallTrsNames[i]<<" whith "<<nbFacesWall<< " faces " <<endl;

    for(CFuint f=0; f<nbFacesWall; ++f, ++ff){
      facesData.idx = f;
      Framework::GeometricEntity *const face = m_wallFaceBuilder.buildGE();
      const CFuint faceGeoID = face->getID();
      const CFuint faceGhostStateID = face->getState(1)->getLocalID();
      CFreal area = faceAreas[faceGeoID];
      area *= (m_isAxi) ? 6.283185307179586*faceCenters[faceGeoID*DIM_2D + YY] : 1.;

//      CFLog(INFO, "face[" << f << "] =>" <<  m_ghostStateInRadPowers[faceGhostStateID] << ", "
//        << m_ghostStateRadPower[faceGhostStateID] << "\n");
      wallRadiativeHeatSource[ff] = (-m_ghostStateInRadPowers[faceGhostStateID]
                                     +m_ghostStateRadPower[faceGhostStateID] ) / area;


      m_wallFaceBuilder.releaseGE();
    }
    //cout<<"CPU "<< Common::PE::GetPE().GetRank()<<" done " <<endl;
  }


  m_postProcess->runPostProcess(gasRadiativeHeatSource);

  // change the sign of the heat source
  for (CFuint i = 0; i < gasRadiativeHeatSource.size(); ++i) {
   gasRadiativeHeatSource[i] *= -1.;
  }



  CFLog(VERBOSE, "RadiativeTransferMonteCarloHSNB computeHeatFlux() END\n");


}

/////////////////////////////////////////////////////////////////////////////

  } // namespace RadiativeTransfer

} // namespace COOLFluiD

/////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_RadiativeTransfer_RadiativeTransferMonteCarloHSNBFVMCC_hh
