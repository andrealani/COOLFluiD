#ifndef COOLFluiD_RadiativeTransfer_RadiativeTransferMonteCarloFVMCC_hh
#define COOLFluiD_RadiativeTransfer_RadiativeTransferMonteCarloFVMCC_hh

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
#include "LagrangianSolver/ParallelVector/ParallelVector.hh"

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

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  
  namespace Numerics {
    namespace FiniteVolume {
      class CellCenterFVMData;
    }
  }

 namespace RadiativeTransfer {
   
//////////////////////////////////////////////////////////////////////////////

struct PhotonData{
    CFreal KS;
    CFreal energyFraction;
    CFreal wavelength;
};

typedef LagrangianSolver::Particle<PhotonData> Photon;

//////////////////////////////////////////////////////////////////////////////

template<class PARTICLE_TRACKING>
class RadiativeTransferMonteCarlo : public Framework::DataProcessingCom                 
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
   RadiativeTransferMonteCarlo(const std::string& name);

  /**
   * Default destructor
   */
   ~RadiativeTransferMonteCarlo();

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
   * ray tracing
   */
  CFuint rayTracing(Photon& photon);
  
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
  
  /// Compute cell rays
  void computePhotons();
  
  bool getCellPhotonData(Photon& ray);

private: 

  LagrangianSolver::LagrangianSolver<PhotonData, PARTICLE_TRACKING> m_lagrangianSolver;
  
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
  
  /// number of dimension
  CFuint m_dim;

  ///my process
  CFuint m_myProcessRank;

  ///number processes
  CFuint m_nbProcesses;

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

  /// number of rays that each element emits
  CFuint m_nbRaysElem;
  CFuint m_nbRaysCycle;
  
  /// maximum number of visited cells
  CFuint m_maxVisitedCells;
  
  /// True if it is an axisymmetric simulation
  bool m_isAxi;
  
  RandomNumberGenerator m_rand;

  CFuint m_sendBufferSize;

  RealVector m_ghostStateInRadPowers;
  
  CFuint m_dim2;

  void getTotalEnergy();

  void printPhoton(Photon photon);

  CFreal m_totalRadPower;

  

  LagrangianSolver::ParallelVector<CFreal> m_stateRadPower;
  LagrangianSolver::ParallelVector<CFreal> m_stateInRadPowers;
  
  vector<CFuint> m_nbPhotonsState;

  vector<CFreal> m_ghostStateRadPower;
  vector<CFuint> m_nbPhotonsGhostState;

  //FIXME: persistent memory for cell and face iterators
  CFuint  m_iphoton_cell_fix, m_istate_cell_fix;
  CFuint  m_iphoton_face_fix, m_igState_face_fix;

  CFreal m_relaxationFactor;

  bool getFacePhotonData(Photon &ray);
}; // end of class RadiativeTransferMonteCarlo

//////////////////////////////////////////////////////////////////////////////

template<class PARTICLE_TRACKING>
void RadiativeTransferMonteCarlo<PARTICLE_TRACKING>::defineConfigOptions
(Config::OptionList& options)
{
  using namespace std;
  using namespace COOLFluiD::Common;
  
  options.addConfigOption< CFuint >("numberOfRays","number of rays sent by each element.");
  options.addConfigOption< CFuint >("MaxNbVisitedCells","Maximum number of visited cells.");
  options.addConfigOption< bool >("Axi","True if it is an axisymmetric simulation.");
  options.addConfigOption< string >("PostProcessName","Name of the post process routine");
  options.addConfigOption< CFuint >("sendBufferSize","Size of the buffer for communication");
  options.addConfigOption< CFuint >("nbRaysCycle","Number of rays to emit before communication step");
  options.addConfigOption< CFreal >("relaxationFactor","Relaxation Factor");
}

//////////////////////////////////////////////////////////////////////////////

template<class PARTICLE_TRACKING>
RadiativeTransferMonteCarlo<PARTICLE_TRACKING>::RadiativeTransferMonteCarlo
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

  m_postProcessName = "PostProcessNull";

  setParameter("PostProcessName",&m_postProcessName);

  setParameter("numberOfRays",&m_nbRaysElem);

  m_maxVisitedCells = 10000;
  setParameter("MaxNbVisitedCells",&m_maxVisitedCells);

  m_isAxi = false;
  setParameter("Axi",&m_isAxi);

  m_sendBufferSize=10000;
  setParameter("sendBufferSize", &m_sendBufferSize);

  m_nbRaysCycle = m_sendBufferSize / 2;
  setParameter("nbRaysCycle", &m_nbRaysCycle);

  m_relaxationFactor = 1.;
  setParameter("relaxationFactor", &m_relaxationFactor);
}

/////////////////////////////////////////////////////////////////////////////

template<class PARTICLE_TRACKING>
RadiativeTransferMonteCarlo<PARTICLE_TRACKING>::~RadiativeTransferMonteCarlo()
{
}

/////////////////////////////////////////////////////////////////////////////

template<class PARTICLE_TRACKING>
void RadiativeTransferMonteCarlo<PARTICLE_TRACKING>::configure(Config::ConfigArgs& args)
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
vector<SafePtr<BaseDataSocketSink> > RadiativeTransferMonteCarlo<PARTICLE_TRACKING>::
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
vector<SafePtr<BaseDataSocketSource> > RadiativeTransferMonteCarlo<PARTICLE_TRACKING>::
providesSockets()
{
  std::vector<Common::SafePtr<Framework::BaseDataSocketSource> > result;
  result.push_back(&socket_qrad);
  result.push_back(&socket_qradFluxWall);
  return result;
}

/////////////////////////////////////////////////////////////////////////////

template<class PARTICLE_TRACKING>
void RadiativeTransferMonteCarlo<PARTICLE_TRACKING>::setup()
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
  SafePtr<TopologicalRegionSet> cells = MeshDataStack::getActive()->getTrs("InnerCells");
  cellData.trs = cells;

  // setting up the wall face builder
  m_wallFaceBuilder.setup();
  m_wallFaceBuilder.getGeoBuilder()->setDataSockets(socket_states, socket_gstates, socket_nodes);
  m_wallFaceBuilder.getDataGE().isBFace = true;

  const CFuint nCells = cellData.trs->getLocalNbGeoEnts();

  m_radiation->configureTRS();
  m_radiation->setupAxiFlag(m_isAxi);

  //MPI_Datatype Userdatatype, particleDatatype;
  MPIStruct Userdatatype;//, particleDatatype;

  PhotonData photonData;
  int counts[3] = {1,1,1};
  MPIStructDef::buildMPIStruct<CFreal,CFreal,CFreal>
          (&photonData.KS, &photonData.energyFraction, &photonData.wavelength, counts , Userdatatype);

  m_lagrangianSolver.setupParticleDatatype( Userdatatype.type );
 // particleDatatype.type = m_lagrangianSolver.getParticleDataType();

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
void RadiativeTransferMonteCarlo<PARTICLE_TRACKING>::executeOnTrs()
{
  using namespace std;
  using namespace COOLFluiD::Framework;
  using namespace COOLFluiD::MathTools;
  using namespace COOLFluiD::Common;
  using namespace COOLFluiD::LagrangianSolver;
  
  CFLog(VERBOSE, "RadiativeTransferMonteCarlo::executeOnTrs() START\n");

  // for testcases == 0, the following must be recomputed everytime
  Stopwatch<WallTime> s;
  s.restart();

  // re-initialize the heat fluxes to 0
  //m_stateInRadPowers=0.;
  //m_RHSGas.assign(m_RHSGas.size(), 0.);
  
  MonteCarlo();
  computeHeatFlux();
  
  CFLog(INFO, "MonteCarlo() took " << s << "s\n");
  
  CFLog(VERBOSE, "RadiativeTransferMonteCarlo::executeOnTrs() END\n");
}

/////////////////////////////////////////////////////////////////////////////

template<class PARTICLE_TRACKING>
void RadiativeTransferMonteCarlo<PARTICLE_TRACKING>::getTotalEnergy()
{
  using namespace std;
  using namespace COOLFluiD::Framework;
  using namespace COOLFluiD::MathTools;
  using namespace COOLFluiD::Common;
  using namespace COOLFluiD::Numerics::FiniteVolume;
  using namespace COOLFluiD::LagrangianSolver;
  
  CFuint nbStates = m_radiation->getNbStates();
  CFuint nbGhostStates = m_radiation->getNbGhostStates();
  
  m_stateRadPower.resize(nbStates);
  m_nbPhotonsState.resize(nbStates);
  
  m_ghostStateRadPower.resize(nbGhostStates);
  m_nbPhotonsGhostState.resize(nbGhostStates);
  
  //CFreal nbPhotons = nbStates ;//* m_nbRaysElem / m_radiation->getNumberLoops();
  //cout<<"nbPhotons: "<< nbPhotons <<endl;
  
  //Get the total radiative power of the Cells
  CFreal stateRadPower =  0;
  for(CFuint state=0; state<nbStates; ++state){
    if ( !m_radiation->isStateNull(state) ){
      //cout<<"is not ghost!"<<endl;
      m_stateRadPower[state]= m_radiation->getCellDistPtr(state)
	->getRadiatorPtr()->getSpectraLoopPower();
      stateRadPower += m_stateRadPower[state];
      //cout<<"state radPower: "<<m_stateRadPower[state]<<endl;
      //cout<<"THE OTER: m_axi Volume: "<<m_axiVolumes[i]<<endl;
      //m_stateRadPower[i] = cellK*sigma*pow(T,4.)*m_axiVolumes[i];
    }
    else{
      //cout<<"is ghost!"<<endl;
      m_stateRadPower[state] = .0;
    }
  }
  
  //Get the total radiative power of the Wall Faces
  CFreal gStateRadPower =  0;
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
  }
  //cout<<" walls rad power : "<< gStateRadPower<<endl;
  
  m_totalRadPower = stateRadPower + gStateRadPower;
  
  //cout<<"statePhotons"<<endl;
  for(CFuint i= 0; i< m_stateRadPower.size(); ++i){
    //CFuint test = CFuint(m_stateRadPower[i]/m_totalRadPower*nbPhotons);
    m_nbPhotonsState[i] =( m_stateRadPower[i] > 0. )? m_nbRaysElem : 0 ;
    //cout<<"nbPhotons: "<<m_nbPhotonsState[i]<< ' '<<m_stateRadPower[i] <<endl;
  }
  
  //cout<<"ghostPhotons"<<endl;
  for(CFuint i= 0; i< m_ghostStateRadPower.size(); ++i){
    //CFuint test = CFuint(m_stateRadPower[i]/m_totalRadPower*nbPhotons);
    m_nbPhotonsGhostState[i] =( m_ghostStateRadPower[i] > 0.)? m_nbRaysElem : 0 ;
    //cout<<"nbPhotons: "<<m_nbPhotonsGhostState[i]<< ' '<<m_ghostStateRadPower[i] <<endl;
  }
}
  
/////////////////////////////////////////////////////////////////////////////

template<class PARTICLE_TRACKING>
bool RadiativeTransferMonteCarlo<PARTICLE_TRACKING>::getCellPhotonData(Photon &ray)
{
  using namespace std;
  using namespace COOLFluiD::Framework;
  using namespace COOLFluiD::MathTools;
  using namespace COOLFluiD::Common;
  using namespace COOLFluiD::Numerics::FiniteVolume;
  using namespace COOLFluiD::LagrangianSolver;
  
  //static CellTrsGeoBuilder::GeoData& cellData = m_cellBuilder.getDataGE();
  const CFuint nbStates = m_nbPhotonsState.size();
  
  for(;m_istate_cell_fix<nbStates; ++m_istate_cell_fix) {
    for(;m_iphoton_cell_fix<m_nbPhotonsState[ m_istate_cell_fix ];) {
      //        cout<<"m_nbPhotonsState[state]: "<<m_nbPhotonsState[state]<<' '<<";state: "<<state<<endl;
      //        cout<<"photon: "<<photon<<endl;
      
      //cellData.idx =  state;
      
      //cellData.idx = 50;
      //GeometricEntity *const cell = m_cellBuilder.buildGE();
      //if( cell->getState(0)->isParUpdatable() ){
      // Calculate the wavelength
      //cout<<"wavelength= "<<ray.userData.wavelength<<endl;
      
      //Get directions
      m_radiation->getCellDistPtr( m_istate_cell_fix )->
	getRadiatorPtr()->getRandomEmission(ray.userData.wavelength, m_direction );
      
      //cout<<"ray directions: ";
      for(CFuint ii=0; ii<m_dim2; ++ii){
	ray.commonData.direction[ii]=m_direction[ii];
      }
      //cout<<endl;
      
      //Get the beam max optical path Ks
      ray.userData.KS = - std::log( m_rand.uniformRand() );
      // ray.actualKS = 0;
      //  cout<<"getCellcenter"<<endl;
      //Get cell center
      static Framework::DataHandle<Framework::State*, Framework::GLOBAL> states
	= socket_states.getDataHandle();
      
      Node& baricenter = (*states[ m_istate_cell_fix ]).getCoordinates();
      
      for(CFuint i=0;i<m_dim;++i){
	//  cout<<baricenter[i]<<' ';
	ray.commonData.currentPoint[i]=baricenter[i];
      }
      for(CFuint i=m_dim;i<m_dim2;++i){
	//  cout<<baricenter[i]<<' ';
	ray.commonData.currentPoint[i] = 0.;
      }
      
      //Get the remaining information
      //const CFuint cellID = cell->getID();
      CFuint cellID = m_radiation->getCurrentCellStateID();
      ray.commonData.cellID = cellID;
      
      ray.userData.energyFraction= m_stateRadPower[ m_istate_cell_fix ]/CFreal(m_nbPhotonsState[ m_istate_cell_fix ]);
      
      //m_cellBuilder.releaseGE();
      ++m_iphoton_cell_fix;
      return true;
    }
    m_iphoton_cell_fix=0;
  }
  return false;
}
  
/////////////////////////////////////////////////////////////////////////////

template<class PARTICLE_TRACKING>
bool RadiativeTransferMonteCarlo<PARTICLE_TRACKING>::getFacePhotonData(Photon &ray)
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
      
      ++m_iphoton_face_fix;
      return true;
    }
    m_iphoton_face_fix=0;
  }
  return false;
}

/////////////////////////////////////////////////////////////////////////////

template<class PARTICLE_TRACKING>
void RadiativeTransferMonteCarlo<PARTICLE_TRACKING>::computePhotons()
{
  using namespace std;
  using namespace COOLFluiD::Framework;
  using namespace COOLFluiD::MathTools;
  using namespace COOLFluiD::Common;
  using namespace COOLFluiD::Numerics::FiniteVolume;
  using namespace COOLFluiD::LagrangianSolver;
    
  CFLog(DEBUG_MAX, "RadiativeTransferMonteCarlo::computeCellRays()\n");
  
  m_rand.seed(time(NULL)*(m_myProcessRank+1));


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
  CFuint totalnbPhotons = toGenerateWallPhotons + toGenerateCellPhotons;
  boost::progress_display* progressBar = NULL;
  if (m_myProcessRank == 0) progressBar = new boost::progress_display(totalnbPhotons);

  //CFuint toGeneratePhotons = totalnbPhotons;

  // Stopwatch<WallTime> s;
  //s.restart();

  vector< Photon > photonStack;
  photonStack.reserve(m_sendBufferSize);
  while( !done ){
    recvSize = photonStack.size();
    //generate and raytrace the inner photons
    

    CFuint nbCellPhotons =
        std::min(std::max(CFint(m_nbRaysCycle) - CFint(recvSize),(CFint)0), CFint(toGenerateCellPhotons) );

    CFuint nbWallPhotons =
        std::min(std::max(CFint(m_nbRaysCycle) - CFint(recvSize) - CFint(nbCellPhotons),(CFint)0), CFint(toGenerateWallPhotons ));

    for(CFuint i=0; i < nbCellPhotons ; ++i ){
       
      if(getCellPhotonData( photon )){
        //CFLog(INFO,"PHOTON: " << photon.cellID<<' '<<photon.userData.KS<<'\n' );
        //printPhoton(photon);
        rayTracing(photon);
      }
      --toGenerateCellPhotons;
      if (m_myProcessRank == 0)  ++*(progressBar);
    }

    for(CFuint i=0; i < nbWallPhotons ; ++i ){
      if(getFacePhotonData( photon )){
	//printPhoton(photon);
        rayTracing(photon);
      }
      -- toGenerateWallPhotons;
      if (m_myProcessRank == 0)  ++*(progressBar);
    }

//    CFLog(INFO, "raytrace the outer photons \n");
    for(CFuint i = 0; i< photonStack.size(); ++i ){
      //photon=photonStack[i];
      //CFLog(INFO,"PHOTON: " << photon.cellID<<' '<<photon.userData.KS<<'\n' );
      //printPhoton(photonStack[i]);
      rayTracing( photonStack[i] );
    }

    //sincronize
    //      CFLog(INFO, "sincronizing\n");
    bool isLastPhoton = (toGenerateCellPhotons + toGenerateWallPhotons == 0);
    done = m_lagrangianSolver.sincronizeParticles(photonStack, isLastPhoton);

    CFLog(VERBOSE,"Received "<<photonStack.size()<< " photons and has generated "<< nbCellPhotons <<" photons\n");
    CFLog(VERBOSE,"Number of photons left: "<< toGenerateCellPhotons <<"\n");
  }
  delete progressBar;

  //CFLog(INFO,"Raytracing took "<<s.readTimeHMS().str()<<'\n');
}

/////////////////////////////////////////////////////////////////////////////

template<class PARTICLE_TRACKING>
void RadiativeTransferMonteCarlo<PARTICLE_TRACKING>::printPhoton(Photon photon)
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
	);
}

/////////////////////////////////////////////////////////////////////////////

template<class PARTICLE_TRACKING>
void RadiativeTransferMonteCarlo<PARTICLE_TRACKING>::MonteCarlo()
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
CFuint RadiativeTransferMonteCarlo<PARTICLE_TRACKING>::rayTracing(Photon& beam)
{
  using namespace std;
  using namespace COOLFluiD::Framework;
  using namespace COOLFluiD::MathTools;
  using namespace COOLFluiD::Common;
  using namespace COOLFluiD::Numerics::FiniteVolume;
  using namespace COOLFluiD::LagrangianSolver;
  
  CFuint nbCrossedCells = 0;
  CFuint nbIter = 0;
  CFint exitCellID = 0;
  CFint exitFaceID = 0; 
  CFint currentCellID = 0;
  
  CFLog(DEBUG_MED, "RadiativeTransferMonteCarlo::rayTracing() => START\n");
  
  m_lagrangianSolver.newParticle(beam);
  
  CFLog(DEBUG_MED, "RadiativeTransferMonteCarlo::rayTracing() => particle ID: "<<beam.commonData.cellID<< "\n");
  
  PhotonData &beamData = m_lagrangianSolver.getUserDataPtr();
  exitCellID=m_lagrangianSolver.getExitCellID();
  
  //bool foundEntity = false;
  //cout<<"Start K= "<<previousK<<endl;
  //cout<<"Beam KS= "<<beam.KS<<endl;
  while(nbIter <= m_maxVisitedCells){
    //CFLog(INFO, "New step!\n");
    nbIter++;
    //GeoEntOut out;
    //if(_Axi) {
    //cout<<"beam.tt before: "<<beam.tt<<endl;
    //cout<<"beam.tt after: "<<beam.tt<<endl;
    
    currentCellID = exitCellID;

    m_lagrangianSolver.trackingStep();
    exitFaceID=m_lagrangianSolver.getExitFaceID();
    exitCellID=m_lagrangianSolver.getExitCellID();

    if(exitFaceID>=0){

      const CFreal stepDistance=m_lagrangianSolver.getStepDistance();
      RealVector null;
      
      //CFLog(INFO, "Absorption IN!\n");
      const CFreal cellK= m_radiation->getCellDistPtr(currentCellID)
          ->getRadiatorPtr()->getAbsorption(beamData.wavelength, null);

      //CFLog(INFO, "Absorption OUT!\n");

      beamData.KS -= stepDistance*cellK;
      //cout<<"beam data: "<<beamData.KS<<endl;
      //cout<<"currentCellID: "<< currentCellID <<endl;

      //_actualPoint = _internalPoint;
      //_intersectionPointOld = _intersectionPoint;
      //cout<<"New K= "<<beamData.KS<<" cellK: "<<cellK<<endl;
      if(beamData.KS <= 0.){ // photon absorbed by a cell
        //entity = INTERNAL_CELL;
        //CFLog(INFO, "PARTICLE ABSORBED!!\n");
        const CFuint gEndId = currentCellID;
        //CFuint gStartId = beam.emittingEntityID;
        const CFreal energyFraction = beamData.energyFraction;
        //        cout<<"sizeBuffer: "<< m_gInRadPowers.size()<<endl;
        //add directly to the in Rad Heat Power vector
        //cout<<"energy added : "<<energyFraction<<endl;
        cf_assert(gEndId < m_stateInRadPowers.size());
        m_stateInRadPowers[gEndId]+=energyFraction;
        //cout<<"new energy: "<<m_gInRadPowers[gEndId]<<endl;
        //foundEntity = true;
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
          //cout<<"ABSORBED!"<<endl;
          //entity = WALL_FACE;
          const CFuint ghostStateID = m_lagrangianSolver.getWallGhotsStateId(exitFaceID);
	  m_ghostStateInRadPowers[ghostStateID] += beamData.energyFraction;
	  CFLog(DEBUG_MIN, "Rad power in ghostStateID[" << ghostStateID << "] = " << 
		m_ghostStateInRadPowers[ghostStateID] << "\n");
	  //foundEntity = true;
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
      //CFLog(INFO,"HERE DOMAIN FACE!!\n");
          m_lagrangianSolver.bufferCommitParticle(exitFaceID);
        return 0;
      }

      nbCrossedCells++;
    }
  else{
    CFLog(VERBOSE, "RadiativeTransferMonteCarlo::rayTracing() => enter negligible\n");
    //entity = NEGLIGIBLE;
    return 0;
  }
  }
  CFLog(INFO, "RadiativeTransferMonteCarlo::rayTracing() => Max number of steps reached! \n");
  return 0;
}

/////////////////////////////////////////////////////////////////////////////

template<class PARTICLE_TRACKING>
void RadiativeTransferMonteCarlo<PARTICLE_TRACKING>::computeHeatFlux()
{
  using namespace std;
  using namespace COOLFluiD::Framework;
  using namespace COOLFluiD::MathTools;
  using namespace COOLFluiD::Common;
  using namespace COOLFluiD::Numerics::FiniteVolume;
  using namespace COOLFluiD::LagrangianSolver;
  
  CFLog(VERBOSE, "RadiativeTransferMonteCarlo computeHeatFlux() START\n");
  
  //  cout<<endl<<"COMPUTE HEAT FLUX"<<endl;
  //  cout<<"delta wave: "<<m_deltaWavs[m_spectralIdx]<<endl;

  DataHandle<State*, GLOBAL> states = socket_states.getDataHandle();
  DataHandle<CFreal> volumes = socket_volumes.getDataHandle();

  DataHandle<CFreal> gasRadiativeHeatSource = socket_qrad.getDataHandle();
  
  m_stateInRadPowers.sincronizeAssign();
  m_stateRadPower.sincronizeAssign();
  //fill qrad socket (loop states)
  for(CFuint istate = 0; istate < states.size(); ++istate){
    //if (states[istate]->isParUpdatable() ){ //is it necessary ?
      CFreal volume = volumes[istate];
      volume *= (m_isAxi ? 6.283185307179586*(*states[istate]).getCoordinates()[YY] : 1.);
      //gasRadiativeHeatSource[istate] = (- m_stateInRadPowers[istate] + m_stateRadPower[istate]) / volume ;
      gasRadiativeHeatSource[istate] = (1.-m_relaxationFactor)*gasRadiativeHeatSource[istate] + 
                                       m_relaxationFactor*(m_stateInRadPowers[istate]-m_stateRadPower[istate])/volume;
      //cout<<"state: "<<m_stateInRadPowers[istate]<<
      //	      ' '<< m_stateRadPower[istate]<<' '<<gasRadiativeHeatSource[istate]<<endl;
    
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
  //  CFLog(INFO, "############# nbWallTrs = " << nbWallTrs << "\n");
 
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
      CFuint faceGeoID = face->getID();
      CFuint faceGhostStateID = face->getState(1)->getLocalID();
      CFreal area = faceAreas[faceGeoID];
      area *= (m_isAxi) ? 6.283185307179586*faceCenters[faceGeoID*DIM_2D + YY] : 1.;
   
      //  CFLog(INFO, "face[" << f << "] =>" <<  m_ghostStateInRadPowers[faceGhostStateID] << ", " 
      //	<< m_ghostStateRadPower[faceGhostStateID] << "\n");
      wallRadiativeHeatSource[ff] = (-m_ghostStateInRadPowers[faceGhostStateID]
                                     +m_ghostStateRadPower[faceGhostStateID] ) / area;
      m_wallFaceBuilder.releaseGE();
    }
    //cout<<"CPU "<< Common::PE::GetPE().GetRank()<<" done " <<endl;
  }
  m_postProcess->runPostProcess(gasRadiativeHeatSource); 
  
  CFLog(VERBOSE, "RadiativeTransferMonteCarlo computeHeatFlux() END\n");
}

/////////////////////////////////////////////////////////////////////////////

  } // namespace RadiativeTransfer

} // namespace COOLFluiD

/////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_RadiativeTransfer_RadiativeTransferMonteCarloFVMCC_hh
