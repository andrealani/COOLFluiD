#ifndef COOLFluiD_RadiativeTransfer_RadiativeTransferMonteCarloFVMCC_ci
#define COOLFluiD_RadiativeTransfer_RadiativeTransferMonteCarloFVMCC_ci

//////////////////////////////////////////////////////////////////////////////
#include <utility>
#include <time.h>
#include <cmath>
#include <algorithm>
#ifdef CF_HAVE_BOOST_1_85
#define BOOST_TIMER_ENABLE_DEPRECATED
#endif
#include <boost/progress.hpp>
#include <boost/random.hpp>

#include "MathTools/MathChecks.hh"
#include "MathTools/MathConsts.hh"
#include "RadiativeTransfer/RadiativeTransferModule.hh"
#include "Framework/DataProcessingData.hh"
#include "Framework/ConvectiveVarSet.hh"
#include "Framework/DataSocketSource.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/GeometricEntity.hh"
#include "Framework/GeometricEntityPool.hh"
#include "Framework/FaceTrsGeoBuilder.hh"
#include "Framework/CellTrsGeoBuilder.hh"
#include "RadiativeTransfer/RadiationLibrary/RadiationLibrary.hh"
#include "Framework/PhysicalConsts.hh"
#include "Common/CFPrintContainer.hh"
#include "Common/MPI/MPIError.hh"
#include "FiniteVolume/CellCenterFVM.hh"
#include "MathTools/MathFunctions.hh"
#include "Common/MPI/MPIStructDef.hh"
#include "Common/SocketBundleSetter.hh"
#include "LagrangianSolver/ParallelVector/ParallelVector.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Numerics::FiniteVolume;
using namespace COOLFluiD::LagrangianSolver;


namespace COOLFluiD {

namespace RadiativeTransfer {


template<class PARTICLE_TRACKING>
void RadiativeTransferMonteCarlo<PARTICLE_TRACKING>::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFuint >("numberOfRays","number of rays sent by each elemet.");
  options.addConfigOption< CFuint >("MaxNbVisitedCells","Maximum number of visited cells.");
  options.addConfigOption< bool >("Axi","True if it is an axisymmetric simulation.");
  options.addConfigOption< string >("PostProcessName","Name of the post process routine");
  options.addConfigOption< CFuint >("sendBufferSize","Size of the buffer for comunication");
}

//////////////////////////////////////////////////////////////////////////////

template<class PARTICLE_TRACKING>
RadiativeTransferMonteCarlo<PARTICLE_TRACKING>::RadiativeTransferMonteCarlo(const std::string& name):
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
  m_radiation(new RadiationPhysicsHandler("RadiationPhysicsHandler"))
{
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
  //cout<<"Configure\n";

  CFAUTOTRACE;
  //cout<<"Configure\n";

  DataProcessingCom::configure(args);
  cf_assert(m_radiation.isNotNull());
  configureNested ( m_radiation.getPtr(), args );
  //cout<<"DONE!\n";

//  m_radLibrary = Environment::Factory<RadiationLibrary>::getInstance().
//      getProvider(m_radLibraryName)->create(m_radLibraryName);
//  cf_assert(m_radLibrary.isNotNull());
//  configureNested ( m_radLibrary.getPtr(), args );

  //cout<<"Configure Postprocess\n";
  m_postProcess = Environment::Factory< PostProcess >::getInstance().
     getProvider(m_postProcessName)->create(m_postProcessName);
  cf_assert(m_postProcess.isNotNull());
  configureNested ( m_postProcess.getPtr(), args );
  //cout<<"DONE!\n";
}

//////////////////////////////////////////////////////////////////////////////
template<class PARTICLE_TRACKING>
vector<SafePtr<BaseDataSocketSink> > RadiativeTransferMonteCarlo<PARTICLE_TRACKING>::needsSockets(){

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
  result.push_back(&socket_faceCenters);

  return result;
}

//////////////////////////////////////////////////////////////////////////////
template<class PARTICLE_TRACKING>
vector<SafePtr<BaseDataSocketSource> > RadiativeTransferMonteCarlo<PARTICLE_TRACKING>::providesSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSource> > result;
  result.push_back(&socket_qrad);
  result.push_back(&socket_qradFluxWall);

  return result;
}

/////////////////////////////////////////////////////////////////////////////
template<class PARTICLE_TRACKING>
void RadiativeTransferMonteCarlo<PARTICLE_TRACKING>::setup()
{

  // MPI parameters
  m_myProcessRank = PE::GetPE().GetRank();
  m_nbProcesses = PE::GetPE().GetProcessorCount();
  m_comm = PE::GetPE().GetCommunicator();
  MPIError::getInstance().init(m_comm, m_myProcessRank);

  // set dimensions
  m_dim = PhysicalModelStack::getActive()->getDim();
  m_dim2=(m_isAxi? 3 : m_dim);

  Common::SocketBundle sockets;
  sockets.states      = socket_states;
  sockets.gstates     = socket_gstates;
  sockets.nodes       = socket_nodes;
  sockets.nstates     = socket_nstates;
  sockets.normals     = socket_normals;
  sockets.isOutward   = socket_isOutward;
  sockets.volumes     = socket_volumes;
  sockets.faceCenters = socket_faceCenters;
  sockets.faceAreas   = socket_faceAreas;

  //m_radLibrary.setDataSockets(sockets);

  // initialize ParticleTracking
  m_lagrangianSolver.setDataSockets(sockets);
  m_lagrangianSolver.setupSendBufferSize(m_sendBufferSize);

  //initialize PostProcessign
  m_postProcess->setDataSockets(sockets);

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

  m_radiation->setupDataSockets(sockets);
  m_radiation->configureTRS();
  m_radiation->setupAxiFlag(m_isAxi);

  //MPI_Datatype Userdatatype, particleDatatype;
  MPIStruct Userdatatype, particleDatatype;

  PhotonData photonData;
  int counts[3] = {1,1,1};
  MPIStructDef::buildMPIStruct<CFreal,CFreal,CFreal>
          (&photonData.KS, &photonData.energyFraction, &photonData.wavelength, counts , Userdatatype);

  m_lagrangianSolver.setupParticleDatatype( Userdatatype.type );
  particleDatatype.type = m_lagrangianSolver.getParticleDataType();

  // resize the storage of qrad in each cell
  socket_qrad.getDataHandle().resize(nCells);
  m_gInRadPowers.resize(nCells,0);
  m_RHSGas.resize(nCells, 0.);

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
  m_RHSWall.resize(nbFaces);
}
/////////////////////////////////////////////////////////////////////////////
template<class PARTICLE_TRACKING>
void RadiativeTransferMonteCarlo<PARTICLE_TRACKING>::executeOnTrs()
{
  CFLog(VERBOSE, "RadiativeTransferMonteCarlo::executeOnTrs() START\n");

  // for testcases == 0, the following must be recomputed everytime
  Stopwatch<WallTime> s;
  s.restart();

  // re-initialize the heat fluxes to 0

  m_gInRadPowers=0.;
  m_RHSGas.assign(m_RHSGas.size(), 0.);

  MonteCarlo();

  computeHeatFlux();

  getRHF();

  //  }

  CFLog(DEBUG_MAX, "MonteCarlo() took " << s << "s\n");


//  //cout<<"HEREASASDDS"<<std::flush<<std::endl;
//  //fill wallQrad socket
//  DataHandle<CFreal> faceAreas = socket_faceAreas.getDataHandle();
//  DataHandle<CFreal> wallRadiativeHeatSource = socket_qradFluxWall.getDataHandle();
//  FaceTrsGeoBuilder::GeoData& wallFacesData = m_wallFaceBuilder.getDataGE();

//  const CFuint nbWallTrs = wallTrsNames.size();
//  for(CFuint i=0; i<nbWallTrs; ++i){
//    SafePtr<TopologicalRegionSet> WallFaces = MeshDataStack::getActive()->getTrs(wallTrsNames[i]);
//    const CFuint nbFacesWall = WallFaces->getLocalNbGeoEnts();
//    wallFacesData.trs = WallFaces;
//    for(CFuint f=0; f<nbFacesWall; ++f){
//      wallFacesData.idx = f;
//      Framework::GeometricEntity *const face = m_wallFaceBuilder.buildGE();
//      CFuint id = face->getID(); //LIKE THIS ??
//      CFuint sID = m_lagrangianSolver.getFaceStateID(id);
//      if( m_lagrangianSolver.getFaceOwnerRank(id) == m_myProcessRank ){
//        //cout<<"id: "<<id<<flush<<endl;
//        RealVector centroid = face->computeCentroid();
//        CFreal area = faceAreas[id] * (m_isAxi ? 6.283185307179586*centroid[YY] : 1.);
//        wallRadiativeHeatSource[f] = m_RHSGas[sID]/area;
//        // CFLog(INFO, " AFTER faceID = " << WallFaces->getLocalGeoID(f) << ", index = " << ff << "\n");
//      }
//      m_wallFaceBuilder.releaseGE();
//    }
//  }
//  //cout<<"HEREASASDDS"<<std::flush<<std::endl;
//  CFLog(VERBOSE, "RadiativeTransferMonteCarlo::executeOnTrs() END\n");
}

/////////////////////////////////////////////////////////////////////////////

template<class PARTICLE_TRACKING>
void RadiativeTransferMonteCarlo<PARTICLE_TRACKING>::getTotalEnergy(){

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
                               ->getRadiatorPtr()->getSpectaLoopPower();
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
                                                  ->getRadiatorPtr()->getSpectaLoopPower();
        //cout<<"THE OTHER: m_axi Volume: "<<m_axiVolumes[i]<<endl;
        //m_stateRadPower[i] = cellK*sigma*pow(T,4.)*m_axiVolumes[i];
        gStateRadPower += m_ghostStateRadPower[gstate];

      }
      else{
       //cout<<"is ghost!"<<endl;
       m_ghostStateRadPower[gstate] = .0;
      }
    }

    m_totalRadPower = stateRadPower+ gStateRadPower;

    //cout<<"statePhotons"<<endl;
    for(CFuint i= 0; i< m_stateRadPower.size(); ++i){
        //CFuint test = CFuint(m_stateRadPower[i]/m_totalRadPower*nbPhotons);
        m_nbPhotonsState[i] =( m_stateRadPower[i] > 0. )? m_nbRaysElem : 0 ;

        cout<<"nbPhotons: "<<m_nbPhotonsState[i]<< ' '<<m_stateRadPower[i] <<endl;
    }

    cout<<"ghostPhotons"<<endl;
    for(CFuint i= 0; i< m_ghostStateRadPower.size(); ++i){
        //CFuint test = CFuint(m_stateRadPower[i]/m_totalRadPower*nbPhotons);
        m_nbPhotonsGhostState[i] =( m_ghostStateRadPower[i] > 0.)? m_nbRaysElem : 0 ;
        cout<<"nbPhotons: "<<m_nbPhotonsGhostState[i]<< ' '<<m_ghostStateRadPower[i] <<endl;
    }

}
/////////////////////////////////////////////////////////////////////////////
template<class PARTICLE_TRACKING>
bool RadiativeTransferMonteCarlo<PARTICLE_TRACKING>::getCellPhotonData(Photon &ray){

    //static CellTrsGeoBuilder::GeoData& cellData = m_cellBuilder.getDataGE();
    static CFuint photon=0, state=0;
    static CFuint nbStates = m_nbPhotonsState.size();

    for(;state<nbStates; ++state)
    {
      for(;photon<m_nbPhotonsState[state];)
      {
//        cout<<"m_nbPhotonsState[state]: "<<m_nbPhotonsState[state]<<' '<<";state: "<<state<<endl;
//        cout<<"photon: "<<photon<<endl;

        //cellData.idx =  state;

        //cellData.idx = 50;
        //GeometricEntity *const cell = m_cellBuilder.buildGE();
        //if( cell->getState(0)->isParUpdatable() ){
        // Calculate the wavelength
        //cout<<"wavelength= "<<ray.userData.wavelength<<endl;

        //Get directions
        RealVector directions(m_dim2);
        m_radiation->getCellDistPtr(state)->
            getRadiatorPtr()->getRandomEmission(ray.userData.wavelength, directions );

        //cout<<"ray directions: ";
        for(CFuint ii=0; ii<m_dim2; ++ii){
          ray.commonData.direction[ii]=directions[ii];
        }
        //cout<<endl;

        //Get the beam max optical path Ks
        ray.userData.KS = - std::log( m_rand.uniformRand() );
        // ray.actualKS = 0;

        //Get cell center
        static Framework::DataHandle<Framework::State*, Framework::GLOBAL> states
                = socket_states.getDataHandle();

        Node& baricenter = (*states[state]).getCoordinates();

        //RealVector baricenter(_dim);
        //m_particleTracking->computeAverage(*cell->getNodes(), cell->nbNodes(), baricenter);

        //cout<<"baricenter: ";

        //cout<<endl;

        for(CFuint i=0;i<m_dim;++i){
          //  cout<<baricenter[i]<<' ';
          ray.commonData.currentPoint[i]=baricenter[i];
        }
        for(CFuint i=m_dim;i<3;++i){
          //  cout<<baricenter[i]<<' ';
          ray.commonData.currentPoint[i]=0.;
        }

        //Get the remaining information
        //const CFuint cellID = cell->getID();
        CFuint cellID = m_radiation->getCurrentCellStateID();
        ray.commonData.cellID = cellID;

        ray.userData.energyFraction= m_stateRadPower[state]/CFreal(m_nbPhotonsState[state]);

        //m_cellBuilder.releaseGE();
        ++photon;
        return true;
      }
       photon=0;
    }
    return false;
}

template<class PARTICLE_TRACKING>
bool RadiativeTransferMonteCarlo<PARTICLE_TRACKING>::getFacePhotonData(Photon &ray){

    //static CellTrsGeoBuilder::GeoData& cellData = m_cellBuilder.getDataGE();
    static CFuint photon=0, gState=0;
    static CFuint nbGstates = m_nbPhotonsGhostState.size();

    for(; gState < nbGstates ; ++gState)
    {
      for(;photon<m_nbPhotonsGhostState[gState];)
      {
//        cout<<"m_nbPhotonsState[state]: "<<m_nbPhotonsState[state]<<' '<<";state: "<<state<<endl;
//        cout<<"photon: "<<photon<<endl;

        //cellData.idx =  state;

        //cellData.idx = 50;
        //GeometricEntity *const cell = m_cellBuilder.buildGE();
        //if( cell->getState(0)->isParUpdatable() ){
        // Calculate the wavelength
        //cout<<"wavelength= "<<ray.userData.wavelength<<endl;

        //Get directions

        RealVector directions(m_dim2);
        m_radiation->getWallDistPtr( gState )->
            getRadiatorPtr()->getRandomEmission(ray.userData.wavelength, directions );

        CFuint faceGeoID = m_radiation->getCurrentWallGeoID();

        //cout<<"ray directions: ";
        for(CFuint ii=0; ii<m_dim2; ++ii){
          ray.commonData.direction[ii]=directions[ii];
        }
        //cout<<endl;

        //Get the beam max optical path Ks
        ray.userData.KS = - std::log( m_rand.uniformRand() );
        // ray.actualKS = 0;

        //Get the face center
        DataHandle<CFreal> faceCenters = socket_faceCenters.getDataHandle();

        //cout<<"baricenter: ";
        for(CFuint i=0;i<m_dim;++i){
            //cout<<faceCenters[m_dim * faceGeoID + i]<<' ';
          ray.commonData.currentPoint[i]=faceCenters[m_dim * faceGeoID + i];
        }
        for(CFuint i=m_dim;i<3;++i){
          //cout<<"0 ";
          ray.commonData.currentPoint[i]=0.;
        }
        //cout<<endl;

        //Get the remaining information
        //const CFuint cellID = cell->getID();
        CFuint cellID = m_lagrangianSolver.getWallStateId( faceGeoID );
        ray.commonData.cellID = cellID;

        ray.userData.energyFraction= m_ghostStateRadPower[gState]/CFreal(m_nbPhotonsGhostState[gState]);

        //m_cellBuilder.releaseGE();
        ++photon;
        return true;
      }
       photon=0;
    }
    return false;
}

/////////////////////////////////////////////////////////////////////////////
template<class PARTICLE_TRACKING>
void RadiativeTransferMonteCarlo<PARTICLE_TRACKING>::computePhotons()
{
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


  CFuint recvSize = 0;
  bool done = false;

  Photon photon;
  CFuint totalnbPhotons = toGenerateWallPhotons + toGenerateCellPhotons;
  boost::progress_display* progressBar = NULL;
  if (m_myProcessRank == 0) progressBar = new boost::progress_display(totalnbPhotons);

  //CFuint toGeneratePhotons = totalnbPhotons;

  Stopwatch<WallTime> s;
  s.restart();

  vector< Photon > photonStack;
  photonStack.reserve(m_sendBufferSize);
  while( !done ){
      recvSize = photonStack.size();
      //generate and raytrace the inner photons
      //CFLog(INFO, recvSize<<' '<<toGeneratePhotons<<'\n');
      CFuint nbCellPhotons = std::min(m_sendBufferSize - recvSize, toGenerateCellPhotons );
      CFuint nbWallPhotons = std::min(m_sendBufferSize - recvSize - nbCellPhotons, toGenerateWallPhotons );

      for(CFuint i=0; i < nbCellPhotons ; ++i ){
        if(getCellPhotonData( photon )){
          //CFLog(INFO,"PHOTON: " << photon.cellID<<' '<<photon.userData.KS<<'\n' );
          //printPhoton(photon);
          rayTracing(photon);
        }
        -- toGenerateCellPhotons;
          if (m_myProcessRank == 0)  ++*(progressBar);
      }

      for(CFuint i=0; i < nbWallPhotons ; ++i ){
        if(getFacePhotonData( photon )){
          //CFLog(INFO,"PHOTON: " << photon.cellID<<' '<<photon.userData.KS<<'\n' );
          //printPhoton(photon);
          rayTracing(photon);
        }
        -- toGenerateWallPhotons;
          if (m_myProcessRank == 0)  ++*(progressBar);
      }

      //CFLog(INFO, "raytrace the outer photons \n");
      for(CFuint i = 0; i< photonStack.size(); ++i ){
          //photon=photonStack[i];
          //CFLog(INFO,"PHOTON: " << photon.cellID<<' '<<photon.userData.KS<<'\n' );
          //printPhoton(photonStack[i]);
          rayTracing( photonStack[i] );
      }

      //sincronize
      bool isLastPhoton = (toGenerateCellPhotons + toGenerateWallPhotons == 0);
      done = m_lagrangianSolver.sincronizeParticles(photonStack, isLastPhoton);

      //CFLog(INFO, "done Sincronizing; "<< "receiv_size= "<< photonStack.size() << " photons to Generate= "<<toGeneratePhotons <<"\n");

//      if(m_myProcessRank == 0){
//       cout<<"Done sincronizing"<<endl<<"CPU "<<m_myProcessRank<< " has received: "<<photonStack.size()<<
//             " photons and has generated "<< totalnbPhotons-toGeneratePhotons <<" photons"<<endl;
//      }
  }
   delete progressBar;

  CFLog(INFO,"Raytracing took "<<s.readTimeHMS().str()<<'\n');
}

/////////////////////////////////////////////////////////////////////////////

template<class PARTICLE_TRACKING>
void RadiativeTransferMonteCarlo<PARTICLE_TRACKING>::printPhoton(Photon photon){
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
    CFuint nbLoops = m_radiation->getNumberLoops();

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
  CFuint nbCrossedCells = 0;
  CFuint nbIter = 0;
  static CFint exitCellID, exitFaceID, currentCellID;
  //CFLog(INFO, "start Raytracing\n");

  //CFLog(DEBUG_MAX, "RadiativeTransferMonteCarlo::rayTracing() => actualCellID = " << actualCellID << "\n");

  //CFLog(INFO,"NEW particle!\n");
  m_lagrangianSolver.newParticle(beam);
  //cout<<"particle ID: "<<beam.commonData.cellID<<endl;
  PhotonData &beamData = m_lagrangianSolver.getUserDataPtr();
  exitCellID=m_lagrangianSolver.getExitCellID();

  bool foundEntity = false;
  //cout<<"Start K= "<<previousK<<endl;
  //cout<<"Beam KS= "<<beam.KS<<endl;
  while(nbIter <= m_maxVisitedCells){
    //CFLog(INFO, "New step!\n");
    nbIter++;
    //GeoEntOut out;
    //if(_Axi){
    //cout<<"beam.tt before: "<<beam.tt<<endl;
    //out=m_particleTracking.myAxiRayTracing(beam,actualCellID);
    //cout<<"beam.tt after: "<<beam.tt<<endl;

    currentCellID = exitCellID;

    m_lagrangianSolver.trackingStep();
    exitFaceID=m_lagrangianSolver.getExitFaceID();
    exitCellID=m_lagrangianSolver.getExitCellID();

    if(exitFaceID>=0){
      static CFreal cellK, stepDistance;

      stepDistance=m_lagrangianSolver.getStepDistance();
      RealVector null;

      //CFLog(INFO, "Absorption IN!\n");
      cellK= m_radiation->getCellDistPtr(currentCellID)
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
        //CFLog(INFO, "PARTICLE ABSORVED!!\n");
        CFuint gEndId = currentCellID;
        //CFuint gStartId = beam.emittingEntityID;
        CFreal energyFraction = beamData.energyFraction;
        //        cout<<"sizeBuffer: "<< m_gInRadPowers.size()<<endl;
        //add directly to the in Rad Heat Power vector
        //cout<<"energy added : "<<energyFraction<<endl;
        cf_assert(gEndId < m_gInRadPowers.size());
        m_gInRadPowers[gEndId]+=energyFraction;
        //cout<<"new energy: "<<m_gInRadPowers[gEndId]<<endl;
        foundEntity = true;
        return currentCellID;
      }

      CFuint faceType = m_lagrangianSolver.getFaceType(exitFaceID);

      if ( faceType == ParticleTracking::WALL_FACE){
        //CFLog(INFO,"HERE WALL !!\n");

        RealVector entryDirection(m_dim2), position(m_dim2);
        for(CFuint i=0; i < m_dim2; ++i){
          entryDirection[i]= beam.commonData.direction[i];
        }
        m_lagrangianSolver.getExitPoint(position);

        RealVector normal(m_dim2);
        //CFuint stateID = m_lagrangianSolver.getWallGhotsStateId(exitFaceID);
        CFuint ghostStateID = m_lagrangianSolver.getWallGhotsStateId(exitFaceID);

        m_lagrangianSolver.getNormals(exitFaceID,position,normal);

        const CFreal wallK = m_radiation->getWallDistPtr(ghostStateID)
            ->getRadiatorPtr()->getAbsorption( beamData.wavelength,entryDirection );

        const CFreal reflectionProbability =  m_rand.uniformRand();

        if(reflectionProbability <= wallK){ // the photon is absorved by the wall
          //cout<<"ABSORVED!"<<endl;
          //entity = WALL_FACE;
          foundEntity = true;
          return exitFaceID;
        }
        else {
          //m_radiation->getDistPtr(stateID)->reflectionDist->getRandomProperties(
          //      beamData.wavelength, exitDirection, normal, entryDirection);

          //m_lagrangianSolver.newDirection(exitDirection);
          return exitFaceID;

        }
      }

      if (faceType == ParticleTracking::BOUNDARY_FACE ){
        //CFLog(INFO,"HERE BOUNDARY !!\n");
        //entity = DISAPPEARED;
        foundEntity = true;
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
    CFLog(INFO, "enter negligle \n");
    //entity = NEGLIGIBLE;
    return 0;
  }
  }
  return 0;
}


template<class PARTICLE_TRACKING>
void RadiativeTransferMonteCarlo<PARTICLE_TRACKING>::computeHeatFlux()
{
  CFLog(VERBOSE, "RadiativeTransferMonteCarlo computeHeatFlux()\n");

  //  cout<<endl<<"COMPUTE HEAT FLUX"<<endl;
  //  cout<<"delta wave: "<<m_deltaWavs[m_spectralIdx]<<endl;

  DataHandle<State*, GLOBAL> states = socket_states.getDataHandle();
  DataHandle<CFreal> volumes = socket_volumes.getDataHandle();

  DataHandle<CFreal> gasRadiativeHeatSource = socket_qrad.getDataHandle();


  //fill qrad socket (loop states)
  for(CFuint istate = 0; istate < states.size(); ++istate){
    if (states[istate]->isParUpdatable() ){ //is it necessary ?
      CFreal volume = volumes[id];
      volume *= (m_isAxi ? 6.283185307179586*(*states[id]).getCoordinates()[YY] : 1.);
      gasRadiativeHeatSource[istate] = (- m_gInRadPowers[c] + m_stateRadPower[c]) / volume ;
    }
  }

  //Fill qradFluxWall socket (loop Trs idxs for each trs)
  //Needs to be consistent with AeroForces definition
  DataHandle<CFreal> faceAreas   = socket_faceAreas.getDataHandle();
  DataHandle<CFreal> faceCenters = socket_faceCenters.getDataHandle();
  DataHandle<CFreal> wallRadiativeHeatSource = socket_qradFluxWall.getDataHandle();

  vector<string> wallTrsNames;
  m_radiation->getWallTRSnames(wallTrsNames);

  const CFuint nbWallTrs = wallTrsNames.size();
  CFuint ff = 0;
  for(CFuint i=0; i<nbWallTrs; ++i){
    SafePtr<TopologicalRegionSet> WallFaces = MeshDataStack::getActive()->getTrs(_wallTrsNames[i]);
    const CFuint nbFacesWall = WallFaces->getLocalNbGeoEnts();
    for(CFuint f=0; f<nbFacesWall; ++f, ++ff){
      //cf_assert(m_wallSurfaces[ff] > 0.);
      wallRadiativeHeatSource[ff] = _RHSWall[ff]/m_wallSurfaces[ff];
      // CFLog(INFO, " AFTER faceID = " << WallFaces->getLocalGeoID(f) << ", index = " << ff << "\n");
    }
  }
}


template<class PARTICLE_TRACKING>
void RadiativeTransferMonteCarlo<PARTICLE_TRACKING>::getRHF(){

    CFreal T = 10000.;
    CFreal cellK = 1.;

    CFreal sigma = PhysicalConsts::StephanBolzmann();
    CFreal NonD = 4. *cellK*sigma*pow(T,4.);
    //CFreal opticalTreshold = -1./cellK *std::log(1. - fraction);

    Framework::DataHandle<Framework::State*, Framework::GLOBAL> states
            = socket_states.getDataHandle();

    DataHandle<CFreal> volumes = socket_volumes.getDataHandle();

    vector<CFreal> Q_out;
    Q_out.reserve(volumes.size() );

    cf_assert(volumes.size() == m_stateRadPower.size());
    cf_assert(volumes.size() == m_RHSGas.size() );

    CFreal volume = 0, value = 0;
    for(CFuint c=0; c<volumes.size(); ++c){
      volume = volumes[c] * (m_isAxi ? 6.283185307179586*(*states[c]).getCoordinates()[YY] : 1.);
      value = ( std::abs(m_RHSGas[c]) > 1e-6 ) ? (m_RHSGas[c]/ NonD /volumes[c]) : 0.;

      Q_out.push_back( value );
    }

    m_postProcess->runPostProcess(Q_out);
}



  } // namespace RadiativeTransfer

} // namespace COOLFluiD

#endif
