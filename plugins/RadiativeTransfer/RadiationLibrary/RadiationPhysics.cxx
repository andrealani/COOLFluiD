#include "RadiativeTransfer/RadiationLibrary/RadiationPhysics.hh"
#include "Framework/TopologicalRegionSet.hh"
#include "Framework/MeshData.hh"
#include "MathTools/RealVector.hh"

/////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

namespace RadiativeTransfer {

/////////////////////////////////////////////////////////////////////////////

void RadiationPhysics::defineConfigOptions(Config::OptionList& options)
{
    options.addConfigOption< std::string >
        ("ApplyTRS","TRS to apply the radiative library to");
    options.addConfigOption< std::string >
        ("TypeTRS","Type of TRS: Wall, Medium or Boundary");
    options.addConfigOption< std::string >
        ("Radiator","Radiator Physics Distribution");
    options.addConfigOption< std::string >
        ("Reflector","Reflector Physics Distribution");
}

/////////////////////////////////////////////////////////////////////////////

void RadiationPhysics::getWallStateIDs(std::vector<CFuint>& statesID,
                                        std::vector<CFuint>& wallGeoIdx)
{
  statesID.clear();
  wallGeoIdx.clear();
  Framework::SocketBundleSetter socketBundle;
  socketBundle.setDataSockets( *(m_radPhysicsHandlerPtr->getDataSockets()) );

  Framework::GeometricEntityPool<Framework::FaceTrsGeoBuilder>* faceTRSBuilder =
      socketBundle.getFaceTrsBuilder();

  Framework::FaceTrsGeoBuilder::GeoData& facesData = faceTRSBuilder->getDataGE();
  faceTRSBuilder->getDataGE().isBFace = true;

  Common::SafePtr<Framework::TopologicalRegionSet> WallFaces =
      Framework::MeshDataStack::getActive()->getTrs( m_TRSname );

  facesData.trs = WallFaces;
  const CFuint nbFaces = WallFaces->getLocalNbGeoEnts();
  statesID.reserve( nbFaces );
  for(CFuint i=0; i<nbFaces; ++i){
    facesData.idx = i;
    Framework::GeometricEntity *const face = faceTRSBuilder->buildGE();
     if( !face->getState(1)->isParUpdatable() && face->getState(0)->isParUpdatable() ){
      statesID.push_back( face->getState(1)->getLocalID() );
      wallGeoIdx.push_back( face->getID() );
    }
    faceTRSBuilder->releaseGE();
  }
}

/////////////////////////////////////////////////////////////////////////////

void RadiationPhysics::getCellStateIDs(std::vector<CFuint>& statesID)
{
  Framework::SocketBundleSetter socketBundle;
  socketBundle.setDataSockets( *(m_radPhysicsHandlerPtr->getDataSockets()) );

  Framework::GeometricEntityPool<Framework::CellTrsGeoBuilder>* cellBuilder =
      socketBundle.getCellTRSbuilder();

  Framework::CellTrsGeoBuilder::GeoData& cellsData = cellBuilder->getDataGE();

  Common::SafePtr<Framework::TopologicalRegionSet> cellTRS =
      Framework::MeshDataStack::getActive()->getTrs( m_TRSname );
  cellsData.trs = cellTRS;
  const CFuint nbCells = cellTRS->getLocalNbGeoEnts();
  statesID.reserve( nbCells );
  for(CFuint i=0; i<nbCells; ++i){
    cellsData.idx = i;
    Framework::GeometricEntity *const cell = cellBuilder->buildGE();
    cf_assert(!cell->getState(0)->isGhost());
    statesID.push_back( cell->getState(0)->getLocalID() );
    cellBuilder->releaseGE();
  }
}

//////////////////////////////////////////////////////////////////////////////

void RadiationPhysics::computeInterpolatedStates()
{
  if(m_TRStypeID == WALL){
    //std::cout<<"isWall!"<<std::endl;

    const CFuint dim = Framework::PhysicalModelStack::getActive()->getDim();

    Framework::SocketBundleSetter socketBundle;
    socketBundle.setDataSockets( *(m_radPhysicsHandlerPtr->getDataSockets()) )  ;

    Framework::DataHandle< RealVector> nstates = socketBundle.getDataSocket()->nstates.getDataHandle();
    Framework::DataHandle< CFreal > faceCenters = socketBundle.getDataSocket()->faceCenters.getDataHandle();

    //std::cout<<"is WALL OR BOUNDARY"<<std::endl;
    Framework::GeometricEntityPool<Framework::FaceTrsGeoBuilder>* faceTRSBuilder =
        socketBundle.getFaceTrsBuilder();

    Framework::FaceTrsGeoBuilder::GeoData& facesData = faceTRSBuilder->getDataGE();
    faceTRSBuilder->getDataGE().isBFace = true;

    Common::SafePtr<Framework::TopologicalRegionSet> WallFaces =
        Framework::MeshDataStack::getActive()->getTrs( m_TRSname );

    facesData.trs = WallFaces;
    const CFuint nbFaces = WallFaces->getLocalNbGeoEnts();
    m_interpolatedStates.resize(nbFaces);
    RealVector nodeCenter(dim);
    
    for(CFuint i=0; i<nbFaces; ++i){
      facesData.idx = i;
      Framework::GeometricEntity *const face = faceTRSBuilder->buildGE();
      const std::vector<Framework::Node*>& nodesInFace = *face->getNodes();
      const CFuint nbNodesInFace = nodesInFace.size();
      m_interpolatedStates[i] = 0.0;
      nodeCenter = 0.;
      for (CFuint node = 0; node < nbNodesInFace; ++node) {
	const Framework::State& nstate = nstates[nodesInFace[node]->getLocalID()];
        m_interpolatedStates[i]  += nstate;
        nodeCenter += *nodesInFace[node]->getData();
      }
      m_interpolatedStates[i] /= nbNodesInFace;
      nodeCenter/= nbNodesInFace;
      for(CFuint d=0;d<dim;++d){
        faceCenters[face->getID()*dim+d] = nodeCenter[d];
      }
      faceTRSBuilder->releaseGE();
    }
  }
}


//////////////////////////////////////////////////////////////////////////////

RadiationPhysics::RadiationPhysics(const std::string& name) :
  Common::OwnedObject(),
  ConfigObject(name)
{
  addConfigOptionsTo(this);

  m_TRSname = "Null";
  m_TRStype = "Null";

  setParameter("ApplyTRS"  , &m_TRSname  );
  setParameter("TypeTRS"  ,  &m_TRStype  );

  m_radiatorName   = "NullRadiator";
  m_reflectorName  = "NullReflector";
  //m_reflectionName = "NullReflection";
  //m_scatteringName = "NullScattering";

  setParameter("Radiator"  , &m_radiatorName  );
  setParameter("Reflector" , &m_reflectorName );
  //setParameter("Scattering", &m_scatteringName);
}

//////////////////////////////////////////////////////////////////////////////

void RadiationPhysics::configure(Config::ConfigArgs& args)
{
  Config::ConfigObject::configure(args);
  
  Common::SafePtr<Radiator::PROVIDER> radiatorProv = CFNULL;
  try {
    radiatorProv = Environment::Factory<Radiator>::getInstance().
      getProvider(m_radiatorName);
  }
  catch (Common::NoSuchValueException& e) {
    CFLog(WARN, "RadiationPhysics::configure() => creating default NullRadiator\n");
    m_radiatorName = "NullRadiator";
    radiatorProv = Environment::Factory<Radiator>::getInstance().
      getProvider(m_radiatorName);
  }
  cf_assert(radiatorProv.isNotNull());
  m_radiator.reset(radiatorProv->create(m_radiatorName));
  cf_assert(m_radiator.isNotNull());
  
  m_reflector = Environment::Factory<Reflector>::getInstance().
      getProvider(m_reflectorName)->create(m_reflectorName);

  //scatteringDist = Environment::Factory<RadiationDistribution>::getInstance().
  //    getProvider(m_scatteringDistName)->create(m_scatteringDistName);

  cf_assert(   m_radiator.isNotNull() );
  cf_assert(  m_reflector.isNotNull() );

  //cf_assert( reflectionDist.isNotNull() );
  //cf_assert( scatteringDist.isNotNull() );

  configureNested ( m_radiator.getPtr(), args );
  configureNested ( m_reflector.getPtr(), args );


  //convert TRS type string to enum
  if (m_TRStype.compare("Wall") == 0){
    m_TRStypeID = WALL;
  } else{
    if (m_TRStype.compare("Medium") == 0){
      m_TRStypeID = MEDIUM;
    } else {
      if ( m_TRStype.compare("Boundary") ==0 ){
        m_TRStypeID = BOUNDARY;
      } else {
        CFLog(INFO, "ERROR : "<<m_TRSname<<" Undefined TRS type\n");
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void RadiationPhysics::setupSpectra(CFreal wavMin,CFreal wavMax)
{
  m_radiator->setupSpectra(wavMin, wavMax);
  m_reflector->setupSpectra(wavMin, wavMax);
  
  m_radiator->computeEmissionCPD();
  m_reflector->computeReflectionCPD();
  
  //scatteringDist ->setupSectra(wavMin, wavMax);
}
 
//////////////////////////////////////////////////////////////////////////////
 
void RadiationPhysics::setup()
{
  m_radiator->setRadPhysicsPtr(this);
  m_radiator->setRadPhysicsHandlerPtr(m_radPhysicsHandlerPtr);
  m_radiator->setup();

  m_reflector->setRadPhysicsPtr(this);
  m_reflector->setRadPhysicsHandlerPtr(m_radPhysicsHandlerPtr);
  m_reflector->setup();

  //absorptionDist ->setupSectra(wavMin, wavMax);
  //reflectionDist ->setupSectra(wavMin, wavMax);
  //scatteringDist ->setupSectra(wavMin, wavMax);
}
  
  //////////////////////////////////////////////////////////////////////////////
  
}
}
