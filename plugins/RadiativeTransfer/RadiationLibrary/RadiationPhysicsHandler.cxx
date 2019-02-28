#include "RadiationPhysicsHandler.hh"
#include "Framework/PhysicalChemicalLibrary.hh"
#include "Framework/MeshData.hh"
#include "MathTools/RealMatrix.hh"

////////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

namespace RadiativeTransfer {

////////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace RadiativeDistTypes;
  
////////////////////////////////////////////////////////////////////////////////

void RadiationPhysicsHandler::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal >
      ("WavelengthMin","Minimum wavelenght to consider [A] ");
  options.addConfigOption< CFreal >
      ("WavelengthMax","Maximum wavelenght to consider [A]");
  options.addConfigOption< CFuint >
      ("NumberLoops","Number of loops to divide the spectra calculations");
  options.addConfigOption< CFint >
      ("TempID","ID of the first temperature state");
  options.addConfigOption< bool >
      ("TwoTemperatureModel","Enables a two temperature model with Tr and Tv");
  options.addConfigOption< CFint >
      ("TempVID","ID of the vibrational temperature state");
  options.addConfigOption< std::vector <std::string> >
      ("RadiationPhysicsNames","Names to attach a Radiation Physics class");
}

//////////////////////////////////////////////////////////////////////////////

RadiationPhysicsHandler::RadiationPhysicsHandler(const std::string& name) :
  Common::OwnedObject(),
  ConfigObject(name),
  Common::NonCopyable<RadiationPhysicsHandler>(),
  m_sockets(),
  m_radiationPhysics(),
  m_statesOwner(),
  m_ghostStatesOwner(),
  m_wallTRSnames(),
  m_boundaryTRSnames(),
  m_mediumTRSnames(),
  m_nbTemps(1),
  m_cellStateID(0),
  m_cellStateOwnerIdx(0),
  m_ghostStateID(0),
  m_ghostStateOwnerIdx(0),
  m_ghostStateWallGeoID(0),
  m_isAxi(false)
{
  addConfigOptionsTo(this);
  
  m_wavMin = -1.0; // default value is on purpose <0
  setParameter("WavelengthMin", &m_wavMin);
  
  m_wavMax = -1.0; // default value is on purpose <0
  setParameter("WavelengthMax", &m_wavMax);
  
  m_nbLoops = 1;
  setParameter("NumberLoops",&m_nbLoops);
  
  m_TempID = -1;
  setParameter("TempID",&m_TempID);

  m_TempVID = -1;
  setParameter("TempVID",&m_TempVID);

  m_useTwoTemps = false;
  setParameter("TwoTemperatureModel",&m_useTwoTemps);

  setParameter("RadiationPhysicsNames",&m_radiationPhysicsNames);
}

//////////////////////////////////////////////////////////////////////////////

void RadiationPhysicsHandler::configure(Config::ConfigArgs& args)
{
  Config::ConfigObject::configure(args);
  m_radiationPhysics.resize( m_radiationPhysicsNames.size() );
  
  for(CFuint i=0 ; i< m_radiationPhysicsNames.size(); ++i){
    m_radiationPhysics[i] = new RadiationPhysics( m_radiationPhysicsNames[i] );
    cf_assert( m_radiationPhysics[i].isNotNull() );
    CFLog(VERBOSE, "RadiationPhysicsHandler::configure() => "<<m_radiationPhysicsNames[i]<<'\n');
    m_radiationPhysics[i]->setRadPhysicsHandlerPtr(this);
    configureNested ( m_radiationPhysics[i].getPtr(), args );
  }
}

//////////////////////////////////////////////////////////////////////////////

void RadiationPhysicsHandler::setup()
{
  //Assume temperature to be the last equation
  Common::SafePtr<Framework::PhysicalChemicalLibrary> library;
  library = Framework::PhysicalModelStack::getActive()->getImplementor()->
      getPhysicalPropertyLibrary< Framework::PhysicalChemicalLibrary >();

  m_nbTemps = 1;
  if (library.isNotNull()) {
    m_nbTemps = 1 + library->getNbTempVib() + library->getNbTe();
  }

  if(m_TempID==-1){
    const CFuint nbEqs = Framework::PhysicalModelStack::getActive()->getNbEq();
    // we assume here 2 (T-Tv) or 3 (T, Tv, Te) temperatures
    // translational temperature ID
    m_TempID = nbEqs - m_nbTemps - 1 ;
    //std::cout<<"TEMPERATURE ID: "<<m_TempID<<' '<<m_nbTemps<<std::endl;
  }

  for(CFuint i=0 ; i< m_radiationPhysics.size(); ++i){
    m_radiationPhysics[i]->setup();
  }
}

//////////////////////////////////////////////////////////////////////////////

RadiationPhysicsHandler::~RadiationPhysicsHandler(){

}

//////////////////////////////////////////////////////////////////////////////

void RadiationPhysicsHandler::setupWavStride(CFuint loop)
{
  cf_assert(m_nbLoops > 0);
  static CFreal wavStride = ( m_wavMax - m_wavMin ) / m_nbLoops;
  CFreal wavMin = m_wavMin + wavStride * loop;
  CFreal wavMax = wavMin + wavStride;
  
  for(CFuint i=0; i<m_radiationPhysics.size();++i) {
    CFLog(VERBOSE, "RadiationPhysicsHandler::setupWavStride() for [" << m_radiationPhysics[i]->getName() << "] START\n");
    m_radiationPhysics[i]->setupSpectra(wavMin, wavMax);
    CFLog(VERBOSE, "RadiationPhysicsHandler::setupWavStride() for [" << m_radiationPhysics[i]->getName() <<"] END\n");
    m_radiationPhysics[i]->computeInterpolatedStates();
    CFLog(VERBOSE, "RadiationPhysicsHandler::computeInterpolatedStates() for [" << m_radiationPhysics[i]->getName() <<" END\n");
  }
}
  
//////////////////////////////////////////////////////////////////////////////

Common::SharedPtr< RadiationPhysics > RadiationPhysicsHandler::getCellDistPtr
(CFuint stateID)
{
  //CFLog(INFO,"Cell; stateID: "<<stateID<<"\n");
  cf_assert(stateID<m_statesOwner.size() );
  cf_assert(m_statesOwner[stateID][0] != -1 );
  m_cellStateID = stateID;
  m_cellStateOwnerIdx = m_statesOwner[stateID][1];

  return m_radiationPhysics[ m_statesOwner[stateID][0] ];
}

//////////////////////////////////////////////////////////////////////////////

Common::SharedPtr< RadiationPhysics > RadiationPhysicsHandler::getWallDistPtr
(CFuint GhostStateID)
{
  using namespace COOLFluiD::Common;
  
  CFLog(DEBUG_MIN, "P" << PE::GetPE().GetRank("Default") << 
	" RadiationPhysicsHandler::getWallDistPtr() => stateID["
	<< GhostStateID << "] has owner [" << m_ghostStatesOwner[GhostStateID][0] << "]\n");
  
  cf_assert(GhostStateID<m_ghostStatesOwner.size() );
  cf_assert(m_ghostStatesOwner[GhostStateID][0] != -1 );
  m_ghostStateID = GhostStateID;
  m_ghostStateOwnerIdx  = m_ghostStatesOwner[GhostStateID][1];
  m_ghostStateWallGeoID = m_ghostStatesOwner[GhostStateID][2];
  
  return m_radiationPhysics[ m_ghostStatesOwner[GhostStateID][0] ];
}

//////////////////////////////////////////////////////////////////////////////

void RadiationPhysicsHandler::configureTRS()
{
  CFLog(VERBOSE, "RadiationPhysicsHandler::configureTRS() => START\n");
  
  Framework::DataHandle < Framework::State* , Framework::GLOBAL > states =
    m_sockets.states.getDataHandle();

  Framework::DataHandle < Framework::State*, Framework::LOCAL> gstates =
      m_sockets.gstates.getDataHandle();
  m_sockets.states.getDataHandle();

  m_statesOwner.resize(states.size(), std::vector<CFint>(2,-1) );
  m_ghostStatesOwner.resize(gstates.size(), std::vector<CFint>(3,-1) );
  vector<CFuint> buffer, buffer2;
  for(CFuint i=0; i<m_radiationPhysics.size();++i){
    TRStypeID bufferType = m_radiationPhysics[i]->getTRStypeID();
    CFLog(VERBOSE, "RadiationPhysicsHandler::configureTRS() => bufferType = " << bufferType << "\n");
    
    if (bufferType == MEDIUM) {
      m_radiationPhysics[i]->getCellStateIDs(buffer);
      CFLog(VERBOSE, "RadiationPhysicsHandler::configureTRS() => MEDIUM has size " << buffer.size() << "\n");
      
      const CFuint sizeBuffer = buffer.size();
      const CFuint sizeStates = states.size();
      m_statesOwner.resize(std::max( sizeBuffer , sizeStates ) );
      //std::cout<<"m_statesOwner: "<<sizeBuffer<<' '<<sizeStates<<std::endl;
      
      for(CFuint j=0; j<buffer.size();++j){
        m_statesOwner[ buffer[j] ][0] = i; //RadPhysics owner id of that state
        m_statesOwner[ buffer[j] ][1] = j; //state idx on the RadPhysics
      }
      
      buffer.clear();
    }
    // else{ // AL: OLD code, buggy in my opinion 
    else if (bufferType == WALL) {
      m_radiationPhysics[i]->getWallStateIDs(buffer, buffer2);
      CFLog(VERBOSE, "RadiationPhysicsHandler::configureTRS() => WALL has size " << buffer.size() << "\n");
      
      for (CFuint j=0; j<buffer.size();++j) {
        m_ghostStatesOwner[ buffer[j] ][0] = i;          //RadPhysics owner id of that g state
        m_ghostStatesOwner[ buffer[j] ][1] = j;          //state idx on the RadPhysics
        m_ghostStatesOwner[ buffer[j] ][2] = buffer2[j]; //Wall geo ID
      }
      
      buffer.clear();
    }
  }
  
/*
  CFLog(INFO,"m_statesOwner: \n");
  for(CFuint i = 0; i< m_statesOwner.size(); ++i){
      CFLog(INFO, m_statesOwner[i][0] <<' '<< m_statesOwner[i][1] <<'\n');
  }

  CFLog(INFO, " m_ghostStatesOwner \n");
  for(CFuint i = 0; i< m_ghostStatesOwner.size(); ++i){
    CFLog(INFO,   m_ghostStatesOwner[i][0] <<' '
               << m_ghostStatesOwner[i][1] <<' '
               << m_ghostStatesOwner[i][2] <<'\n');
  }
*/

  std::string bufferName;
  TRStypeID bufferType;
  for(CFuint i=0; i<m_radiationPhysics.size();++i){
    bufferName=m_radiationPhysics[i]->getTRSname();
    bufferType=m_radiationPhysics[i]->getTRStypeID();
    //std::cout<<"TRSNAME: "<<bufferName<<" TRSTypeID: "<<bufferType<<std::endl;

    switch (bufferType) {
    case WALL: // 0
      m_wallTRSnames.push_back(bufferName);
      break;
    case MEDIUM: // 1
      m_mediumTRSnames.push_back(bufferName);
      break;
    case BOUNDARY: // 2
      m_boundaryTRSnames.push_back(bufferName);
      break;
    default:
      CFLog(ERROR, "ERROR HANDLER: "<<bufferName<<": Undefined TRS type\n");
      break;
    }
  }
  
  CFLog(VERBOSE, "RadiationPhysicsHandler::configureTRS() => END\n");
}
  
//////////////////////////////////////////////////////////////////////////////

}
}
