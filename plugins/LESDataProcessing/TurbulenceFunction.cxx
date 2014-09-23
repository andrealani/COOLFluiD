
#include "TurbulenceFunction.hh"
#include "Common/NotImplementedException.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {
	
		namespace LESDataProcessing {

//////////////////////////////////////////////////////////////////////////////

void TurbulenceFunction::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< bool >("ExtrapolateToNodal","Flag that tells if the turbulence function will be extrapolated to nodal values");
  options.addConfigOption< bool >("FirstTimeCreation","Flag that tells to create source socket instead of sink sockets. (default=false)");
}

//////////////////////////////////////////////////////////////////////////////

TurbulenceFunction::TurbulenceFunction(const std::string& name) :
  Framework::MethodStrategy<LESProcessingData>(name),
  m_socket(CFNULL)
{
  addConfigOptionsTo(this);
  // Setting default configurations here.
  m_extrapolate = false;
  setParameter("ExtrapolateToNodal",&m_extrapolate);
  
  m_firstTimeCreation = true;
  setParameter("FirstTimeCreation",&m_firstTimeCreation);

  
}

//////////////////////////////////////////////////////////////////////////////


TurbulenceFunction::~TurbulenceFunction()
{
}

//////////////////////////////////////////////////////////////////////////////

void TurbulenceFunction::computeAverage(const CFreal& oldWeight, const CFreal& newWeight, const RealVector& states, std::vector<RealVector>& gradients, const CFuint& iCell)
{
  throw Common::NotImplementedException (FromHere(),"TurbulenceFunction::computeAverage() needs to be overloaded");
}

//////////////////////////////////////////////////////////////////////////////

void TurbulenceFunction::configure ( Config::ConfigArgs& args )
{
  Framework::MethodStrategy<LESProcessingData>::configure(args);
  setSocketName();
    
  if (m_firstTimeCreation) {
    m_sockets.createSocketSource<CFreal>(m_socketName);
    CFLog (INFO,"   +++ Created source socket for " << m_socketName << "\n");
  } 
  // Make sink socket if source socket already exists
  else {
    m_sockets.createSocketSink<CFreal>(m_socketName);
    CFLog (INFO,"   +++ Created sink socket for " << m_socketName << "\n");
  }
  
  
  if (extrapolate()) {  // Store socket in nodal states
    m_sockets.createSocketSink<RealVector>("nstates");
  }
  else { // Store socket in states
    m_globalSockets.createSocketSink<Framework::State*>("states");
  }
}

//////////////////////////////////////////////////////////////////////////////

void TurbulenceFunction::setup()
{
  CFAUTOTRACE;
  Framework::MethodStrategy<LESProcessingData>::setup();
  if (extrapolate()) {
    Framework::DataHandle<RealVector> nodalStates = m_sockets.getSocketSink<RealVector>("nstates")->getDataHandle();
    m_nbStates = nodalStates.size();
  }
  else {
    Framework::DataHandle<Framework::State*,Framework::GLOBAL> states = m_globalSockets.getSocketSink<Framework::State*>("states")->getDataHandle();
    m_nbStates = states.size();
  }
    
  if (m_firstTimeCreation) {
    m_socket = m_sockets.getSocketSource<CFreal>(m_socketName)->getDataHandle();  
  } 
  // Use sink socket if source socket already exists
  else {
    m_socket = m_sockets.getSocketSink<CFreal>(m_socketName)->getDataHandle();  
  }

}

//////////////////////////////////////////////////////////////////////////////
	
	  } // end of namespace LESDataProcessing
		
  } // namespace Numerics

} // namespace COOLFluiD
