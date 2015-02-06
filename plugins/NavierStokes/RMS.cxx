#include <iostream>
#include <string>
#include <fstream>

#include "Common/BadValueException.hh"
#include "Common/SafePtr.hh"

#include "Environment/FileHandlerOutput.hh"
#include "Environment/SingleBehaviorFactory.hh"
#include "Environment/DirPaths.hh"

#include "Framework/PhysicalModel.hh"
#include "Framework/PathAppender.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/NamespaceSwitcher.hh"
#include "Framework/ConvectiveVarSet.hh"

#include "NavierStokes/NavierStokes.hh"
#include "NavierStokes/RMS.hh"
#include "NavierStokes/EulerTerm.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Environment;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace NavierStokes {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<RMS, DataProcessingData, NavierStokesModule>
rmsProvider("RMS");

//////////////////////////////////////////////////////////////////////////////

void RMS::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< std::string >("OutputFile","Name of Output File.");
  options.addConfigOption< CFuint >("SaveRate","Rate for saving the output file with aerodynamic coefficients.");
  options.addConfigOption< CFuint >("CompRate","Rate for saving the output file with aerodynamic coefficients.");
  options.addConfigOption< bool >("AppendTime","Append time to file name.");
  options.addConfigOption< bool >("AppendIter","Append Iteration# to file name.");
  options.addConfigOption< bool >("Restart","Should the initial rms be read from the CFmesh");
  options.addConfigOption< CFuint >("nbSteps","Number of steps already known (this option imply that restart=true)");
}

//////////////////////////////////////////////////////////////////////////////

RMS::RMS( const std::string& name) :
  DataProcessingCom(name),
  m_fhandle(),
  socket_states("states"),
  socket_rms("rms")
{
  m_fileRMS = Environment::SingleBehaviorFactory<Environment::FileHandlerOutput>::getInstance().create();

  addConfigOptionsTo(this);
  m_nameOutputFileRMS = "RMS_save.dat";
  setParameter("OutputFile",&m_nameOutputFileRMS);

  m_compRateRMS = 1;
  setParameter("CompRate",&m_compRateRMS);
 
  m_saveRateRMS = 1;
  setParameter("SaveRate",&m_saveRateRMS);

  m_appendIter = false;
  setParameter("AppendIter",&m_appendIter);

  m_appendTime = false;
  setParameter("AppendTime",&m_appendTime);
 
  m_restart = false;
  setParameter("Restart",&m_restart);

  m_InitSteps = 0;
  setParameter("nbSteps",&m_InitSteps);
}

//////////////////////////////////////////////////////////////////////////////

RMS::~RMS()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> > RMS::needsSockets()
{ 
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;
  result.push_back(&socket_states); 
  result.push_back(&socket_rms);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

void RMS::setup()
{
  CFAUTOTRACE;
 DataProcessingCom::setup(); // first call setup of parent class

  m_varSet->setup();

  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  const CFuint nbstates = states.size();
  
  /// physical data array 
  PhysicalModelStack::getActive()->getImplementor()->getConvectiveTerm()->
    resizePhysicalData(m_physicalData);

  
  m_rhobar.resize(nbstates,0.0);
  m_Ubar.resize(nbstates,0.0);
  m_Vbar.resize(nbstates,0.0);
  m_Wbar.resize(nbstates,0.0);
  m_pbar.resize(nbstates,0.0);
  m_rhoEbar.resize(nbstates,0.0);
  m_devU.resize(nbstates,0.0);

  m_nbStep = m_InitSteps;


}
//////////////////////////////////////////////////////////////////////////////

void RMS::execute()
{
  Common::SafePtr<SubSystemStatus> ssys_status = SubSystemStatusStack::getActive();

  const CFuint iter = ssys_status->getNbIter();

   if (iter == 0){
     if (m_restart){
        DataHandle<CFreal> rms = socket_rms.getDataHandle();
        DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
        m_nbStep = m_InitSteps;
         const CFuint nbstates = states.size();
         for (CFuint iState = 0; iState < nbstates; ++iState) {
       
         m_rhobar[iState] = rms[7*iState];
   	m_Ubar[iState] =  rms[7*iState+1];
   	m_Vbar[iState] =  rms[7*iState+2];
   	m_Wbar[iState] = rms[7*iState+3];
   	m_pbar[iState] = rms[7*iState+4] ;
   	m_rhoEbar[iState] = rms[7*iState+5];
   	m_devU[iState] = rms[7*iState+6];
     }
   }
   }

  // compute global integrated coefficients
  if(!(iter % m_compRateRMS)){
    if(!(iter % m_saveRateRMS))  { computeRMS(true); }
    else{ computeRMS(false); }
  }
  
  
}


//////////////////////////////////////////////////////////////////////////////

void RMS::computeRMS(bool save)
{
  Common::SafePtr<SubSystemStatus> subSysStatus = SubSystemStatusStack::getActive();
  const CFreal time = subSysStatus->getCurrentTime();
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle<CFreal> rms = socket_rms.getDataHandle();
 m_nbStep += 1;
 const CFuint nbstates = states.size();

 for (CFuint iState = 0; iState < nbstates; ++iState) {
    const CFreal x = (states[iState]->getCoordinates())[XX];
    const CFreal y = (states[iState]->getCoordinates())[YY];
    const CFreal z = (states[iState]->getCoordinates())[ZZ];
    State& curr_state = *states[iState];
    
    m_varSet->computePhysicalData(curr_state, m_physicalData);
    
    CFreal rho = (*states[iState])[0];
    CFreal rhou = (*states[iState])[1];
    CFreal rhov = (*states[iState])[2];
    CFreal rhow = (*states[iState])[3];
    CFreal rhoE = (*states[iState])[4];
    
    CFreal u = m_physicalData[EulerTerm::VX];
    CFreal v = rhov/rho;
    CFreal w = rhow/rho;
    CFreal p = m_physicalData[EulerTerm::P];
      
  m_rhobar[iState] += rho;
  m_Ubar[iState] += u ;
  m_Vbar[iState] += v ;
  m_Wbar[iState] += w ;
  m_pbar[iState] += p ;
  m_rhoEbar[iState] += rhoE;
  m_devU[iState] += (u - m_Ubar[iState]/m_nbStep)*(u - m_Ubar[iState]/m_nbStep);

  rms[7*iState] = m_rhobar[iState];
   rms[7*iState+1] = m_Ubar[iState];
   rms[7*iState+2] = m_Vbar[iState];
   rms[7*iState+3] = m_Wbar[iState];
   rms[7*iState+4] = m_pbar[iState];
   rms[7*iState+5] = m_rhoEbar[iState];
   rms[7*iState+6] = m_devU[iState];
 }


  if (save) 
    { 
      prepareOutputFileRMS();
       cf_assert(m_fileRMS->isopen());
      ofstream& fout = m_fileRMS->get();
      CFreal theta;
      for (CFuint iState = 0; iState < nbstates; ++iState) {
	const CFreal x = (states[iState]->getCoordinates())[XX];
	const CFreal y = (states[iState]->getCoordinates())[YY];
	const CFreal z = (states[iState]->getCoordinates())[ZZ];
	

	fout
	  << x << " "
	  << y << " "
	  << z << " "
	  << m_rhobar[iState]/m_nbStep << " "
	  << m_Ubar[iState]/m_nbStep   << " "
	  << m_Vbar[iState]/m_nbStep << " "
	  <<  m_Wbar[iState]/m_nbStep << " "
	  <<  m_pbar[iState]/m_nbStep << " "
	  <<  m_rhoEbar[iState]/m_nbStep << " "
	  <<  sqrt(m_devU[iState])/m_nbStep << " "
	  << m_nbStep    << "\n ";
 
    
      }
      m_fileRMS->close();
    }

}
  

//////////////////////////////////////////////////////////////////////////////

void RMS::unsetup()
{ 
  DataHandle<CFreal> rms = socket_rms.getDataHandle();

  rms.resize(0);


}

//////////////////////////////////////////////////////////////////////////////
void RMS::prepareOutputFileRMS()
{

  using boost::filesystem::path;

  cf_assert (!m_fileRMS->isopen());
  path file = Environment::DirPaths::getInstance().getResultsDir() / path (m_nameOutputFileRMS);
  file = PathAppender::getInstance().appendAllInfo(file, m_appendIter, m_appendTime );


  ofstream& fout = m_fileRMS->open(file);


  fout << "TITLE  =  RMS in the whole field" << "\n";
  fout << "VARIABLES = x y z rhobar Ubar Vbar Wbar pbar rhoEbar devU m_nbStep " << "\n";
  
}

//////////////////////////////////////////////////////////////////////////////

void RMS::configure ( Config::ConfigArgs& args )
{
  DataProcessingCom::configure( args );
  std::string name = getMethodData().getNamespace();
  Common::SafePtr<Namespace> nsp = NamespaceSwitcher::getInstance().getNamespace(name);
  
  Common::SafePtr<PhysicalModel> physModel = PhysicalModelStack::getInstance().getEntryByNamespace(nsp);
  std::string varSetName = "Euler3DCons";  // THIS NAME HAS TO COME AS AN OPTION
  m_varSet.reset((Environment::Factory<ConvectiveVarSet>::getInstance().getProvider(varSetName)->
		  create(physModel->getImplementor()->getConvectiveTerm())));
  
  cf_assert(m_varSet.isNotNull());
}

//////////////////////////////////////////////////////////////////////////////

 } // namespace NavierStokes

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
