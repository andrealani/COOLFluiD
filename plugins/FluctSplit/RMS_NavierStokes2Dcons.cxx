#include <iostream>
#include <string>
#include <fstream>

#include "Common/BadValueException.hh"
#include "Common/SafePtr.hh"
#include "Common/PE.hh"

#include "Environment/FileHandlerOutput.hh"
#include "Environment/SingleBehaviorFactory.hh"
#include "Environment/DirPaths.hh"

#include "Framework/PathAppender.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/NamespaceSwitcher.hh"

#include "FluctSplit/RMS_NavierStokes2Dcons.hh"
#include "FluctSplit/FluctSplitSpaceTime.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Environment;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

using namespace COOLFluiD::Physics::NavierStokes;
//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {


    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<RMS_NavierStokes2Dcons,
                      DataProcessingData,
		      FluctSplitSpaceTimeModule>
rmsns2dconsProvider("RMS_NavierStokes2Dcons");

//////////////////////////////////////////////////////////////////////////////

void RMS_NavierStokes2Dcons::defineConfigOptions(Config::OptionList& options)
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

RMS_NavierStokes2Dcons::RMS_NavierStokes2Dcons( const std::string& name) :
  DataProcessingCom(name),
  m_fhandle(),
  socket_nodes("nodes"),
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

RMS_NavierStokes2Dcons::~RMS_NavierStokes2Dcons()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
RMS_NavierStokes2Dcons::needsSockets()
{ 
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;
  result.push_back(&socket_nodes);
  result.push_back(&socket_states); 
  result.push_back(&socket_rms);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

void RMS_NavierStokes2Dcons::setup()
{
  CFAUTOTRACE;
 DataProcessingCom::setup(); // first call setup of parent class

  m_varSet->setup();

  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  const CFuint nbstates = states.size();
  m_rhobar.resize(nbstates,0.0);
  m_Ubar.resize(nbstates,0.0);
  m_Vbar.resize(nbstates,0.0);
  m_UVbar.resize(nbstates,0.0);
  m_UUbar.resize(nbstates,0.0) ;
  m_VVbar.resize(nbstates,0.0);
  m_pbar.resize(nbstates,0.0);
  m_rhoEbar.resize(nbstates,0.0);
  m_rms.resize(nbstates,0.0);
  m_devRho.resize(nbstates,0.0);
  m_nbStep = m_InitSteps;


}
//////////////////////////////////////////////////////////////////////////////

void RMS_NavierStokes2Dcons::execute()
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
       
         m_rhobar[iState] = rms[10*iState];
   	m_Ubar[iState] =  rms[10*iState+1];
   	m_Vbar[iState] =  rms[10*iState+2];
   	m_UVbar[iState] = rms[10*iState+3];
   	m_UUbar[iState] =  rms[10*iState+4];
   	m_VVbar[iState] = rms[10*iState+5];
   	m_pbar[iState] = rms[10*iState+6] ;
   	m_rhoEbar[iState] = rms[10*iState+7];
   	m_rms[iState] = rms[10*iState+8];
   	m_devRho[iState] = rms[10*iState+9];
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
void RMS_NavierStokes2Dcons::computeRMS(bool save){


   Common::SafePtr<SubSystemStatus> subSysStatus = SubSystemStatusStack::getActive();
  const CFreal time = subSysStatus->getCurrentTime();
 DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle<CFreal> rms = socket_rms.getDataHandle();
 m_nbStep += 1;
 const CFuint nbstates = states.size();

 for (CFuint iState = 0; iState < nbstates; ++iState) {
    const CFreal x = (states[iState]->getCoordinates())[XX];
    const CFreal y = (states[iState]->getCoordinates())[YY];
    State& curr_state = *states[iState];
   
    CFreal rho = (*states[iState])[0];
    CFreal rhou = (*states[iState])[1];
    CFreal rhov = (*states[iState])[2];
    CFreal rhoE = (*states[iState])[3];
   
  CFreal gamma = m_varSet->getModel()->getGamma();
  CFreal u = rhou/rho;
  CFreal v = rhov/rho;
  CFreal p = (gamma-1.0)*(rhoE-0.5*rho*(u*u+v*v)) ;

  m_rhobar[iState] += rho;
  m_Ubar[iState] += u ;
  m_Vbar[iState] += v ;
  m_UVbar[iState] += u*v ;
  m_UUbar[iState] += u*u  ;
  m_VVbar[iState] += v*v ;
  m_pbar[iState] += p ;
  m_rhoEbar[iState] += rhoE;
  m_rms[iState] += (p - m_pbar[iState]/m_nbStep)*(p - m_pbar[iState]/m_nbStep);
  m_devRho[iState] += (rho - m_rhobar[iState]/m_nbStep)*(rho - m_rhobar[iState]/m_nbStep);

  rms[10*iState] = m_rhobar[iState];
   rms[10*iState+1] = m_Ubar[iState];
   rms[10*iState+2] = m_Vbar[iState];
   rms[10*iState+3] = m_UVbar[iState];
   rms[10*iState+4] = m_UUbar[iState];
   rms[10*iState+5] = m_VVbar[iState];
   rms[10*iState+6] = m_pbar[iState];
   rms[10*iState+7] = m_rhoEbar[iState];
   rms[10*iState+8] = m_rms[iState];
   rms[10*iState+9] = m_devRho[iState];
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
	
	if ((x > 0.0)&&(y>0.0))
	  theta = atan(y/x);
	if ((x<0.0)&&(y>0.0))
	  theta = 3.14159-atan(-x/y);
	if ((x<0.0)&&(y<0.0))
	  theta = 3.14159+atan(x/y);
	if ((x>0.0)&&(y>0.0))
	  theta = -atan(x/(-y));

	fout
	  << x << " "
	  << y << " "
	  << m_rhobar[iState]/m_nbStep << " "
	  << m_Ubar[iState]/m_nbStep   << " "
	  << m_Vbar[iState]/m_nbStep << " "
	  <<  m_UVbar[iState]/m_nbStep << " "
	  <<  m_UUbar[iState]/m_nbStep << " "
	  <<  m_VVbar[iState]/m_nbStep << " "
	  <<  m_pbar[iState]/m_nbStep << " "
	  <<  m_rhoEbar[iState]/m_nbStep << " "
	  <<  sqrt(m_rms[iState])/m_nbStep << " "
	  <<  sqrt(m_devRho[iState])/m_nbStep << " "
	  << theta << " "
	  << m_nbStep    << "\n ";
 
    
      }
      m_fileRMS->close();
    }

}
  

//////////////////////////////////////////////////////////////////////////////

void RMS_NavierStokes2Dcons::unsetup()
{ 
  DataHandle<CFreal> rms = socket_rms.getDataHandle();

  rms.resize(0);


}

//////////////////////////////////////////////////////////////////////////////
void RMS_NavierStokes2Dcons::prepareOutputFileRMS()
{

  using boost::filesystem::path;

  cf_assert (!m_fileRMS->isopen());
  path file = Environment::DirPaths::getInstance().getResultsDir() / path (m_nameOutputFileRMS);
  file = PathAppender::getInstance().appendAllInfo(file, m_appendIter, m_appendTime );


  ofstream& fout = m_fileRMS->open(file);


  fout << "TITLE  =  RMS in the whole field" << "\n";
  fout << "VARIABLES = x y rhobar Ubar Vbar UVbar UUbar VVbar pbar rhoEbar prms devRho theta m_nbStep " << "\n";
  
}

//////////////////////////////////////////////////////////////////////////////

void RMS_NavierStokes2Dcons::configure ( Config::ConfigArgs& args )
{
  DataProcessingCom::configure( args );
  std::string name = getMethodData().getNamespace();
  Common::SafePtr<Namespace> nsp = NamespaceSwitcher::getInstance().getNamespace(name);

  Common::SafePtr<PhysicalModel> physModel = PhysicalModelStack::getInstance().getEntryByNamespace(nsp);
  std::string varSetName = "Euler2DCons";
  m_varSet.reset((Environment::Factory<ConvectiveVarSet>::getInstance().getProvider(varSetName)->
    create(physModel->getImplementor()->getConvectiveTerm())).d_castTo<Physics::NavierStokes::Euler2DCons>());

  cf_assert(m_varSet.isNotNull());


}

//////////////////////////////////////////////////////////////////////////////

    } //namespace FluctSplit

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

