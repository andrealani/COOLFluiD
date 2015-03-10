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
rmsProvider("RMS"); //it is the name that used for the registration. It is "connected" to the pointer which shows this class

//////////////////////////////////////////////////////////////////////////////

void RMS::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< std::string >("OutputFile","Name of Output File.");
  options.addConfigOption< CFuint >("SaveRate","Rate for saving the output file with aerodynamic coefficients.");
  //options.addConfigOption< CFuint >("CompRate","Rate for saving the output file with aerodynamic coefficients.");
  options.addConfigOption< bool >("AppendTime","Append time to file name.");
  options.addConfigOption< bool >("AppendIter","Append Iteration# to file name.");
  options.addConfigOption< bool >("Restart","Should the initial rms be read from the CFmesh");
  options.addConfigOption< CFuint >("nbSteps","Number of steps already known (this option imply that restart=true)");
  
  //We add these options to let the user define which values will be time averaged
  options.addConfigOption< std::vector<std::string> >("rmsVars", "The list of variable names which will be time averaged.");
  
/*
  options.addConfigOption< bool >("DoAvgHeatFlux","(Default = False) If true it will compute the time averaged Heat Flux.");
  options.addConfigOption< bool >("DoAvgSkinFric","(Default = False) If true it will compute the time averaged Skin Friction.");
*/
}

//////////////////////////////////////////////////////////////////////////////

// The constructor initialize the parameters either by default or by taking the value from the CFcase file
RMS::RMS(const std::string& name) 
  : DataProcessingCom (name) , 
  m_fhandle(), // open a file and give access to it for writing the results in the whole domain
  socket_states("states"), // this will give us access to the states through a pointer
  socket_rms("rms") // with this will store our results through a pointer
  
{
  addConfigOptionsTo(this);
  
  // open a file and give access to it for writing the results on the wall
  m_fileRMS = Environment::SingleBehaviorFactory<Environment::FileHandlerOutput>::getInstance().create();

  m_nameOutputFileRMS = "RMS_save.dat";
  setParameter("OutputFile",&m_nameOutputFileRMS);

  //m_compRateRMS = 1;
  //setParameter("CompRate",&m_compRateRMS);
 
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
  
  // Initialization of the user choises
  
  m_rmsOpt = std::vector<std::string>();
  setParameter("rmsVars",&m_rmsOpt);
}

//////////////////////////////////////////////////////////////////////////////

RMS::~RMS() //default destructor
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> > RMS::needsSockets()
{ 
  std::vector<Common::SafePtr<BaseDataSocketSink> > result; // storage the pointers of the local data
  result.push_back(&socket_states); // vector with pointers showing the states P[rho U V W p...]
  return result;
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSource> > RMS::providesSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSource> > result;
  result.push_back(&socket_rms); // vector with pointers showing the rms
  return result;
}
      
//////////////////////////////////////////////////////////////////////////////

// it runs during the phase of global setup and allocate memory for this class, namely rms
void RMS::setup()
{
  CFAUTOTRACE;
  DataProcessingCom::setup(); // first call setup of parent class
  
  // with this command the variable states gets access to the whole states
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  const CFuint nbOpts =  m_rmsOpt.size(); // number of options of the user
  const CFuint nbstates = states.size(); // number of cells 
  const CFuint totstates = nbstates*nbOpts;
  
  /// it is the vector where our results will be stored
  DataHandle<CFreal> rms = socket_rms.getDataHandle();
  
  /// physical data array 
  PhysicalModelStack::getActive()->getImplementor()->getConvectiveTerm()->
    resizePhysicalData(m_physicalData);

  rms.resize(totstates); // it throws an error, perhaps it need resize(totstates) only
  m_rmsresult.resize(totstates,0.0); // it allocates the correct size for the vector of the result
  
  m_nbStep = m_InitSteps; //initialization of steps
  
  if (m_restart){ // this is executed when we restast an application
      DataHandle<CFreal> rms = socket_rms.getDataHandle();
      DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
      const CFuint nbOpts =  m_rmsOpt.size();
      const CFuint nbstates = states.size();

      for (CFuint iState = 0; iState < nbstates; ++iState) {
	for(CFuint iOpt = 0; iOpt < nbOpts; iOpt++) {
	  const CFuint totistates = nbstates*nbOpts + iOpt;
	  m_rmsresult[totistates] = rms[totistates];
	}
      }
    }
    
    resetIter = 0;
}
      
//////////////////////////////////////////////////////////////////////////////

void RMS::execute()
{
  CFAUTOTRACE;
  
  if (SubSystemStatusStack::getActive()->getNbIter() >= getMethodData().getStartIter()) {
    
    resetIter = 0;
    computeRMS();
  
  }
  
  else // we assume that if this check fails then we restart the computation and we need to reset the values
  {
    resetIter += 1;
    
    
    if (resetIter == 1) // to run just the first time
    {
    DataHandle<CFreal> rms = socket_rms.getDataHandle();
    DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
    
    m_nbStep = 0;
    
    const CFuint nbstates = states.size(); // number of cells
    const CFuint nbOpts =  m_rmsOpt.size(); // number of user options
    
    for (CFuint iState = 0; iState < nbstates; ++iState) {
      
      for(CFuint iOpt = 0; iOpt < nbOpts; iOpt++) 
	{m_rmsresult[iState*nbOpts+iOpt] = 0.0;
	rms[iState*nbOpts+iOpt] = 0.0;
	}
      }
    
    }
  }
  
}
      
      
//////////////////////////////////////////////////////////////////////////////

void RMS::computeRMS()
{
  Common::SafePtr<SubSystemStatus> subSysStatus = SubSystemStatusStack::getActive();
  const CFreal time = subSysStatus->getCurrentTime(); 
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle<CFreal> rms = socket_rms.getDataHandle();
  
  m_nbStep += 1; // number of rms iterations

  const CFuint nbstates = states.size(); // number of cells
  const CFuint nbOpts =  m_rmsOpt.size(); // number of user options
  
  SafePtr<ConvectiveVarSet> varSet = getMethodData().getUpdateVarSet();
  
  for (CFuint iState = 0; iState < nbstates; ++iState) {
    State& curr_state = *states[iState]; // a pointer to current state
   
    varSet->computePhysicalData(curr_state, m_physicalData); // compute the [P U V W..] of the current state
    
    CFreal rho = m_physicalData[EulerTerm::RHO];
    CFreal u = m_physicalData[EulerTerm::VX];
    CFreal v = m_physicalData[EulerTerm::VY];
    CFreal w = m_physicalData[EulerTerm::VZ];
    CFreal p = m_physicalData[EulerTerm::P];
    CFreal E = m_physicalData[EulerTerm::E];
    
    for(CFuint iOpt = 0; iOpt < nbOpts; iOpt++) {
      
     if (m_rmsOpt[iOpt] == "rho")
     {m_rmsresult[iState*nbOpts+iOpt] += rho;
      rms[iState*nbOpts+iOpt] = m_rmsresult[iState*nbOpts+iOpt];
      }
     
     else if (m_rmsOpt[iOpt] == "u")
     {m_rmsresult[iState*nbOpts+iOpt] += u;
      rms[iState*nbOpts+iOpt] = m_rmsresult[iState*nbOpts+iOpt];
      }
     
     else if (m_rmsOpt[iOpt] == "v")
     {m_rmsresult[iState*nbOpts+iOpt] += v;
      rms[iState*nbOpts+iOpt] = m_rmsresult[iState*nbOpts+iOpt];
      }
     
     else if (m_rmsOpt[iOpt] == "w")
     {m_rmsresult[iState*nbOpts+iOpt] += w;
      rms[iState*nbOpts+iOpt] = m_rmsresult[iState*nbOpts+iOpt];
      }
     
     else if (m_rmsOpt[iOpt] == "p")
       {
	 m_rmsresult[iState*nbOpts+iOpt] += p;
	 rms[iState*nbOpts+iOpt] = m_rmsresult[iState*nbOpts+iOpt];
       }
     
     else if (m_rmsOpt[iOpt] == "rhoE")
     {m_rmsresult[iState*nbOpts+iOpt] += (rho*E);
      cf_assert(iState*nbOpts+iOpt<rms.size());
      rms[iState*nbOpts+iOpt] = m_rmsresult[iState*nbOpts+iOpt];
      }
     
     else if (m_rmsOpt[iOpt] == "devU")
     {m_rmsresult[iState*nbOpts+iOpt] += u; // decouple the need to calculate also the ubar
      m_rmsresult[iState*nbOpts+iOpt] += (u - m_rmsresult[iState*nbOpts+iOpt]/m_nbStep)*(u - m_rmsresult[iState*nbOpts+iOpt]/m_nbStep);
      rms[iState*nbOpts+iOpt] = m_rmsresult[iState*nbOpts+iOpt];
      }
     
     //else ///put an error here
      
    } // iOpt for
 } // istate for 

 
    Common::SafePtr<SubSystemStatus> ssys_status = SubSystemStatusStack::getActive();
    const CFuint iter = ssys_status->getNbIter();

  if (!(iter % m_saveRateRMS))  { 
    prepareOutputFileRMS();
    cf_assert(m_fileRMS->isopen());
    ofstream& fout = m_fileRMS->get();
    const CFuint dim = PhysicalModelStack::getActive()->getDim(); // dimension of the problem
    CFreal theta;
    for (CFuint iState = 0; iState < nbstates; ++iState) {
      const CFreal x = (states[iState]->getCoordinates())[XX];
      const CFreal y = (states[iState]->getCoordinates())[YY];
      
      fout << x << " " << y << " ";
      if (dim ==3) {
	const CFreal z = (states[iState]->getCoordinates())[ZZ]; 
	fout << z << " ";
      }
      
      for(CFuint iOpt = 0; iOpt < nbOpts; iOpt++) {
	fout << m_rmsresult[iState*nbOpts+iOpt]/m_nbStep << " "; 
      }
      
      fout << m_nbStep    << "\n ";
    }
    m_fileRMS->close();
  } 
}  

//////////////////////////////////////////////////////////////////////////////

void RMS::unsetup()
{ 
  DataHandle<CFreal> rms = socket_rms.getDataHandle();
  rms.resize(0);
  m_rmsresult.resize(0);
  
  DataProcessingCom::unsetup(); // last call setup of parent class
}

//////////////////////////////////////////////////////////////////////////////

void RMS::prepareOutputFileRMS()
{
  using boost::filesystem::path;
  
  cf_assert (!m_fileRMS->isopen());
  path file = Environment::DirPaths::getInstance().getResultsDir() / path (m_nameOutputFileRMS);
  file = PathAppender::getInstance().appendAllInfo(file, m_appendIter, m_appendTime );
  const CFuint dim = PhysicalModelStack::getActive()->getDim(); // dimension of the problem
  const CFuint nbOpts =  m_rmsOpt.size(); // number of user options
  
  ofstream& fout = m_fileRMS->open(file);
  
  fout << "TITLE  =  RMS in the whole field" << "\n";
  fout << "VARIABLES = x y ";
  if (dim == 3) { fout << "z ";}
  
  for(CFuint iOpt = 0; iOpt < nbOpts; iOpt++) {
    fout << m_rmsOpt[iOpt] << " "; 
  }
  
  fout << "m_nbStep " << "\n";
  
  //fout    rhobar Ubar Vbar Wbar pbar rhoEbar devU m_nbStep " << "\n";
}

//////////////////////////////////////////////////////////////////////////////
      
void RMS::configure ( Config::ConfigArgs& args )
{
  DataProcessingCom::configure( args );
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NavierStokes

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
