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
#include "NavierStokes/RMSTurb.hh"
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

MethodCommandProvider<RMSTurb, DataProcessingData, NavierStokesModule>
rmsTurbProvider("RMSTurb"); //it is the name that used for the registration. It is "connected" to the pointer which shows this class

//////////////////////////////////////////////////////////////////////////////

void RMSTurb::defineConfigOptions(Config::OptionList& options)
{
  //options.addConfigOption< std::string >("OutputFile","Name of Output File.");
  //options.addConfigOption< CFuint >("SaveRate","Rate for saving the output file with RMSTurb.");
  //options.addConfigOption< bool >("AppendTime","Append time to file name.");
  //options.addConfigOption< bool >("AppendIter","Append Iteration# to file name.");
  options.addConfigOption< std::string >
  ("InitRMSTurbSocket","name of the input RMSTurb socket");
  options.addConfigOption< bool >
  ("ComNetVelocity","If true it will also compute the values derived by net velocity  (default = false)");
  options.addConfigOption< bool >
  ("ComAllTurbValues","If true it will also compute Turbulent intensities and energy (default = false)");
  options.addConfigOption< bool >("Restart","Should the initial rms be read from the CFmesh");
  options.addConfigOption< CFuint >
  ("nbSteps","Number of steps already known (this option implies that restart=true)");
  
}

//////////////////////////////////////////////////////////////////////////////

// The constructor initialize the parameters either by default or by taking the value from the CFcase file
RMSTurb::RMSTurb(const std::string& name) 
  : DataProcessingCom (name) , 
  //m_fhandle(), // open a file and give access to it for writing the results in the whole domain
  socket_states("states"), // this will give us access to the states through a pointer
  socket_rmsturb("rms"), // with this will store our results through a pointer
  m_dynamicSockets()
  
{
  addConfigOptionsTo(this);
  
  // open a file and give access to it for writing the results on the wall
  //m_fileRMS = Environment::SingleBehaviorFactory<Environment::FileHandlerOutput>::getInstance().create();

  //m_nameOutputFileRMS = "RMS_save.dat";
  //setParameter("OutputFile",&m_nameOutputFileRMS);
 
  //m_saveRateRMS = 1;
  //setParameter("SaveRate",&m_saveRateRMS);

  //m_appendIter = false;
  //setParameter("AppendIter",&m_appendIter);

  //m_appendTime = false;
  //setParameter("AppendTime",&m_appendTime);
  
  m_net = false;
  setParameter("ComNetVelocity",&m_net);
  
  m_comp = false;
  setParameter("ComAllTurbValues",&m_comp);
 
  m_restart = false;
  setParameter("Restart",&m_restart);

  m_InitSteps = 0;
  setParameter("nbSteps",&m_InitSteps);
  
  m_rmsInitSocketName = "Null";
   setParameter("InitRMSTurbSocket",&m_rmsInitSocketName);
  
}

//////////////////////////////////////////////////////////////////////////////

RMSTurb::~RMSTurb() //default destructor
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSource> >
RMSTurb::providesSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSource> > result;

  result.push_back(&socket_rmsturb);// vector with pointers showing the RMSTurb values

  return result;
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> > RMSTurb::needsSockets()
{ 
  // storage the pointers of the local data
  std::vector<Common::SafePtr<BaseDataSocketSink> > result = m_dynamicSockets.getAllSinkSockets();
  result.push_back(&socket_states); // vector with pointers showing the states P[rho U V W p...]

  return result;
  
  
}

//////////////////////////////////////////////////////////////////////////////

// it runs during the phase of global setup and allocate memory for this class, namely rms
void RMSTurb::setup()
{
  CFAUTOTRACE;
  // first call setup of parent class
  DataProcessingCom::setup(); 
  
  // with this command the variable states get access to the whole states
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  
  // this is the vector of the rms results
  DataHandle<CFreal> rmsturb = socket_rmsturb.getDataHandle();
  
    ///Initialization
  
    //number of the RMS values
    nbOpts = 0; 
    //number of the auxiliary additive values
    nbrms = 0; 

    // dimension of the problem
  const CFuint dim = PhysicalModelStack::getActive()->getDim(); 
  switch (dim) {
    case DIM_2D:
      if (m_comp)
      {
      nbOpts = 15;
      nbrms = 8; 
      // for the moment it needs to be true also if m_comp is true
      m_net = true;
      }
      else{
      nbOpts = 11;
      nbrms = 8; 
      }
      break;
    case DIM_3D:
      if (m_comp)
      {
      nbOpts = 19;
      nbrms = 10;
      // for the moment it needs to be true also if m_comp is true
      m_net = true;
      }
      else{
      nbOpts = 14;
      nbrms = 10;
      }
      break;
  } 
  
  // it substracts the number of the NetVelocity Variables
  if (!(m_net))
  {
    nbOpts -= 3;
    nbrms -= 2;
  }
  cf_assert(nbOpts>0);
  cf_assert(nbrms>0);
  
  // number of cells 
  const CFuint nbstates = states.size();
  
  // number of the auxiliary variables
  const CFuint totrms = nbstates*nbrms; 
  
  // number of extra variables for all the states
  const CFuint totstates = nbstates*nbOpts;
  
  /// physical data array 
  PhysicalModelStack::getActive()->getImplementor()->getConvectiveTerm()->
    resizePhysicalData(m_physicalData);
    
    // it allocates the correct size for the auxiliary variables
  rmsturb.resize(totstates); 

    // it allocates the correct size for the auxiliary variables
  m_rmsresult.resize(totrms,0.0); 
  
   //initialization of steps
  m_nbStep = m_InitSteps;

  // this is executed when we restart an application
  if (m_restart && m_dynamicSockets.sinkSocketExists(m_rmsInitSocketName))
  { 
    
    DataHandle<CFreal> initrms = m_dynamicSockets.getSocketSink<CFreal>(m_rmsInitSocketName)->getDataHandle();
    cf_assert(initrms.size() == rmsturb.size());
    for (CFuint i = 0; i < totstates; ++i) {
      //initialization of the rmsturb with previous stored data
      rmsturb[i] = initrms[i]; 
    }
    
     // common part (all except Urms Vrms Wrms and Utotrms) of the rmsturb and the m_rmsresults
     CFuint nbofCommonargs = 0;
     // common part (Urms Vrms Wrms and Utotrms) of the rmsturb and the m_rmsresults
     CFuint nbofCommonrmsargs = 0;
     // In case of m_comp = true it will also store the I and k variables
     CFuint nbofIntensities = 0;
     
    
    switch (dim) {
    case DIM_2D:
      nbofCommonargs = 5;
      nbofCommonrmsargs = 3;
      nbofIntensities = 4;
      
      break;
    case DIM_3D:
      nbofCommonargs = 6;
      nbofCommonrmsargs = 4;
      nbofIntensities = 5;
      break;
		}

  if (!(m_comp))
  {
    nbofIntensities = 0;
  }
  
  if (!(m_net))
  {
    nbofCommonargs -= 1;
    nbofCommonrmsargs -= 1;
  }
  cf_assert((nbofCommonargs + nbofCommonrmsargs)==nbrms);

      //Cause rmsturb are saved more data than m_rmsresult there must be a jump 
      const CFuint iJump = nbOpts - nbrms - nbofIntensities; 
      
      // this says how many jumbs should be done. One
      CFuint nbofparts = 0;  
   
      // this is the step for the accumulated variables of the RMS class
      CFuint jState = 0; 
      
      for (CFuint iState = 0; iState < totstates; ++iState, ++jState) {
	
	  m_rmsresult[jState] = rmsturb[iState]*m_nbStep; //
	  
	  if ( iState == ((nbofCommonargs - 1)+nbofparts*nbOpts))
	  {
	    
	    ++nbofparts;
	    iState += iJump;
	    
	    for (CFuint irms = 0; irms < nbofCommonrmsargs; ++irms) {
	      ++jState;
	      ++iState;
	      m_rmsresult[jState] =  rmsturb[iState]*rmsturb[iState]*m_nbStep;
	    }
	    
	    iState += nbofIntensities;
	    
	    cf_assert(jState == (iState-nbofparts*(iJump+nbofIntensities)));
	    
	    cf_assert((jState+1) == (nbrms*nbofparts));
	    cf_assert((iState+1) == (nbOpts*nbofparts));
	    
	  }
	}
	    cf_assert(nbofparts==nbstates);
      }
    
    
    //initialization of the reset parameter
    resetIter = 0;
    }
      
//////////////////////////////////////////////////////////////////////////////

void RMSTurb::execute()
{
  CFAUTOTRACE;
  cf_assert(nbOpts>0);
  cf_assert(nbrms>0);
  
  // if it should run the rms
  if (SubSystemStatusStack::getActive()->getNbIter() >= getMethodData().getStartIter()) {
    
    resetIter = 0;
    computeRMS();
  
  }
  
  else // we assume that if this check fails then we restart the computation and we need to reset the values
  {
    ++resetIter;
    
    // to run just the first time
    if (resetIter == 1) 
    {
    DataHandle<CFreal> rmsturb = socket_rmsturb.getDataHandle();
    DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
    
    m_nbStep = 0;
    
    const CFuint nbstates = states.size(); // number of cells
     
    for (CFuint iState = 0; iState < nbstates; ++iState) {
      
      for(CFuint iOpt = 0; iOpt < nbOpts; ++iOpt) 
	{const CFuint totistates = iState*nbOpts + iOpt;
	rmsturb[totistates] = 0.0;
	}
	for(CFuint irms = 0; irms < nbrms; ++irms) 
	{const CFuint totrms = iState*nbrms + irms;
	m_rmsresult[totrms] = 0.0;
	}
      }
    
    }
  }
  
}
      
      
//////////////////////////////////////////////////////////////////////////////

void RMSTurb::computeRMS()
{
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle<CFreal> rmsturb = socket_rmsturb.getDataHandle();
  
  // number of rms iterations
  ++m_nbStep;

  // number of cells
  const CFuint nbstates = states.size();
  
  // to compute the [P U V W..] of the current state
  SafePtr<ConvectiveVarSet> varSet = getMethodData().getUpdateVarSet();
  
  // dimension of the problem
  const CFuint dim = PhysicalModelStack::getActive()->getDim(); 
  
  for (CFuint iState = 0; iState < nbstates; ++iState) {
    // a pointer to current state
    State& curr_state = *states[iState];
   
    varSet->computePhysicalData(curr_state, m_physicalData);
    
    CFreal rho = m_physicalData[EulerTerm::RHO];
    CFreal Utot = m_physicalData[EulerTerm::V];
    CFreal u = m_physicalData[EulerTerm::VX];
    CFreal v = m_physicalData[EulerTerm::VY];
    CFreal w = m_physicalData[EulerTerm::VZ];
    CFreal p = m_physicalData[EulerTerm::P];
    //CFreal E = m_physicalData[EulerTerm::E];
    
    if (dim == 2) {
      
      CFuint totnbOpts = iState*nbOpts;
      CFuint totnbrms = iState*nbrms;
      
      // Calculation of rho mean
      m_rmsresult[totnbrms] += rho;
      rmsturb[totnbOpts] = m_rmsresult[totnbrms]/m_nbStep;
      
      // Calculation of p mean
      m_rmsresult[totnbrms + 1] += p;
      rmsturb[totnbOpts + 1] = m_rmsresult[totnbrms + 1]/m_nbStep;
      
      // Calculation of u mean
      m_rmsresult[totnbrms + 2] += u;
      rmsturb[totnbOpts + 2] = m_rmsresult[totnbrms + 2]/m_nbStep;
      
      // Calculation of v mean
      m_rmsresult[totnbrms + 3] += v;
      rmsturb[totnbOpts + 3] = m_rmsresult[totnbrms + 3]/m_nbStep;
      
      if (m_net)
      {
      // Calculation of Utot mean
      m_rmsresult[totnbrms + 4] += Utot;
      rmsturb[totnbOpts + 4] = m_rmsresult[totnbrms + 4]/m_nbStep;
      }
      else{
	totnbOpts -= 1; //we substract so as not to write out of the bound of the vector
	totnbrms -= 1;
      }
      
      // Calculation of uprim (u')
      rmsturb[totnbOpts + 5] = u - rmsturb[totnbOpts + 2];
      m_rmsresult[totnbrms + 5] += (rmsturb[totnbOpts + 5]*rmsturb[totnbOpts + 5]);
      
      // Calculation of vprim (v')
      rmsturb[totnbOpts + 6] = v - rmsturb[totnbOpts + 3];
      m_rmsresult[totnbrms + 6] += (rmsturb[totnbOpts + 6]*rmsturb[totnbOpts + 6]);
      
      if (m_net)
      {
      // Calculation of Utotprim (Utot')
      rmsturb[totnbOpts + 7] = Utot - rmsturb[totnbOpts + 4];
      m_rmsresult[totnbrms + 7] += (rmsturb[totnbOpts + 7]*rmsturb[totnbOpts + 7]);
      }
      else{
	totnbOpts -= 1; //we substract one so as not to write out of the bound of the vector
	cf_assert((totnbrms+7)%nbrms==0);
      }
      
      // Calculation of urms 
      rmsturb[totnbOpts + 8] = std::sqrt(m_rmsresult[totnbrms + 5]/m_nbStep);
      
      // Calculation of vrms 
      rmsturb[totnbOpts + 9] = std::sqrt(m_rmsresult[totnbrms + 6]/m_nbStep);
      
      if (m_net)
      {
      // Calculation of Utotrms 
      rmsturb[totnbOpts + 10] = std::sqrt(m_rmsresult[totnbrms + 7]/m_nbStep);
      }
      else{
	totnbOpts -= 1; //we substract so as not to write out of the bound of the vector
      }
      
      if (m_comp)
      {
      // Turbulence Intensities
      
      // In x-direction
      rmsturb[totnbOpts + 11] = rmsturb[totnbOpts + 8]/rmsturb[totnbOpts + 2];
      
      // In y-direction
      rmsturb[totnbOpts + 12] = rmsturb[totnbOpts + 9]/rmsturb[totnbOpts + 3];
      
      // Total
      rmsturb[totnbOpts + 13] = rmsturb[totnbOpts + 10]/rmsturb[totnbOpts + 4];
      
      // Turbulent kinetic energy
      rmsturb[totnbOpts + 14] = 0.5*(rmsturb[totnbOpts + 8]*rmsturb[totnbOpts + 8]+rmsturb[totnbOpts + 9]*rmsturb[totnbOpts + 9]);
      }
    }
    
    else if (dim == 3) {
      
      CFuint totnbOpts = iState*nbOpts;
      CFuint totnbrms = iState*nbrms;
      
      // Calculation of rho mean
      m_rmsresult[totnbrms] += rho;
      rmsturb[totnbOpts] = m_rmsresult[totnbrms]/m_nbStep;
      
      // Calculation of p mean
      m_rmsresult[totnbrms + 1] += p;
      rmsturb[totnbOpts + 1] = m_rmsresult[totnbrms + 1]/m_nbStep;
      
      // Calculation of u mean
      m_rmsresult[totnbrms + 2] += u;
      rmsturb[totnbOpts + 2] = m_rmsresult[totnbrms + 2]/m_nbStep;
      
      // Calculation of v mean
      m_rmsresult[totnbrms + 3] += v;
      rmsturb[totnbOpts + 3] = m_rmsresult[totnbrms + 3]/m_nbStep;
      
      // Calculation of w mean
      m_rmsresult[totnbrms + 4] += w;
      rmsturb[totnbOpts + 4] = m_rmsresult[totnbrms + 4]/m_nbStep;
      
      if (m_net)
      {
      // Calculation of Utot mean
      m_rmsresult[totnbrms + 5] += Utot;
      rmsturb[totnbOpts + 5] = m_rmsresult[totnbrms + 5]/m_nbStep;
      }
      else{
	totnbOpts -= 1; //we substract so as not to write out of the bound of the vector
	totnbrms -= 1;
      }
      
      // Calculation of uprim (u')
      rmsturb[totnbOpts + 6] = u - rmsturb[totnbOpts + 2];
      m_rmsresult[totnbrms + 6] += (rmsturb[totnbOpts + 6]*rmsturb[totnbOpts + 6]);
      
      // Calculation of vprim (v')
      rmsturb[totnbOpts + 7] = v - rmsturb[totnbOpts + 3];
      m_rmsresult[totnbrms + 7] += (rmsturb[totnbOpts + 7]*rmsturb[totnbOpts + 7]);
      
      // Calculation of vprim (w')
      rmsturb[totnbOpts + 8] = w - rmsturb[totnbOpts + 4];
      m_rmsresult[totnbrms + 8] += (rmsturb[totnbOpts + 8]*rmsturb[totnbOpts + 8]);
      
      
      if (m_net)
      {
      // Calculation of Utotprim (Utot')
      rmsturb[totnbOpts + 9] = Utot - rmsturb[totnbOpts + 5];
      m_rmsresult[totnbrms + 9] += (rmsturb[totnbOpts + 9]*rmsturb[totnbOpts + 9]);
      }
       else{
	totnbOpts -= 1; //we substract so as not to write out of the bound of the vector
      }
      
      // Calculation of urms 
      rmsturb[totnbOpts + 10] = std::sqrt(m_rmsresult[totnbrms + 6]/m_nbStep);
      
      // Calculation of vrms 
      rmsturb[totnbOpts + 11] = std::sqrt(m_rmsresult[totnbrms + 7]/m_nbStep);
      
      // Calculation of wrms 
      rmsturb[totnbOpts + 12] = std::sqrt(m_rmsresult[totnbrms + 8]/m_nbStep);
      
      if (m_net)
      {
      // Calculation of Utotrms 
      rmsturb[totnbOpts + 13] = std::sqrt(m_rmsresult[totnbrms + 9]/m_nbStep);
      }
       else{
	totnbOpts -= 1; //we substract so as not to write out of the bound of the vector
      }
      
      if (m_comp)
      {
      // Turbulence Intensities
      
      // In x-direction
      rmsturb[totnbOpts + 14] = rmsturb[totnbOpts + 10]/rmsturb[totnbOpts + 2];
      
      // In y-direction
      rmsturb[totnbOpts + 15] = rmsturb[totnbOpts + 11]/rmsturb[totnbOpts + 3];
      
      // In z-direction
      rmsturb[totnbOpts + 16] = rmsturb[totnbOpts + 12]/rmsturb[totnbOpts + 4];
      
      // Total
      rmsturb[totnbOpts + 17] = rmsturb[totnbOpts + 13]/rmsturb[totnbOpts + 5];
      
      // Turbulent kinetic energy
      rmsturb[totnbOpts + 18] = 0.5*(rmsturb[totnbOpts + 10]*rmsturb[totnbOpts + 10]+
				     rmsturb[totnbOpts + 11]*rmsturb[totnbOpts + 11]+
				     rmsturb[totnbOpts + 12]*rmsturb[totnbOpts + 12]);
      }
    }
  }
}

 
///@Note By default it will save the results in the CFmesh files and in the .plt with the other states values
/// If it needs to be stored in other file uncomment the below comments
  
//     Common::SafePtr<SubSystemStatus> ssys_status = SubSystemStatusStack::getActive();
//     const CFuint iter = ssys_status->getNbIter();
// 
//   if (!(iter % m_saveRateRMS))  { 
//     prepareOutputFileRMS();
//     cf_assert(m_fileRMS->isopen());
//     ofstream& fout = m_fileRMS->get();
// 
//  
//     for (CFuint iState = 0; iState < nbstates; ++iState) {
//       const CFreal x = (states[iState]->getCoordinates())[XX];
//       const CFreal y = (states[iState]->getCoordinates())[YY];
//       
//       fout << x << " " << y << " ";
//       if (dim ==3) {
// 	const CFreal z = (states[iState]->getCoordinates())[ZZ]; 
// 	fout << z << " ";
//       }
//       
//       for(CFuint iOpt = 0; iOpt < nbOpts; iOpt++) {
// 	const CFuint totistates = iState*nbOpts + iOpt;
// 	fout << rms[totistates] << " "; 
//       }
//       
//       fout << m_nbStep    << "\n ";
//     }
//     m_fileRMS->close();
//   } 
// }  

//////////////////////////////////////////////////////////////////////////////

void RMSTurb::unsetup()
{ 
  DataHandle<CFreal> rmsturb = socket_rmsturb.getDataHandle();
  rmsturb.resize(0);
  m_rmsresult.resize(0);
  
  // last call setup of parent class
  DataProcessingCom::unsetup(); 
}

//////////////////////////////////////////////////////////////////////////////

///@Note By default it will save the results in the CFmesh files and in the .plt with the other states values
/// If it needs to be stored in other file uncomment the below comments

// void RMS::prepareOutputFileRMS()
// {
//   using boost::filesystem::path;
//   
//   cf_assert (!m_fileRMS->isopen());
//   path file = Environment::DirPaths::getInstance().getResultsDir() / path (m_nameOutputFileRMS);
//   file = PathAppender::getInstance().appendAllInfo(file, m_appendIter, m_appendTime );
//   const CFuint dim = PhysicalModelStack::getActive()->getDim(); // dimension of the problem
//   const CFuint nbOpts =  m_rmsOpt.size(); // number of user options
//   
//   ofstream& fout = m_fileRMS->open(file);
//   
//   fout << "TITLE  =  RMS in the whole field" << "\n";
//   fout << "VARIABLES = x y ";
//   if (dim == 3) { fout << "z ";}
//   
//   for(CFuint iOpt = 0; iOpt < nbOpts; iOpt++) {
//     fout << m_rmsOpt[iOpt] << " "; 
//   }
//   
//   fout << "m_nbStep " << "\n";
//   
//   //fout    rhobar Ubar Vbar Wbar pbar rhoEbar devU m_nbStep " << "\n";
// }

//////////////////////////////////////////////////////////////////////////////
      
void RMSTurb::configure ( Config::ConfigArgs& args )
{
  DataProcessingCom::configure( args );
  
  if (m_rmsInitSocketName != "Null") {
    m_dynamicSockets.createSocketSink<CFreal>(m_rmsInitSocketName, true);
  }
  
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NavierStokes

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
