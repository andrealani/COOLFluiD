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
  //options.addConfigOption< CFuint >("SaveRate","Rate for saving the output file with aerodynamic coefficients.");
  //options.addConfigOption< CFuint >("CompRate","Rate for saving the output file with aerodynamic coefficients.");
  //options.addConfigOption< bool >("AppendTime","Append time to file name.");
  //options.addConfigOption< bool >("AppendIter","Append Iteration# to file name.");
  //options.addConfigOption< bool >("ComNetVelocity","If true it will also compute the values derived by net velocity  (default = false)");
  options.addConfigOption< bool >("ComAllTurbValues","If true it will also compute Turbulent intensities and energy (default = false)");
  options.addConfigOption< bool >("Restart","Should the initial rms be read from the CFmesh");
  options.addConfigOption< CFuint >("nbSteps","Number of steps already known (this option imply that restart=true)");
  
  //We add these options to let the user define which values will be time averaged
  //options.addConfigOption< std::vector<std::string> >("rmsVars", "The list of variable names which will be time averaged.");
}

//////////////////////////////////////////////////////////////////////////////

// The constructor initialize the parameters either by default or by taking the value from the CFcase file
RMSTurb::RMSTurb(const std::string& name) 
  : DataProcessingCom (name) , 
  //m_fhandle(), // open a file and give access to it for writing the results in the whole domain
  socket_states("states"), // this will give us access to the states through a pointer
  socket_rms("rms") // with this will store our results through a pointer
  
{
  addConfigOptionsTo(this);
  
  // open a file and give access to it for writing the results on the wall
  //m_fileRMS = Environment::SingleBehaviorFactory<Environment::FileHandlerOutput>::getInstance().create();

  //m_nameOutputFileRMS = "RMS_save.dat";
  //setParameter("OutputFile",&m_nameOutputFileRMS);

  //m_compRateRMS = 1;
  //setParameter("CompRate",&m_compRateRMS);
 
  //m_saveRateRMS = 1;
  //setParameter("SaveRate",&m_saveRateRMS);

  //m_appendIter = false;
  //setParameter("AppendIter",&m_appendIter);

  //m_appendTime = false;
  //setParameter("AppendTime",&m_appendTime);
  
  m_comp = false;
  setParameter("ComAllTurbValues",&m_comp);
 
  m_restart = false;
  setParameter("Restart",&m_restart);

  m_InitSteps = 0;
  setParameter("nbSteps",&m_InitSteps);
  
  // Initialization of the user choises
  
  //m_rmsOpt = std::vector<std::string>();
  //setParameter("rmsVars",&m_rmsOpt);
}

//////////////////////////////////////////////////////////////////////////////

RMSTurb::~RMSTurb() //default destructor
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> > RMSTurb::needsSockets()
{ 
  std::vector<Common::SafePtr<BaseDataSocketSink> > result; // storage the pointers of the local data
  result.push_back(&socket_states); // vector with pointers showing the states P[rho U V W p...]
  return result;
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSource> > RMSTurb::providesSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSource> > result;
  result.push_back(&socket_rms); // vector with pointers showing the rms
  return result;
}
      
//////////////////////////////////////////////////////////////////////////////

// it runs during the phase of global setup and allocate memory for this class, namely rms
void RMSTurb::setup()
{
  CFAUTOTRACE;
  DataProcessingCom::setup(); // first call setup of parent class
  
  // with this command the variable states get access to the whole states
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  //const CFuint nbOpts =  m_rmsOpt.size(); // number of options of the user
  
   //CFuint nbOpts = 0;
   //CFuint nbrms = 0;

  const CFuint dim = PhysicalModelStack::getActive()->getDim(); // dimension of the problem
  switch (dim) {
    case DIM_2D:
      if (m_comp)
      {
      nbOpts = 15;//it is the number of the turbulent variables that will be saved
      nbrms = 8; // it is the number of the auxiliary variables used for the calculation of the turbulent variables
      }
      else{
      nbOpts = 11;//it is the number of the turbulent variables that will be saved
      nbrms = 8; // it is the number of the auxiliary variables used for the calculation of the turbulent variables
      }
      break;
    case DIM_3D:
      if (m_comp)
      {
      nbOpts = 19;//it is the number of the turbulent variables that will be saved
      nbrms = 10; // it is the number of the auxiliary variables used for the calculation of the turbulent variables
      }
      else{
      nbOpts = 14;//it is the number of the turbulent variables that will be saved
      nbrms = 10; // it is the number of the auxiliary variables used for the calculation of the turbulent variables
      }
      break;
  }
  
  const CFuint nbstates = states.size(); // number of cells 
  const CFuint totstates = nbstates*nbOpts; // number of extra variables for the states
  const CFuint totrms = nbstates*nbrms; // number of the auxiliary variables
  
  /// it is the vector where our results will be stored
  DataHandle<CFreal> rms = socket_rms.getDataHandle();
  
  /// physical data array 
  PhysicalModelStack::getActive()->getImplementor()->getConvectiveTerm()->
    resizePhysicalData(m_physicalData);

  rms.resize(totstates); 
  m_rmsresult.resize(totrms,0.0); // it allocates the correct size for the auxiliary variables
  
  m_nbStep = m_InitSteps; //initialization of steps
  
  if (m_restart){ // this is executed when we restart an application
    
     CFuint nbofCommonargs = 0;
    
    switch (dim) {
    case DIM_2D:
      nbofCommonargs = 5;
      break;
    case DIM_3D:
      nbofCommonargs = 6;
      break;
		}

      const CFuint iJump = nbOpts - nbofCommonargs;
      const CFuint jJump = nbrms - nbofCommonargs;
      CFuint nbofparts = 0;
      
      cf_assert((nbofCommonargs+iJump)==nbOpts);
      cf_assert((nbofCommonargs+jJump)==nbrms);
      
      CFuint jState = 0;
      for (CFuint iState = 0; iState < totstates; ++iState) {
	  ++jState;
	  m_rmsresult[jState] = rms[iState]*m_nbStep;
	  if ( iState == ((nbofCommonargs - 1)+nbofparts*nbOpts))
	  {
	    ++nbofparts;
	    iState += iJump;
	    jState += jJump;
      cf_assert(iState<totstates);
	  }
	}
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
    
    
    if (resetIter == 1) // to run just the first time
    {
    DataHandle<CFreal> rms = socket_rms.getDataHandle();
    DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
    
    m_nbStep = 0;
    
    const CFuint nbstates = states.size(); // number of cells
     //CFuint nbOpts = 0;
     //CFuint nbrms = 0;
 /*   const CFuint dim = PhysicalModelStack::getActive()->getDim(); // dimension of the problem
  switch (dim) {
    case DIM_2D:
      nbOpts = 15;//it is the number of the turbulent variables that will be saved
      nbrms = 8; // it is the number of the auxiliary variables used for the calculation of the turbulent variables
      break;
    case DIM_3D:
      nbOpts = 19;
      nbrms = 10;
      break;
  }
*/
    for (CFuint iState = 0; iState < nbstates; ++iState) {
      
      for(CFuint iOpt = 0; iOpt < nbOpts; ++iOpt) 
	{const CFuint totistates = iState*nbOpts + iOpt;
	rms[totistates] = 0.0;
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
  //Common::SafePtr<SubSystemStatus> subSysStatus = SubSystemStatusStack::getActive();
  //const CFreal time = subSysStatus->getCurrentTime(); 
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle<CFreal> rms = socket_rms.getDataHandle();
  
  ++m_nbStep; // number of rms iterations

  const CFuint nbstates = states.size(); // number of cells
  //const CFuint nbOpts =  m_rmsOpt.size(); // number of user options
  
  SafePtr<ConvectiveVarSet> varSet = getMethodData().getUpdateVarSet();
  
  const CFuint dim = PhysicalModelStack::getActive()->getDim(); // dimension of the problem
  
  for (CFuint iState = 0; iState < nbstates; ++iState) {
    State& curr_state = *states[iState]; // a pointer to current state
   
    varSet->computePhysicalData(curr_state, m_physicalData); // compute the [P U V W..] of the current state
    
    CFreal rho = m_physicalData[EulerTerm::RHO];
    CFreal Utot = m_physicalData[EulerTerm::V];
    CFreal u = m_physicalData[EulerTerm::VX];
    CFreal v = m_physicalData[EulerTerm::VY];
    CFreal w = m_physicalData[EulerTerm::VZ];
    CFreal p = m_physicalData[EulerTerm::P];
    //CFreal E = m_physicalData[EulerTerm::E];
    
    if (dim == 2) {
      
      //const CFuint nbOpts = 15;//it is the number of the turbulent variables that will be saved
      //const CFuint nbrms = 8; // it is the number of the auxiliary variables used for the calculation of the turbulent variables

      const CFuint totnbOpts = iState*nbOpts;
      const CFuint totnbrms = iState*nbrms;
      
      // Calculation of rho mean
      m_rmsresult[totnbrms] += rho;
      rms[totnbOpts] = m_rmsresult[totnbrms]/m_nbStep;
      
      // Calculation of p mean
      m_rmsresult[totnbrms + 1] += p;
      rms[totnbOpts + 1] = m_rmsresult[totnbrms + 1]/m_nbStep;
      
      // Calculation of u mean
      m_rmsresult[totnbrms + 2] += u;
      rms[totnbOpts + 2] = m_rmsresult[totnbrms + 2]/m_nbStep;
      
      // Calculation of v mean
      m_rmsresult[totnbrms + 3] += v;
      rms[totnbOpts + 3] = m_rmsresult[totnbrms + 3]/m_nbStep;
      
      // Calculation of Utot mean
      m_rmsresult[totnbrms + 4] += Utot;
      rms[totnbOpts + 4] = m_rmsresult[totnbrms + 4]/m_nbStep;
      
      // Calculation of uprim (u')
      rms[totnbOpts + 5] = u - rms[totnbOpts + 2];
      m_rmsresult[totnbrms + 5] += (rms[totnbOpts + 5]*rms[totnbOpts + 5]);
      
      // Calculation of vprim (v')
      rms[totnbOpts + 6] = v - rms[totnbOpts + 3];
      m_rmsresult[totnbrms + 6] += (rms[totnbOpts + 6]*rms[totnbOpts + 6]);
      
      // Calculation of Utotprim (Utot')
      //rms[7] = std::sqrt(0.5*(rms[5]*rms[5] + rms[6]*rms[6]));
      rms[totnbOpts + 7] = Utot - rms[totnbOpts + 4];
      m_rmsresult[totnbrms + 7] += (rms[totnbOpts + 7]*rms[totnbOpts + 7]);
      
      // Calculation of urms 
      rms[totnbOpts + 8] = std::sqrt(m_rmsresult[totnbrms + 5]/m_nbStep);
      
      // Calculation of vrms 
      rms[totnbOpts + 9] = std::sqrt(m_rmsresult[totnbrms + 6]/m_nbStep);
      
      // Calculation of Utotrms 
      rms[totnbOpts + 10] = std::sqrt(m_rmsresult[totnbrms + 7]/m_nbStep);
      
      if (m_comp)
      {
      // Turbulence Intensities
      
      // In x-direction
      rms[totnbOpts + 11] = rms[totnbOpts + 8]/rms[totnbOpts + 2];
      
      // In y-direction
      rms[totnbOpts + 12] = rms[totnbOpts + 9]/rms[totnbOpts + 3];
      
      // Total
      rms[totnbOpts + 13] = rms[totnbOpts + 10]/rms[totnbOpts + 4];
      
      // Turbulent kinetic energy
      rms[totnbOpts + 14] = 0.5*(rms[totnbOpts + 8]*rms[totnbOpts + 8]+rms[totnbOpts + 9]*rms[totnbOpts + 9]);
      }
    }
    
    else if (dim == 3) {
      
      //const CFuint nbOpts = 19;//it is the number of the turbulent variables that will be saved
      //const CFuint nbrms = 10; // it is the number of the auxiliary variables used for the calculation of the turbulent variables

      const CFuint totnbOpts = iState*nbOpts;
      const CFuint totnbrms = iState*nbrms;
      
      // Calculation of rho mean
      m_rmsresult[totnbrms] += rho;
      rms[totnbOpts] = m_rmsresult[totnbrms]/m_nbStep;
      
      // Calculation of p mean
      m_rmsresult[totnbrms + 1] += p;
      rms[totnbOpts + 1] = m_rmsresult[totnbrms + 1]/m_nbStep;
      
      // Calculation of u mean
      m_rmsresult[totnbrms + 2] += u;
      rms[totnbOpts + 2] = m_rmsresult[totnbrms + 2]/m_nbStep;
      
      // Calculation of v mean
      m_rmsresult[totnbrms + 3] += v;
      rms[totnbOpts + 3] = m_rmsresult[totnbrms + 3]/m_nbStep;
      
      // Calculation of w mean
      m_rmsresult[totnbrms + 4] += w;
      rms[totnbOpts + 4] = m_rmsresult[totnbrms + 4]/m_nbStep;
      
      // Calculation of Utot mean
      m_rmsresult[totnbrms + 5] += Utot;
      rms[totnbOpts + 5] = m_rmsresult[totnbrms + 5]/m_nbStep;
      
      // Calculation of uprim (u')
      rms[totnbOpts + 6] = u - rms[totnbOpts + 2];
      m_rmsresult[totnbrms + 6] += (rms[totnbOpts + 6]*rms[totnbOpts + 6]);
      
      // Calculation of vprim (v')
      rms[totnbOpts + 7] = v - rms[totnbOpts + 3];
      m_rmsresult[totnbrms + 7] += (rms[totnbOpts + 7]*rms[totnbOpts + 7]);
      
      // Calculation of vprim (w')
      rms[totnbOpts + 8] = w - rms[totnbOpts + 4];
      m_rmsresult[totnbrms + 8] += (rms[totnbOpts + 8]*rms[totnbOpts + 8]);
      
      // Calculation of Utotprim (Utot')
      //rms[7] = std::sqrt(0.5*(rms[5]*rms[5] + rms[6]*rms[6]));
      rms[totnbOpts + 9] = Utot - rms[totnbOpts + 5];
      m_rmsresult[totnbrms + 9] += (rms[totnbOpts + 9]*rms[totnbOpts + 9]);
      
      // Calculation of urms 
      rms[totnbOpts + 10] = std::sqrt(m_rmsresult[totnbrms + 6]/m_nbStep);
      
      // Calculation of vrms 
      rms[totnbOpts + 11] = std::sqrt(m_rmsresult[totnbrms + 7]/m_nbStep);
      
      // Calculation of wrms 
      rms[totnbOpts + 12] = std::sqrt(m_rmsresult[totnbrms + 8]/m_nbStep);
      
      // Calculation of Utotrms 
      rms[totnbOpts + 13] = std::sqrt(m_rmsresult[totnbrms + 9]/m_nbStep);
      
      if (m_comp)
      {
      // Turbulence Intensities
      
      // In x-direction
      rms[totnbOpts + 14] = rms[totnbOpts + 10]/rms[totnbOpts + 2];
      
      // In y-direction
      rms[totnbOpts + 15] = rms[totnbOpts + 11]/rms[totnbOpts + 3];
      
      // In z-direction
      rms[totnbOpts + 16] = rms[totnbOpts + 12]/rms[totnbOpts + 4];
      
      // Total
      rms[totnbOpts + 17] = rms[totnbOpts + 13]/rms[totnbOpts + 5];
      
      // Turbulent kinetic energy
      rms[totnbOpts + 18] = 0.5*(rms[totnbOpts + 10]*rms[totnbOpts + 10]+
				 rms[totnbOpts + 11]*rms[totnbOpts + 11]+
				 rms[totnbOpts + 12]*rms[totnbOpts + 12]);
      }
    }
  }
}

 
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
  DataHandle<CFreal> rms = socket_rms.getDataHandle();
  rms.resize(0);
  m_rmsresult.resize(0);
  
  DataProcessingCom::unsetup(); // last call setup of parent class
}

//////////////////////////////////////////////////////////////////////////////

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
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NavierStokes

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
