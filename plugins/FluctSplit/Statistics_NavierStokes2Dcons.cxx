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

#include "FluctSplit/Statistics_NavierStokes2Dcons.hh"
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

MethodCommandProvider<Statistics_NavierStokes2Dcons,
                      DataProcessingData,
		      FluctSplitSpaceTimeModule>
rmsns2dconsProvider("Statistics_NavierStokes2Dcons");

//////////////////////////////////////////////////////////////////////////////

void Statistics_NavierStokes2Dcons::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< std::string >("OutputFile","Name of Output File.");
  options.addConfigOption< CFuint >("SaveRate","Rate for saving the output file with.");
  options.addConfigOption< bool >("Restart","Continue averaging.");
}

//////////////////////////////////////////////////////////////////////////////

Statistics_NavierStokes2Dcons::Statistics_NavierStokes2Dcons( const std::string& name) :
  DataProcessingCom(name),
  m_fhandle(),
  socket_nodes("nodes"),
  socket_states("states"),
  socket_meanflowNS("meanflowNS")
{
  addConfigOptionsTo(this);

  m_OutputFile = "Probe.dat";
  setParameter("OutputFile",&m_OutputFile);
  
  m_saveRate = 1;
  setParameter("SaveRate",&m_saveRate);
  
  m_Restart = false;
  setParameter("Restart",&m_Restart);
}

//////////////////////////////////////////////////////////////////////////////

Statistics_NavierStokes2Dcons::~Statistics_NavierStokes2Dcons()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
Statistics_NavierStokes2Dcons::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;
  result.push_back(&socket_nodes);
  result.push_back(&socket_states);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSource> >
Statistics_NavierStokes2Dcons::providesSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSource> > result;
  result.push_back(&socket_meanflowNS);
  return result;
}
//////////////////////////////////////////////////////////////////////////////

void Statistics_NavierStokes2Dcons::setup()
{
  CFAUTOTRACE;

  DataProcessingCom::setup(); // first call setup of parent class
  
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  const CFuint nbstates = states.size();

  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  
  DataHandle<RealVector> meanflow = socket_meanflowNS.getDataHandle();
  meanflow.resize(nbstates);

  for (CFuint i = 0 ; i < nbstates; i++){
    meanflow[i].resize(nbEqs);
   for (CFuint j = 0 ; j < nbEqs; j++) {
     meanflow[i][j] = 0.0;
   }
  }
  
  m_varSet->setup();

  m_rhobar.resize(nbstates) ;
  m_Ubar.resize(nbstates) ;
  m_Vbar.resize(nbstates) ;
  m_UVbar.resize(nbstates) ;
  m_UUbar.resize(nbstates) ;
  m_VVbar.resize(nbstates) ;
  m_pbar.resize(nbstates) ;
  m_rhoEbar.resize(nbstates) ;
  m_rms.resize(nbstates) ;
  m_devRho.resize(nbstates) ;

  
  for (CFuint i = 0 ; i < nbstates; i++){
    m_rhobar[i] = 0.0 ;
    m_Ubar[i] = 0.0 ;
    m_Vbar[i] = 0.0 ;
    m_UVbar[i] = 0.0 ;
    m_UUbar[i] = 0.0 ;
    m_VVbar[i] = 0.0 ;
    m_pbar[i] = 0.0 ;
    m_rhoEbar[i] = 0.0 ;
    m_nbStep = 0 ;
    m_rms[i] = 0.0 ;
    m_devRho[i] = 0.0 ;
  }
  
  for (CFuint i = 0 ; i < nbstates; i++){
    mean_rhobar[i] = 0.0 ;
    mean_Ubar[i] = 0.0 ;
    mean_Vbar[i] = 0.0 ;
    mean_UVbar[i] = 0.0 ;
    mean_UUbar[i] = 0.0 ;
    mean_VVbar[i] = 0.0 ;
    mean_pbar[i] = 0.0 ;
    mean_rhoEbar[i] = 0.0 ;
    mean_rms[i] = 0.0 ;
    mean_devRho[i] = 0.0 ;
    mean_nbStep = 0;
  }



//   if(m_Restart) {
//   
//     m_fhandle = Environment::SingleBehaviorFactory<Environment::FileHandlerOutput>::getInstance().create();
// 
//     cf_assert (!m_fhandle->isopen());
// 
//     boost::filesystem::path m_path_name = PathAppender::getInstance().appendParallel(m_OutputFile);
//     ifstream& myfile = m_fhandle->open(m_path_name);  
// 
//     std::string dummy;
//     CFreal data;
//  
//     for (CFuint i=0; i<4; i++)
//       getline (myfile,dummy);
// 
//     for (CFuint i=0; i<nbstates; i++) {
//       myfile >> data;
//       myfile >> data;
//       myfile >> mean_rhobar[i];
//       myfile >> mean_Ubar[i];
//       myfile >> mean_Vbar[i];
//       myfile >>  mean_UVbar[i];
//       myfile >> mean_UUbar[i];
//       myfile >> mean_VVbar[i];
//       myfile >> mean_pbar[i];
//       myfile >> mean_rhoEbar[i];
//       myfile >>  mean_rms[i];
//       myfile >> mean_devRho[i];
//       myfile >> mean_nbStep;
//     }
//     m_fhandle->close();
//   }
  
}

//////////////////////////////////////////////////////////////////////////////

void Statistics_NavierStokes2Dcons::execute()
{
  
  
  Common::SafePtr<SubSystemStatus> subSysStatus = SubSystemStatusStack::getActive();
  const CFuint iter = subSysStatus->getNbIter();
  
  m_nbStep = iter;
  
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle<RealVector> meanflow = socket_meanflowNS.getDataHandle();
  const CFuint nbstates = states.size();

  m_fhandle = Environment::SingleBehaviorFactory<Environment::FileHandlerOutput>::getInstance().create();

  boost::filesystem::path m_path_name = PathAppender::getInstance().appendParallel(m_OutputFile);
  ofstream& myfile = m_fhandle->open(m_path_name);

  myfile << "TITLE = \"RMS\"" << "\n";
  myfile << "VARIABLES = \"x0\" \"x1\" \"Rhobar\" \"Ubar\" \"Vbar\" \"UVbar\" \"UUbar\" \"VVbar\" \"Pbar\" \"rhoEbar\" \"Prms\" \"Rhorms\" \"nbStep\"" << "\n";
  myfile << "ZONE T = " << " ";
  myfile << m_OutputFile << "\n";
  myfile << "DATAPACKING = POINT" << "\n";

//   SAVE RATE !!!
  
 
  for (CFuint i = 0 ; i < nbstates; i++){
    const CFreal x_p = (states[i]->getCoordinates())[XX];
    const CFreal y_p = (states[i]->getCoordinates())[YY];
    const CFreal rho = (*states[i])[0];
    const CFreal rhou = (*states[i])[1];
    const CFreal rhov = (*states[i])[2];
    const CFreal rhoE = (*states[i])[3];
    CFreal gamma = m_varSet->getModel()->getGamma();
    CFreal u = rhou/rho;
    CFreal v = rhov/rho;
    CFreal p = (gamma-1.0)*(rhoE-0.5*rho*(u*u+v*v)) ;
    m_rhobar[i] += rho;
    m_Ubar[i] += u ;
    m_Vbar[i] += v ;
    m_UVbar[i] += u*v ;
    m_UUbar[i] += u*u  ;
    m_VVbar[i] += v*v ;
    m_pbar[i] += p ;
    m_rhoEbar[i] += rhoE;
    m_rms[i] += (p - m_pbar[i])*(p - m_pbar[i]);
    m_devRho[i] += (rho - m_rhobar[i])*(rho - m_rhobar[i]);


    myfile
      << x_p << " "
      << y_p << " "
      << m_rhobar[i]/m_nbStep << " "
      << m_Ubar[i]/m_nbStep   << " "
      << m_Vbar[i]/m_nbStep << " "
      <<  m_UVbar[i]/m_nbStep << " "
      <<  m_UUbar[i]/m_nbStep << " "
      <<  m_VVbar[i]/m_nbStep << " "
      <<  m_pbar[i]/m_nbStep << " "
      <<  m_rhoEbar[i]/m_nbStep << " "
      <<  sqrt(m_rms[i])/m_nbStep << " "
      <<  sqrt(m_devRho[i])/m_nbStep << " "
      << m_nbStep    << "\n ";



    meanflow[i][0] =  m_rhobar[i]/m_nbStep;
    meanflow[i][1] =  m_Ubar[i]/m_nbStep;  
    meanflow[i][2] =  m_Vbar[i]/m_nbStep;
    meanflow[i][3] =  m_rhoEbar[i]/m_nbStep;
  }

  m_fhandle->close();

}

//////////////////////////////////////////////////////////////////////////////

void Statistics_NavierStokes2Dcons::unsetup()
{
  
  DataHandle<RealVector> meanflow = socket_meanflowNS.getDataHandle();
  meanflow.resize(0);
  
  DataProcessingCom::unsetup(); // at last call setup of parent class
}

//////////////////////////////////////////////////////////////////////////////

void Statistics_NavierStokes2Dcons::configure ( Config::ConfigArgs& args )
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

