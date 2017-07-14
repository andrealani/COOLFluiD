
#include <sstream>

#include "Common/PE.hh"
#include "Common/ProcessInfo.hh"
#include "Framework/MeshData.hh"
#include "Common/CFLog.hh"
#include "Framework/CFL.hh"
#include "Framework/MeshCreator.hh"
#include "Framework/MeshAdapterMethod.hh"
#include "Framework/ErrorEstimatorMethod.hh"
#include "Framework/CouplerMethod.hh"
#include "Framework/ConvergenceMethod.hh"
#include "Framework/SpaceMethod.hh"
#include "Framework/DataProcessing.hh"
#include "Framework/OutputFormatter.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/LinearSystemSolver.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/StopConditionController.hh"
#include "Framework/InteractiveParamReader.hh"
#include "Framework/CommandsToTRSMapper.hh"
#include "Common/NullPointerException.hh"
#include "Environment/DirPaths.hh"
#include "Environment/ObjectProvider.hh"
#include "Environment/CFEnvVars.hh"
#include "Framework/MeshCreator.hh"
#include "Common/Stopwatch.hh"

#include "RemeshMeandr/RemeshMeandr.hh"
#include "RemeshMeandr/AdaptSubSystem.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Environment;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace RemeshMeandros {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<AdaptSubSystem, SubSystem, RemeshMeandrModule, 1> adaptSubSysProvider("AdaptSubSystem");

//////////////////////////////////////////////////////////////////////////////

void AdaptSubSystem::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< CFuint >("NumOfEstimate","Number of error estimations during solution used for time dependent adaptation.");
}

//////////////////////////////////////////////////////////////////////////////

AdaptSubSystem::AdaptSubSystem(const std::string& name)
  : StandardSubSystem(name)
 {
   addConfigOptionsTo(this);
  m_estimateNum = 3;
   setParameter("NumOfEstimate",&m_estimateNum);



 }

//////////////////////////////////////////////////////////////////////////////

AdaptSubSystem::~AdaptSubSystem()
{
}

//////////////////////////////////////////////////////////////////////////////

void AdaptSubSystem::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

  StandardSubSystem::configure(args);

  // builds a stop condition
//  StopCondition::PROVIDER* const stopCondProv = Environment::Factory<StopCondition>::getInstance().getProvider( m_stopConditionStr);
//  SelfRegistPtr<StopCondition> stopCondition = stopCondProv->create( stopCondProv->getName());

//  configureNested ( stopCondition.getPtr(), args );

  m_pMaxTimeStopCond = m_stopCondControler->getStopCondition().d_castTo<MaxTimeCondition>();

}


//////////////////////////////////////////////////////////////////////////////

void AdaptSubSystem::run()
{
  CFAUTOTRACE;

  CFout << "###### Processing Phase ######\n";

  const bool isParallel = PE::GetPE().IsParallel ();

  Common::SafePtr<SubSystemStatus> subSysStatus =
    SubSystemStatusStack::getActive();

  vector <Common::SafePtr<SubSystemStatus> > subSysStatusVec =
    SubSystemStatusStack::getInstance().getAllEntries();

  bool bEstimEnd = false;
  CFreal t_estim = 0.;
  CFreal dt_estim = 0.;

  if ( m_pMaxTimeStopCond.isNotNull() )
  {
	CFout << "#####--- MaxTime = " << m_pMaxTimeStopCond->getMaxTime() << "\n";
	CFout << "--- NumEstimate = " << m_estimateNum << "\n";

	if ( m_estimateNum < 2)
	{
      m_estimateNum = 2;
	  t_estim = m_pMaxTimeStopCond->getMaxTime();
	}

	dt_estim = m_pMaxTimeStopCond->getMaxTime() / ( m_estimateNum - 1);
  }
  else
  {
    bEstimEnd = true;
  }

    cf_assert(isConfigured());

    // reset to 0 (or Initial Value) the number of iterations and the residual
    subSysStatus->setNbIter(m_initialIter);
    subSysStatus->startWatch();

    // reset to 0 the number of iterations and the residual
    subSysStatus->resetResidual();

    m_duration.set(0.);

    DataHandle<State*, GLOBAL>* states = 0;

    if (isParallel) {
      states = new DataHandle<State*, GLOBAL>
	(MeshDataStack::getActive()->getDataStorage()->getGlobalData<State*>("states"));
      
      if (CFEnv::getInstance().getVars()->SyncAlgo != "Old") {
	states->synchronize();
      } 
      else {
        states->beginSync();
        states->endSync();
      }
    }

    // set the handle to the global states in the
    // convergence method
    //  m_convergenceMethod[0]->setGlobalStates(states);

    Stopwatch<WallTime> stopTimer;
    stopTimer.start();

    // 	// estimation for initial soultion
    //     m_errorEstimatorMethod.apply(mem_fun<void,ErrorEstimatorMethod>(&ErrorEstimatorMethod::estimate));

    for ( ; !m_stopCondControler->isAchieved(SubSystemStatusStack::getActive()->getConvergenceStatus()); )
    {
	  //CFout << " glob t " << subSysStatus->getCurrentTimeDim() << "\n";

      // read the interactive parameters
      getInteractiveParamReader()->readFile();

      // PreProcess the data
      m_dataPreProcessing.apply(mem_fun<void,DataProcessingMethod>
                               (&DataProcessingMethod::processData));

      // here be careful : in implicit it will be already updated....
      // CFL::getInstance().update();

      // Transfer of the data for the subsystems coupling
      // reads the data from other subsystems
      m_couplerMethod.apply(mem_fun<void,CouplerMethod>
                           (&CouplerMethod::dataTransferRead));


      m_convergenceMethod[0]->takeStep();

      // Transfer of the data for the subsystems coupling
      // here, write the data to the other subsystems
      m_couplerMethod.apply(mem_fun<void,CouplerMethod>
                           (&CouplerMethod::dataTransferWrite));

      writeConvergenceOnScreen();

      // error estimation
      //CFout << "--- t_estim = " << t_estim << "\n";

      if ( subSysStatus->getCurrentTimeDim() >= t_estim && !bEstimEnd )
      {
        m_errorEstimatorMethod.apply(mem_fun<void,ErrorEstimatorMethod>(&ErrorEstimatorMethod::estimate));
        writeSolution(true);
		    t_estim += dt_estim;
      }

      // write solution to file
      writeSolution(false);

      // PostProcess the data
      m_dataPostProcessing.apply(mem_fun<void,DataProcessingMethod>
                            (&DataProcessingMethod::processData));

    }

     writeSolution(true);


    if ( bEstimEnd )
      m_errorEstimatorMethod.apply(mem_fun<void,ErrorEstimatorMethod>(&ErrorEstimatorMethod::estimate));


    // adaptation using glaobal_metric
    m_meshAdapterMethod.apply( mem_fun<void,MeshAdapterMethod>( &MeshAdapterMethod::adaptMesh ) );


    subSysStatus->stopWatch();
    stopTimer.stop ();

    CFLog(NOTICE, "SubSystem WallTime: " << stopTimer << "s\n");
    m_duration = subSysStatus->readWatchHMS();

}

//////////////////////////////////////////////////////////////////////////////

void AdaptSubSystem::registActionListeners()
{
  CFAUTOTRACE;

  // always call the parent class
  StandardSubSystem::registActionListeners();
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace RemeshMeandros

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
