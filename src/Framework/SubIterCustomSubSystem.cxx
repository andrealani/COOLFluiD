// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Environment/ObjectProvider.hh"

#include "Framework/StopConditionController.hh"
#include "Framework/MeshCreator.hh"
#include "Framework/MeshAdapterMethod.hh"
#include "Framework/ErrorEstimatorMethod.hh"
#include "Framework/CouplerMethod.hh"
#include "Framework/ConvergenceMethod.hh"
#include "Framework/SpaceMethod.hh"
#include "Framework/DataProcessingMethod.hh"
#include "Framework/OutputFormatter.hh"
#include "Framework/LinearSystemSolver.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/SimulationStatus.hh"
#include "Framework/InteractiveParamReader.hh"
#include "Framework/Framework.hh"
#include "Framework/NamespaceSwitcher.hh"
#include "Framework/SubIterCustomSubSystem.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<SubIterCustomSubSystem,
               SubSystem,
               FrameworkLib,
               1>
subIterCustomSubSystemProvider("SubIterCustomSubSystem");

//////////////////////////////////////////////////////////////////////////////

void SubIterCustomSubSystem::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< std::vector<std::string> >("RunSequence","Sequence for running.");
  options.addConfigOption< CFuint >("NbSubIterations","Number of SubIterations.");
  options.addConfigOption< std::string >("SubIterStopCondition","The stop condition to control the subIteration procedure.");
  options.addConfigOption< std::string >("SubIterStopConditionName","The name to assign to the stop condition to control the subIteration procedure.");

}

//////////////////////////////////////////////////////////////////////////////

SubIterCustomSubSystem::SubIterCustomSubSystem(const std::string& name)
  : StandardSubSystem(name),
  _endSubIteration(false)
{

  addConfigOptionsTo(this);

  m_runSequenceStr = std::vector<std::string>();
  setParameter("RunSequence",&m_runSequenceStr);

  m_nbSubIterations = 1;
  setParameter("NbSubIterations",&m_nbSubIterations);

  m_subIterStopConditionStr = "MaxNumberSteps";
  setParameter("SubIterStopCondition",&m_subIterStopConditionStr);

  m_subIterStopConditionNameStr = "SubIterMaxNumberSteps";
  setParameter("SubIterStopConditionName",&m_subIterStopConditionNameStr);

}

//////////////////////////////////////////////////////////////////////////////

SubIterCustomSubSystem::~SubIterCustomSubSystem()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void SubIterCustomSubSystem::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

  StandardSubSystem::configure(args);

  // builds a stop condition
  Common::SafePtr<StopCondition::PROVIDER> stopCondProv =
    Environment::Factory<StopCondition>::getInstance().getProvider(m_subIterStopConditionStr);
  SelfRegistPtr<StopCondition> subIterStopCondition =
    stopCondProv->create(m_subIterStopConditionNameStr);

  configureNested ( subIterStopCondition.getPtr(), args );

  // Give the stopcondition to the StopConditionController
  // StopConditionController::Create makes a new StopConditionController
  // object, adapted to the current PE mode
  m_subIterStopCondControler.reset(StopConditionController::Create(subIterStopCondition));

}

//////////////////////////////////////////////////////////////////////////////

void SubIterCustomSubSystem::runSequenceID(const CFuint ID)
{

  CFAUTOTRACE;

// std::cout << "Sequence: " << m_runSequenceStr[ID] << std::endl;
  CFchar separator = QualifiedName::separator();
  const CFuint size = StringOps::getWords(m_runSequenceStr[ID],separator).size();
  cf_assert(size > 1);

  if(size == 2)
  {
    const std::string methodName = StringOps::getWords(m_runSequenceStr[ID],separator)[0];
    const std::string functionName = StringOps::getWords(m_runSequenceStr[ID],separator)[1];
    MethodRegistry::getInstance().getMethod(methodName)->run_function(functionName);
  }

  // Loop over one method function
  if(size == 3)
  {
    const std::string methodName = StringOps::getWords(m_runSequenceStr[ID],separator)[0];
    const std::string functionName = StringOps::getWords(m_runSequenceStr[ID],separator)[1];
    const CFuint functionIter = Common::StringOps::from_str<CFuint>(StringOps::getWords(m_runSequenceStr[ID],separator)[2]);
    //for subitrating you can only do one iter
    cf_assert(functionIter == 1);
    Common::SafePtr<Framework::Method> method = MethodRegistry::getInstance().getMethod(methodName);
    for(CFuint iRun = 0; iRun < functionIter; iRun++)
    {
      method->run_function(functionName);
    }
  }

  // Loop over one method function
  if(size > 3)
  {
    ///@todo cf_assert(size is uneven);
    const CFuint functionIter = StringOps::from_str<CFuint>(StringOps::getWords(m_runSequenceStr[ID],separator)[size-1]);
    //for subitrating you can only do one iter
    cf_assert(functionIter == 1);

    const CFuint nbMethods = (size - 1)/2;
    std::vector<std::string> methodName(nbMethods);
    std::vector<std::string> functionName(nbMethods);

    for(CFuint iMtd = 0; iMtd < nbMethods; iMtd++)
    {
      methodName[iMtd] = StringOps::getWords(m_runSequenceStr[ID],separator)[(iMtd*2)];
      functionName[iMtd] = StringOps::getWords(m_runSequenceStr[ID],separator)[(iMtd*2) + 1];
    }

    for(CFuint iRun = 0; iRun < functionIter; iRun++)
    {
      for(CFuint iMtd = 0; iMtd < nbMethods; iMtd++)
      {
        Common::SafePtr<Framework::Method> method =
          MethodRegistry::getInstance().getMethod(methodName[iMtd]);

        method->run_function(functionName[iMtd]);
      }
    }
  }

}

//////////////////////////////////////////////////////////////////////////////

void SubIterCustomSubSystem::run()
{
  CFAUTOTRACE;

  cf_assert(isConfigured());

  vector <Common::SafePtr<SubSystemStatus> > subSysStatusVec =
    SubSystemStatusStack::getInstance().getAllEntries();

  for(CFuint i = 0; i<subSysStatusVec.size();i++){
    // reset to 0 (or Initial Value) the number of iterations and the residual
    subSysStatusVec[i]->setNbIter(m_initialIter);
    // reset to 0 the number of iterations and the residual
    subSysStatusVec[i]->resetResidual();
    //here we say that we will be using subiterations
    subSysStatusVec[i]->setSubIterationFlag(true);

    subSysStatusVec[i]->startWatch();
  }

  m_duration.set(0.);
  Stopwatch<WallTime> stopTimer;
  stopTimer.start();

  // Transfer of the data for the subsystems coupling
  // here, write the data to the other subsystems
  m_couplerMethod.apply(mem_fun<void,CouplerMethod>
                      (&CouplerMethod::dataTransferWrite));

  for ( ; !m_stopCondControler->isAchieved(SubSystemStatusStack::getActive()->getConvergenceStatus()); )
  {
    _endSubIteration = false;
    for(CFuint iSub=0;!_endSubIteration;++iSub)
    {

      /// @todo all this should go into the ConvergenceStatus
      for(CFuint i = 0; i<NamespaceSwitcher::getInstance().getAllNamespaces().size();i++)
      {
        const std::string nspName = NamespaceSwitcher::getInstance().getAllNamespaces()[i]->getName();
        NamespaceSwitcher::getInstance().pushNamespace(nspName);
        Common::SafePtr<SubSystemStatus> subSysStatus = SubSystemStatusStack::getActive();

        subSysStatus->setSubIter(iSub);
        subSysStatus->setIsSubIterationLastStep(false);

        //if stopCond is achieved in any namespace...
        if(m_subIterStopCondControler->isAchieved(SubSystemStatusStack::getActive()->getConvergenceStatus())){
          subSysStatus->setIsSubIterationLastStep(true);
          _endSubIteration = true;
        }

        NamespaceSwitcher::getInstance().popNamespace();
      }

      if(_endSubIteration == true){
        for(CFuint i = 0; i<NamespaceSwitcher::getInstance().getAllNamespaces().size();i++)
        {
          const std::string nspName = NamespaceSwitcher::getInstance().getAllNamespaces()[i]->getName();
          NamespaceSwitcher::getInstance().pushNamespace(nspName);
          Common::SafePtr<SubSystemStatus> subSysStatus = SubSystemStatusStack::getActive();
          subSysStatus->setIsSubIterationLastStep(true);
          NamespaceSwitcher::getInstance().popNamespace();
        }
      }

      // read the interactive parameters
      getInteractiveParamReader()->readFile();

      for(CFuint i=0;i < m_runSequenceStr.size();i++){
        runSequenceID(i);
      }

      writeConvergenceOnScreen();

      // write solution to file
      writeSolution(false);

    } //end for the subiteration loop

  } // end for convergence loop

    // estimate the error one final time
    m_errorEstimatorMethod.apply(mem_fun<void,ErrorEstimatorMethod>(&ErrorEstimatorMethod::estimate));

  for(CFuint i = 0; i<subSysStatusVec.size();i++){
    subSysStatusVec[i]->stopWatch();
  }

  stopTimer.stop ();

  CFLog(NOTICE, "SubSystem WallTime: " << stopTimer << "s\n");
  m_duration = subSysStatusVec[0]->readWatchHMS();
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
