// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/CustomSubSystem.hh"
#include "Framework/NamespaceSwitcher.hh"
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
#include "Framework/InteractiveParamReader.hh"
#include "Framework/Framework.hh"
#include "Environment/ObjectProvider.hh"
#include "Common/PEFunctions.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<CustomSubSystem,
               SubSystem,
               FrameworkLib,
               1>
sequentialMFSubSysProvider("CustomSubSystem");

//////////////////////////////////////////////////////////////////////////////

void CustomSubSystem::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< std::vector<std::string> >("RunSequence","Sequence for running.");
}

//////////////////////////////////////////////////////////////////////////////

CustomSubSystem::CustomSubSystem(const std::string& name)
  : StandardSubSystem(name)
{
  CFAUTOTRACE;
  addConfigOptionsTo(this);

  m_runSequenceStr = std::vector<std::string>();
  setParameter("RunSequence",&m_runSequenceStr);
}

//////////////////////////////////////////////////////////////////////////////

CustomSubSystem::~CustomSubSystem()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void CustomSubSystem::runSequenceID(const CFuint ID)
{
  CFAUTOTRACE;

//std::cout << "Sequence: " << m_runSequenceStr[ID] << std::endl;
  CFchar separator = QualifiedName::separator();
  const CFuint size = StringOps::getWords(m_runSequenceStr[ID],separator).size();
  cf_assert(size > 1);

  if(size == 2)
  {
    const std::string methodName   = StringOps::getWords(m_runSequenceStr[ID],separator)[0];
    const std::string functionName = StringOps::getWords(m_runSequenceStr[ID],separator)[1];
    executeMethodFunction(methodName, functionName);
  }

  // Loop over one method function
  if(size == 3)
  {
    const std::string methodName   = StringOps::getWords(m_runSequenceStr[ID],separator)[0];
    const std::string functionName = StringOps::getWords(m_runSequenceStr[ID],separator)[1];
    const CFuint functionIter   = StringOps::from_str<CFuint>(StringOps::getWords(m_runSequenceStr[ID],separator)[2]);
    for(CFuint iRun = 0; iRun < functionIter; ++iRun)
    {
      executeMethodFunction(methodName, functionName);
    }
  }

  // Loop over one method function
  if(size > 3)
  {
    /// @todo cf_assert(size is uneven);
    const CFuint functionIter = StringOps::from_str<CFuint>(StringOps::getWords(m_runSequenceStr[ID],separator)[size-1]);
    const CFuint nbMethods = (size - 1)/2;
    std::vector<std::string> methodName(nbMethods);
    std::vector<std::string> functionName(nbMethods);

    for(CFuint iMtd = 0; iMtd < nbMethods; iMtd++)
    {
      methodName[iMtd]   = StringOps::getWords(m_runSequenceStr[ID],separator)[(iMtd*2)];
      functionName[iMtd] = StringOps::getWords(m_runSequenceStr[ID],separator)[(iMtd*2) + 1];
    }

    for(CFuint iRun = 0; iRun < functionIter; ++iRun)
    {
      for(CFuint iMtd = 0; iMtd < nbMethods; iMtd++)
      {
        executeMethodFunction(methodName[iMtd],functionName[iMtd]);
      }
    }
  }

}

//////////////////////////////////////////////////////////////////////////////

void CustomSubSystem::executeMethodFunction(const std::string& method_name, const std::string& function_name)
{
  SafePtr<Method> method_ptr = MethodRegistry::getInstance().getMethod(method_name);
  const int rank = PE::GetPE().GetRank("Default");
  if (PE::GetPE().isRankInGroup(rank, method_ptr->getNamespace())) {
    CFLog(INFO, "Executing function [" << function_name << "] in method ["<<  method_name << "]\n");
    method_ptr->run_function(function_name);
    CFLog(VERBOSE, "Finished  function [" << function_name << "] in method ["<<  method_name << "]\n");
  }
}
    
//////////////////////////////////////////////////////////////////////////////

void CustomSubSystem::run()
{
  CFAUTOTRACE;

  cf_assert(isConfigured());

  vector <Common::SafePtr<SubSystemStatus> > subSysStatusVec =
    SubSystemStatusStack::getInstance().getAllEntries();

  for(CFuint i = 0; i<subSysStatusVec.size();i++)
  {
    // reset to 0 (or Initial Value) the number of iterations and the residual
    subSysStatusVec[i]->setNbIter(m_initialIter);
    // reset to 0 the number of iterations and the residual
    subSysStatusVec[i]->resetResidual();
    subSysStatusVec[i]->startWatch();
  }
  
  NamespaceSwitcher& nsw = NamespaceSwitcher::getInstance(SubSystemStatusStack::getCurrentName());
  const string sssName = nsw.getName(mem_fun<string,Namespace>(&Namespace::getSubSystemStatusName), true);
  SafePtr<SubSystemStatus> currSSS = SubSystemStatusStack::getInstance().getEntry(sssName);
  cf_assert(currSSS.isNotNull());
  
  m_duration.set(0.);
  Stopwatch<WallTime> stopTimer;
  stopTimer.start();

  // Transfer of the data for the subsystems coupling
  // here, write the data to the other subsystems
  m_couplerMethod.apply(mem_fun<void,CouplerMethod>
                      (&CouplerMethod::dataTransferWrite));
  
  const string ssGroupName = SubSystemStatusStack::getCurrentName();
  for ( ; (!m_stopCondControler->isAchieved(currSSS->getConvergenceStatus())) && (!m_forcedStop); ) {
    
    // read the interactive parameters
    runSerial<void, InteractiveParamReader, &InteractiveParamReader::readFile>
      (&*getInteractiveParamReader(), ssGroupName);
    
    for (CFuint i=0; i < m_runSequenceStr.size(); ++i) { runSequenceID(i); }
    
    writeConvergenceOnScreen();
    
    // write solution to file
    const bool force_write = false;
    writeSolution (force_write);
    
  } // end for convergence loop
  
  // estimate the error one final time
  m_errorEstimatorMethod.apply(mem_fun<void,ErrorEstimatorMethod>(&ErrorEstimatorMethod::estimate));
  
  for(CFuint i = 0; i < subSysStatusVec.size() ; i++) {
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
