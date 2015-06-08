// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include <sstream>


#include "Common/PE.hh"
#include "Common/NullPointerException.hh"
#include "Common/CFLog.hh"
#include "Common/ProcessInfo.hh"

#include "Environment/FileHandlerOutput.hh"
#include "Environment/DirPaths.hh"
#include "Environment/ObjectProvider.hh"
#include "Environment/SingleBehaviorFactory.hh"

#include "Framework/PrePostProcessingSubSystem.hh"
#include "Framework/MeshData.hh"
#include "Framework/CFL.hh"
#include "Framework/MeshCreator.hh"
#include "Framework/MeshAdapterMethod.hh"
#include "Framework/ErrorEstimatorMethod.hh"
#include "Framework/CouplerMethod.hh"
#include "Framework/ConvergenceMethod.hh"
#include "Framework/SpaceMethod.hh"
#include "Framework/DataProcessingMethod.hh"
#include "Framework/OutputFormatter.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/LinearSystemSolver.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/StopConditionController.hh"
#include "Framework/InteractiveParamReader.hh"
#include "Framework/CommandsToTRSMapper.hh"
#include "Framework/PathAppender.hh"
#include "Framework/MeshCreator.hh"
#include "Framework/NamespaceSwitcher.hh"
#include "Framework/Namespace.hh"
#include "Framework/Framework.hh"
#include "Framework/SimulationStatus.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace boost::filesystem;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<PrePostProcessingSubSystem,
               SubSystem,
               FrameworkLib,
               1>
aPrePostProcessingSubSysProvider("PrePostProcessingSubSystem");

//////////////////////////////////////////////////////////////////////////////

PrePostProcessingSubSystem::PrePostProcessingSubSystem(const std::string& name)
  : StandardSubSystem(name)
{
}

//////////////////////////////////////////////////////////////////////////////

PrePostProcessingSubSystem::~PrePostProcessingSubSystem()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void PrePostProcessingSubSystem::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

  StandardSubSystem::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void PrePostProcessingSubSystem::setup()
{
  CFAUTOTRACE;
  StandardSubSystem::setup();
}

//////////////////////////////////////////////////////////////////////////////

void PrePostProcessingSubSystem::run()
{
  CFAUTOTRACE;
  cf_assert(isConfigured());

  m_forcedStop = false;

  std::vector <Common::SafePtr<SubSystemStatus> > subSysStatusVec = SubSystemStatusStack::getInstance().getAllEntries();

  ///@todo this should not be here!!
  for(CFuint i = 0; i<subSysStatusVec.size();i++)
  {
    // reset to 0 (or Initial Value) the number of iterations and the residual
    subSysStatusVec[i]->setNbIter(m_initialIter);
    subSysStatusVec[i]->setCurrentTimeDim(m_initialTime);
    // reset to 0 the number of iterations and the residual
    subSysStatusVec[i]->resetResidual();
    subSysStatusVec[i]->startWatch();
  }

  m_duration.set(0.);
  Stopwatch<WallTime> stopTimer;
  stopTimer.start();
  
  // read the interactive parameters
  getInteractiveParamReader()->readFile();

  // pre-process the data
  m_dataPreProcessing.apply(mem_fun<void,DataProcessingMethod>
      (&DataProcessingMethod::processData));

  // post-process the data
  m_dataPostProcessing.apply(mem_fun<void,DataProcessingMethod>
            (&DataProcessingMethod::processData));

  writeConvergenceOnScreen();

  // write solution to file
  bool dontforce = false;
  writeSolution(dontforce);

  for(CFuint i = 0; i < subSysStatusVec.size() ; i++)
  {
    subSysStatusVec[i]->stopWatch();
  }
  stopTimer.stop ();

  CFLog(NOTICE, "SubSystem WallTime: " << stopTimer << "s\n");
  m_duration = subSysStatusVec[0]->readWatchHMS();

}

//////////////////////////////////////////////////////////////////////////////

void PrePostProcessingSubSystem::unsetup()
{
  CFAUTOTRACE;
  StandardSubSystem::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

vector<Method*> PrePostProcessingSubSystem::getMethodList()
{
  vector<Method*>  mList = StandardSubSystem::getMethodList();
  return mList;
}

//////////////////////////////////////////////////////////////////////////////

void PrePostProcessingSubSystem::setCollaborators()
{
  CFAUTOTRACE;
  StandardSubSystem::setCollaborators();
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
