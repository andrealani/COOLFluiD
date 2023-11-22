// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.



#include "Environment/ObjectProvider.hh"
#include "Common/CFLog.hh"

#include "Framework/MaxTimeCondition.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/Framework.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<MaxTimeCondition,
               StopCondition,
               FrameworkLib,
               1>
maxTimeConditionProvider("MaxTime");

//////////////////////////////////////////////////////////////////////////////

void MaxTimeCondition::defineConfigOptions(Config::OptionList& options)
{
 options.addConfigOption< CFreal, Config::DynamicOption<> >("maxTime","Maximum Physical Time for the computation.");
}

//////////////////////////////////////////////////////////////////////////////

MaxTimeCondition::MaxTimeCondition(const std::string& name)
 : StopCondition(name)
{
  addConfigOptionsTo(this);
  
  m_maxTime = -1.; // AL: WARNING: changed default, used to be =1
  setParameter("maxTime",&m_maxTime);
}

//////////////////////////////////////////////////////////////////////////////

MaxTimeCondition::~MaxTimeCondition()
{
}

//////////////////////////////////////////////////////////////////////////////

bool MaxTimeCondition::IsGlobal () const
{
  return false;
}

//////////////////////////////////////////////////////////////////////////////

bool MaxTimeCondition::isAchieved (const ConvergenceStatus& status)
{
  if (SubSystemStatusStack::getActive()->getStopSimulation()) return true;

  CFLog(VERBOSE, "MaxTimeCondition::isAchieved() => status.time = " << status.time << "\n");
  
  Common::SafePtr<SubSystemStatus> subSysStatus = SubSystemStatusStack::getActive();

  // if the user has set a max time use it
  if (m_maxTime >= 0.) {
    subSysStatus->setMaxTimeDim(m_maxTime);
  }
  else { // otherwise use the one set by some other object (e.g. unsteady BC)
    m_maxTime = subSysStatus->getMaxTimeDim();
  }
  cf_assert(m_maxTime >= 0.);
  
  static bool isFirstIter = true;
  if (isFirstIter)
  {
    isFirstIter = false;
    return false;
  }

  bool timeOver = (status.time >= m_maxTime);

  /// @todo explain this or remove it
  if (status.res < -30.0) timeOver = true;
  
  return timeOver;
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
