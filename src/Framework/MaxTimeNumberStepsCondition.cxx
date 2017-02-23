// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/MaxTimeNumberStepsCondition.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/Framework.hh"
#include "Framework/NamespaceSwitcher.hh"
#include "Environment/ObjectProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<MaxTimeNumberStepsCondition,
			    StopCondition,
			    FrameworkLib,
			    1>
maxTimeNbStepsConditionProvider("MaxTimeNumberSteps");

//////////////////////////////////////////////////////////////////////////////

void MaxTimeNumberStepsCondition::defineConfigOptions(Config::OptionList& options)
{
 options.addConfigOption< CFreal, Config::DynamicOption<> >
   ("maxTime","Maximum Physical Time for the computation.");
}

//////////////////////////////////////////////////////////////////////////////

MaxTimeNumberStepsCondition::MaxTimeNumberStepsCondition(const std::string& name)
  : MaxNumberStepsCondition(name)
{
  addConfigOptionsTo(this);
  
  m_maxTime = 1.;
  setParameter("maxTime",&m_maxTime);
}

//////////////////////////////////////////////////////////////////////////////

MaxTimeNumberStepsCondition::~MaxTimeNumberStepsCondition()
{
}

//////////////////////////////////////////////////////////////////////////////

bool MaxTimeNumberStepsCondition::IsGlobal () const
{
  return false;
}

//////////////////////////////////////////////////////////////////////////////

bool MaxTimeNumberStepsCondition::isAchieved (const ConvergenceStatus& status)
{
  if (SubSystemStatusStack::getActive()->getStopSimulation()) return true;
  if (MaxNumberStepsCondition::isAchieved(status)) return true;
  
  // set the maximum time
  SubSystemStatusStack::getActive()->setMaxTimeDim(m_maxTime);
  return status.time >= m_maxTime;
}
    
//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
