// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.



#include "Framework/MaxNumberSubIterCondition.hh"
#include "Framework/SubSystemStatus.hh"
#include "Environment/ObjectProvider.hh"
#include "Framework/Framework.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<MaxNumberSubIterCondition,
	       StopCondition,
               FrameworkLib,
	       1>
maxNbSubIterConditionProvider("MaxNumberSubIter");

//////////////////////////////////////////////////////////////////////////////

void MaxNumberSubIterCondition::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< CFuint >("nbSteps","Maximum number of steps to compute.");
}

//////////////////////////////////////////////////////////////////////////////

MaxNumberSubIterCondition::MaxNumberSubIterCondition(const std::string& name)
 : StopCondition(name)
{
  addConfigOptionsTo(this);
  _maxNbSteps = 1;
  setParameter("nbSteps",&_maxNbSteps);
}

//////////////////////////////////////////////////////////////////////////////

MaxNumberSubIterCondition::~MaxNumberSubIterCondition()
{
}

//////////////////////////////////////////////////////////////////////////////

bool MaxNumberSubIterCondition::IsGlobal () const
{
  return false;
}

//////////////////////////////////////////////////////////////////////////////

bool MaxNumberSubIterCondition::isAchieved (const ConvergenceStatus& status)
{
  if (SubSystemStatusStack::getActive()->getStopSimulation()) return true;
  return status.subiter >= _maxNbSteps;
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
