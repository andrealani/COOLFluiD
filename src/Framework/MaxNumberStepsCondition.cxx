// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.



#include "Framework/MaxNumberStepsCondition.hh"
#include "Framework/SubSystemStatus.hh"
#include "Environment/ObjectProvider.hh"
#include "Framework/Framework.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<MaxNumberStepsCondition,
	       StopCondition,
               FrameworkLib,
	       1>
maxNbStepsConditionProvider("MaxNumberSteps");

//////////////////////////////////////////////////////////////////////////////

void MaxNumberStepsCondition::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< CFuint >("nbSteps","Maximum number of steps to compute.");
}

//////////////////////////////////////////////////////////////////////////////

MaxNumberStepsCondition::MaxNumberStepsCondition(const std::string& name)
 : StopCondition(name)
{
  addConfigOptionsTo(this);
  _maxNbSteps = 100;
  setParameter("nbSteps",&_maxNbSteps);
}

//////////////////////////////////////////////////////////////////////////////

MaxNumberStepsCondition::~MaxNumberStepsCondition()
{
}

//////////////////////////////////////////////////////////////////////////////

bool MaxNumberStepsCondition::IsGlobal () const
{
  return false;
}

//////////////////////////////////////////////////////////////////////////////

bool MaxNumberStepsCondition::isAchieved (const ConvergenceStatus& status)
{
  return status.iter >= _maxNbSteps;
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
