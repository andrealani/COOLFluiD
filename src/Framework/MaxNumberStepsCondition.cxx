// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/MaxNumberStepsCondition.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/Framework.hh"
#include "Framework/NamespaceSwitcher.hh"
#include "Environment/ObjectProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;

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
  options.addConfigOption< vector<CFuint> >("nbSteps","Maximum number of steps to compute.");
}

//////////////////////////////////////////////////////////////////////////////

MaxNumberStepsCondition::MaxNumberStepsCondition(const std::string& name)
  : StopCondition(name)
{
  addConfigOptionsTo(this);
  _maxNbSteps = vector<CFuint>();
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
  if (SubSystemStatusStack::getActive()->getStopSimulation()) return true;

  NamespaceSwitcher& nsw = 
    NamespaceSwitcher::getInstance(SubSystemStatusStack::getCurrentName());
  const CFuint nbNsp = nsw.getAllNamespaces().size();
  
  if (_maxNbSteps.size() < nbNsp) {
    if (_maxNbSteps.size() == 0) {
      _maxNbSteps.resize(nbNsp, 100);
    }
    else if (_maxNbSteps.size() == 1) {
      const CFuint value = _maxNbSteps[0];
      _maxNbSteps.resize(nbNsp, value);
    }
    else {
      CFLog(ERROR, "MaxNumberStepsCondition::isAchieved() needs one entry per Namespace for option \"nbSteps\"\n");
      abort();
    }
  }
  cf_assert(_maxNbSteps.size() == nbNsp);
  
  const CFuint namespaceID = nsw.getID(true);
  cf_assert(namespaceID < _maxNbSteps.size());
  
  return status.iter >= _maxNbSteps[namespaceID];
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
