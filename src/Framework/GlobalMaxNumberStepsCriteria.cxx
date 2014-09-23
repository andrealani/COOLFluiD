// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Environment/ObjectProvider.hh"

#include "Framework/GlobalMaxNumberStepsCriteria.hh"
#include "Framework/SimulationStatus.hh"
#include "Framework/Framework.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<GlobalMaxNumberStepsCriteria, GlobalStopCriteria,FrameworkLib, 1>
globalMaxNbStepsCriteriaProvider("GlobalMaxNumberSteps");

//////////////////////////////////////////////////////////////////////////////

void GlobalMaxNumberStepsCriteria::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< CFuint >("nbSteps","Maximum number of overall steps to perform.");
}

//////////////////////////////////////////////////////////////////////////////

GlobalMaxNumberStepsCriteria::GlobalMaxNumberStepsCriteria(const std::string& name)
 : GlobalStopCriteria(name)
{
  addConfigOptionsTo(this);
  m_maxNbSteps = 1;
  setParameter("nbSteps",&m_maxNbSteps);
}

//////////////////////////////////////////////////////////////////////////////

GlobalMaxNumberStepsCriteria::~GlobalMaxNumberStepsCriteria()
{
}

//////////////////////////////////////////////////////////////////////////////

bool GlobalMaxNumberStepsCriteria::isSatisfied ()
{
  return SimulationStatus::getInstance().getNbIter() >= m_maxNbSteps;
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
