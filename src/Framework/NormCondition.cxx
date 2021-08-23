// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/NormCondition.hh"
#include "Framework/SubSystemStatus.hh"
#include "Environment/ObjectProvider.hh"
#include "Framework/Framework.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<NormCondition, StopCondition, FrameworkLib, 1>
normCondProvider("Norm");

//////////////////////////////////////////////////////////////////////////////

void NormCondition::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< CFreal, Config::DynamicOption<> >("valueNorm","Target value of the residual norm.");
}

//////////////////////////////////////////////////////////////////////////////

NormCondition::NormCondition(const std::string& name) :
  StopCondition(name)
{
   addConfigOptionsTo(this);
  _valueNorm = -8.0;
   setParameter("valueNorm",&_valueNorm);
}

//////////////////////////////////////////////////////////////////////////////

NormCondition::~NormCondition()
{
}

//////////////////////////////////////////////////////////////////////////////

bool NormCondition::IsGlobal () const
{
  return false;
}

//////////////////////////////////////////////////////////////////////////////

bool NormCondition::isAchieved (const ConvergenceStatus& status)
{
  if (SubSystemStatusStack::getActive()->getStopSimulation()) return true;
  static bool isFirstIter = true;
  if (isFirstIter) {
    isFirstIter = false;
    return false;
  }
  CFLog(VERBOSE, "NormCondition::isAchieved() => " << status.res << " < " << _valueNorm << "\n");
  return (status.res < _valueNorm);
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
