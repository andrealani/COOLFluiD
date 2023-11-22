// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Environment/ObjectProvider.hh"

#include "Framework/AbsoluteNormAndMaxIter.hh"
#include "Framework/SubSystemStatus.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<AbsoluteNormAndMaxIter, StopCondition, FrameworkLib, 1>
aAbsoluteNormAndMaxIterProvider("AbsoluteNormAndMaxIter");

//////////////////////////////////////////////////////////////////////////////

void AbsoluteNormAndMaxIter::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< CFreal, Config::DynamicOption<> >("AbsNorm","Target value of the relative residual norm.");
   options.addConfigOption< CFuint, Config::DynamicOption<> >("MaxIter","Maximum number of iterations to perform.");
}

//////////////////////////////////////////////////////////////////////////////

AbsoluteNormAndMaxIter::AbsoluteNormAndMaxIter(const std::string& name) :
  StopCondition(name)
{
  addConfigOptionsTo(this);

  m_norm = -10.;
  setParameter("AbsNorm",&m_norm);

  m_maxIter = 1;
  setParameter("MaxIter",&m_maxIter);
}

//////////////////////////////////////////////////////////////////////////////

AbsoluteNormAndMaxIter::~AbsoluteNormAndMaxIter()
{
}

//////////////////////////////////////////////////////////////////////////////

bool AbsoluteNormAndMaxIter::IsGlobal () const
{
  return false;
}

//////////////////////////////////////////////////////////////////////////////

bool AbsoluteNormAndMaxIter::isAchieved(const ConvergenceStatus& status)
{
  CFLog(VERBOSE, "AbsoluteNormAndMaxIter::isAchieved() => status.res = "
	<< status.res << ", status.iter = " << m_maxIter << "\n");
  if (SubSystemStatusStack::getActive()->getStopSimulation()) return true;
  // if no iteration has been done, dont evalute the residual
  if (status.iter == 0) { return (status.iter >= m_maxIter); }
  return ( ( status.res < m_norm ) || (status.iter >= m_maxIter) );
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
