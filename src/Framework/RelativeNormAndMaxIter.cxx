// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "MathTools/MathConsts.hh"
#include "Environment/ObjectProvider.hh"

#include "Framework/RelativeNormAndMaxIter.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/Framework.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::MathTools;

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<RelativeNormAndMaxIter, StopCondition, FrameworkLib, 1>
aRelativeNormAndMaxIterProvider("RelativeNormAndMaxIter");

//////////////////////////////////////////////////////////////////////////////

void RelativeNormAndMaxIter::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< CFreal, Config::DynamicOption<> >("RelativeNorm","Target value of the relative residual norm.");
   options.addConfigOption< CFuint, Config::DynamicOption<> >("MaxIter","Maximum number of iterations to perform.");
   options.addConfigOption< bool   >("Warn","Warn when max iterations is reached.");
}

//////////////////////////////////////////////////////////////////////////////

RelativeNormAndMaxIter::RelativeNormAndMaxIter(const std::string& name) :
  StopCondition(name)
{
  addConfigOptionsTo(this);

  m_relNorm = -8.0;
  setParameter("RelativeNorm",&m_relNorm);

  m_firstNorm = MathTools::MathConsts::CFrealMax();

  m_maxIter = std::numeric_limits<CFuint>::max();
  setParameter("MaxIter",&m_maxIter);

  m_warn = false;
  setParameter("Warn",&m_warn);
}

//////////////////////////////////////////////////////////////////////////////

RelativeNormAndMaxIter::~RelativeNormAndMaxIter()
{
}

//////////////////////////////////////////////////////////////////////////////

bool RelativeNormAndMaxIter::IsGlobal () const
{
  return false;
}

//////////////////////////////////////////////////////////////////////////////

bool RelativeNormAndMaxIter::isAchieved(const ConvergenceStatus& status)
{
  if (SubSystemStatusStack::getActive()->getStopSimulation()) return true;
  const CFuint start = 0; // first residual comes after first iteration
  const CFuint first = 1; // first residual comes after first iteration


  if (status.iter == start) { return (status.iter >= m_maxIter); }

  if (status.iter == first)
  {
    m_firstNorm = status.res;
  }

  if ( status.iter >= m_maxIter && m_warn )
    CFLog(INFO, "!!! Max number iterations [" << m_maxIter << "] reached !!!\n");

  const CFreal diff = status.res - m_firstNorm;
  bool return_value = ( (diff < m_relNorm) || (status.iter >= m_maxIter) );
  if ( return_value && m_warn )
    CFLog(INFO, "Stop condition reached at iteration [" << status.iter << "]\n");

  return return_value;
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
