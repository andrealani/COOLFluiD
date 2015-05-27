// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/StateInterpolator.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

void StateInterpolator::defineConfigOptions(Config::OptionList& options)
{
  //options.addConfigOption< CFuint >("MonitoredVarID","ID of the variable whose residual will be monitored.");
}

//////////////////////////////////////////////////////////////////////////////

StateInterpolator::StateInterpolator(const std::string& name) :
  NullableObject(),
  NumericalStrategy(name)
{
  addConfigOptionsTo(this);
  
  //  m_compute_var_id = std::vector<CFuint>();
  // setParameter("ComputedVarID",&m_compute_var_id);
}


//////////////////////////////////////////////////////////////////////////////

StateInterpolator::~StateInterpolator()
{
}

//////////////////////////////////////////////////////////////////////////////

void StateInterpolator::setup()
{
  CFAUTOTRACE;
  
  NumericalStrategy::setup();
}

//////////////////////////////////////////////////////////////////////////////

void StateInterpolator::unsetup()
{
  CFAUTOTRACE;

  NumericalStrategy::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<Framework::BaseDataSocketSink> >
StateInterpolator::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;
  return result;
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
