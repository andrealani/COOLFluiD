// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "MathTools/MathChecks.hh"
#include "Environment/ObjectProvider.hh"

#include "Framework/InteractiveComputeCFL.hh"
#include "Framework/CFL.hh"
#include "Framework/Framework.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<InteractiveComputeCFL,
	                          ComputeCFL,
                            FrameworkLib,
	                          1>
aInteractiveComputeCFLProvider("Interactive");

//////////////////////////////////////////////////////////////////////////////

void InteractiveComputeCFL::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal, Config::DynamicOption<> >("CFL","Value for the CFL");
}

//////////////////////////////////////////////////////////////////////////////

InteractiveComputeCFL::InteractiveComputeCFL(const std::string& name) :
  ComputeCFL(name),
  m_curr_cfl(0.0),
  m_prev_cfl(0.0)
{
  addConfigOptionsTo(this);
  setParameter("CFL",&m_curr_cfl);
}

//////////////////////////////////////////////////////////////////////////////

InteractiveComputeCFL::~InteractiveComputeCFL()
{
}

//////////////////////////////////////////////////////////////////////////////

void InteractiveComputeCFL::configure ( Config::ConfigArgs& args )
{
  ComputeCFL::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void InteractiveComputeCFL::operator() (const ConvergenceStatus& m_cstatus)
{
  if( MathTools::MathChecks::isNotEqual(m_curr_cfl, m_prev_cfl) )
  {
    m_prev_cfl = m_curr_cfl;
    _cfl->setCFLValue(m_curr_cfl);
  }
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
