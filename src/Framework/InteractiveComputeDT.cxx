// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Environment/ObjectProvider.hh"

#include "Framework/SubSystemStatus.hh"
#include "Framework/ComputeDT.hh"
#include "Framework/InteractiveComputeDT.hh"
#include "Framework/Framework.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<InteractiveComputeDT,
               ComputeDT,
               FrameworkLib,
               1>
interactiveComputeDTProvider("Interactive");

//////////////////////////////////////////////////////////////////////////////

void InteractiveComputeDT::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal, Config::DynamicOption<> >("DT","Value for the physical time step");
}

//////////////////////////////////////////////////////////////////////////////

InteractiveComputeDT::InteractiveComputeDT(const std::string& name) :
  ComputeDT(name),
  m_dt_value(0.0)
{
  addConfigOptionsTo(this);
  setParameter("DT",&m_dt_value);
}

//////////////////////////////////////////////////////////////////////////////

InteractiveComputeDT::~InteractiveComputeDT()
{
}

//////////////////////////////////////////////////////////////////////////////

void InteractiveComputeDT::configure ( Config::ConfigArgs& args )
{
  ComputeDT::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void InteractiveComputeDT::operator() ()
{
  if (m_dt_value != SubSystemStatusStack::getActive()->getDTDim())
  {
    SubSystemStatusStack::getActive()->setDTDim(m_dt_value);
  }
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
