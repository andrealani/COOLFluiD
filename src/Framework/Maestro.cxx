// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Environment/Factory.hh"

#include "Framework/Maestro.hh"
#include "Framework/GlobalStopCriteria.hh"
#include "Framework/SimulationStatus.hh"
#include "Framework/Simulator.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

void Maestro::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< std::string >("GlobalStopCriteria","The stop Criteria to control the iteration procedure.");
   options.addConfigOption< bool >("RestartFromPreviousSolution","Restart each subsystem from the previous iteration.");
   options.addConfigOption< bool >("AppendIter","Save each iteration to different Tecplot file with suffix _GlobalIter#.");
}

//////////////////////////////////////////////////////////////////////////////

Maestro::Maestro(const std::string& name) :
  OwnedObject(),
  ConfigObject(name)
{
  addConfigOptionsTo(this);

  m_stopcriteria_str = "GlobalMaxNumberSteps";
  setParameter("GlobalStopCriteria",&m_stopcriteria_str);

  m_append_iter = false;
  setParameter("AppendIter",&m_append_iter);

  m_restart_from_previous = false;
  setParameter("RestartFromPreviousSolution",&m_restart_from_previous);
}

//////////////////////////////////////////////////////////////////////////////

Maestro::~Maestro() {}

//////////////////////////////////////////////////////////////////////////////

void Maestro::configure ( Config::ConfigArgs& args )
{
  ConfigObject::configure(args);

  // builds the global stop Criteria
  Common::SelfRegistPtr<GlobalStopCriteria>* stopcriteria = FACTORY_GET_PROVIDER
    (getFactoryRegistry(), GlobalStopCriteria, m_stopcriteria_str)->
    createPtr(m_stopcriteria_str);
  
  m_stopcriteria = *stopcriteria;
  m_stopcriteria->setFactoryRegistry(getFactoryRegistry());
  configureNested ( m_stopcriteria.getPtr(), args );

  SimulationStatus::getInstance().setAppendIter(m_append_iter);
  
  const bool restart = SimulationStatus::getInstance().isRestart();
  if(!restart) SimulationStatus::getInstance().setRestart(m_restart_from_previous);
 
  delete stopcriteria;
}

//////////////////////////////////////////////////////////////////////////////

void Maestro::manage ( Common::SharedPtr<Simulator> sim )
{
  m_sim.reset ( sim );
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

