// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include <boost/algorithm/string.hpp>

#include "Common/EventHandler.hh"

#include "Environment/ObjectProvider.hh"

#include "Framework/SMaestro.hh"
#include "Framework/SimulationStatus.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/GlobalStopCriteria.hh"
#include "Framework/Framework.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<SMaestro, Maestro, FrameworkLib, 1>
sMaestroProvider("SimpleMaestro");

//////////////////////////////////////////////////////////////////////////////

SMaestro::SMaestro(const std::string& name) : Maestro(name)
{
#if defined CF_HAVE_BOOST_1_76 || defined CF_HAVE_BOOST_1_79
   create_signal ( "control" , "Take full control of the simulation" )->connect( boost::bind ( &SMaestro::control, this, std::placeholders::_1 ) );
#else
   create_signal ( "control" , "Take full control of the simulation" )->connect( boost::bind ( &SMaestro::control, this, _1 ) );
#endif
}

//////////////////////////////////////////////////////////////////////////////

SMaestro::~SMaestro()
{
}

//////////////////////////////////////////////////////////////////////////////

Common::Signal::return_t SMaestro::control ( Common::Signal::arg_t input )
{
  CFLog(VERBOSE, "SMaestro::control() start\n");
  
  /// @todo change in input to xml
  ///       now we assume name type format
  boost::trim(input);
  std::vector< string > split_strs;
  boost::split( split_strs, input, boost::is_any_of(" "), boost::token_compress_on );

  cf_assert ( split_strs.size()>0 );
  cf_assert ( ! (split_strs.size()%2) );

  std::vector< string > subsysnames;
  std::vector< string > subsystypes;
  for ( CFuint i = 0; i < split_strs.size(); ++i, ++i )
  {
    subsysnames.push_back ( split_strs[i] );
    subsystypes.push_back ( split_strs[i+1] );
  }

  cf_assert(subsysnames.size() == subsystypes.size());

  Common::SafePtr<EventHandler> event_handler = Environment::CFEnv::getInstance().getEventHandler();

  std::vector<std::string>::const_iterator subSysName = subsysnames.begin();
  std::vector<std::string>::const_iterator subSysType = subsystypes.begin();
  for ( ;subSysName!=subsysnames.end(); ++subSysName, ++subSysType)
  {
    SimulationStatus& simStatus = SimulationStatus::getInstance();
    simStatus.resetAll();
    simStatus.setSubSystems(subsysnames);

    Common::Signal::arg_t msg;

    msg += *subSysName + "\n";
    msg += *subSysType + "\n";

    SubSystemStatusStack::setCurrentName(*subSysName);
    cf_assert(*subSysName == SubSystemStatusStack::getCurrentName());
    
    CFLog(INFO, "#\n###### STARTING SUBSYSTEM [" << *subSysName << "] ######\n#\n");
    event_handler->call_signal (event_handler->key("", "CF_ON_MAESTRO_BUILDSUBSYSTEM"), msg );
    
    CFLog(INFO, "#\n###### CONFIG PHASE #################\n#\n");
    event_handler->call_signal (event_handler->key("", "CF_ON_MAESTRO_CONFIGSUBSYSTEM"), msg );

    CFLog(INFO, "#\n###### SOCKETS PLUG PHASE ###########\n#\n");
    event_handler->call_signal (event_handler->key(*subSysName, "CF_ON_MAESTRO_PLUGSOCKETS"), msg );

    CFLog(INFO, "#\n###### BUILD PHASE ##################\n#\n");
    event_handler->call_signal (event_handler->key(*subSysName, "CF_ON_MAESTRO_BUILDPHYSICALMODEL"), msg );
    event_handler->call_signal (event_handler->key(*subSysName, "CF_ON_MAESTRO_BUILDMESHDATA"), msg );
    
    CFLog(INFO, "#\n###### SETUP PHASE ##################\n#\n");
    event_handler->call_signal (event_handler->key(*subSysName, "CF_ON_MAESTRO_SETUP"), msg );

    CFLog(INFO, "#\n###### RUN PHASE ####################\n#\n");
    //for ( ; !m_stopcriteria->isSatisfied(); )
   // {
      simStatus.incrementNbIter();
      event_handler->call_signal (event_handler->key(*subSysName, "CF_ON_MAESTRO_RUN"), msg );
   // }

    CFLog(INFO, "#\n###### UNSETUP PHASE ################\n#\n");
    event_handler->call_signal (event_handler->key(*subSysName, "CF_ON_MAESTRO_UNSETUP"), msg );

    CFLog(INFO, "#\n###### SOCKETS UNPLUG PHASE #########\n#\n");
    event_handler->call_signal (event_handler->key(*subSysName, "CF_ON_MAESTRO_UNPLUGSOCKETS"), msg );

    CFLog(INFO, "#\n###### DESTRUCTION SUBSYSTEM PHASE #########\n#\n");
    event_handler->call_signal (event_handler->key("", "CF_ON_MAESTRO_DESTROYSUBSYSTEM"), msg );
  }
  
  CFLog(VERBOSE, "SMaestro::control() end\n");
  
  return Common::Signal::return_t();
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

