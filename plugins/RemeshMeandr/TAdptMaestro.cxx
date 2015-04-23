#include <boost/algorithm/string.hpp>

#include "Common/EventHandler.hh"
#include "Common/OSystem.hh"
#include "Environment/DirPaths.hh"
#include "Environment/ObjectProvider.hh"
#include "Framework/SimulationStatus.hh"
#include "Framework/GlobalStopCriteria.hh"
#include "Framework/SubSystemStatus.hh"

#include "RemeshMeandr/TAdptMaestro.hh"
#include "RemeshMeandr/RemeshMeandr.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace boost::filesystem;

namespace COOLFluiD {

  using namespace Framework;
  using namespace Common;
  using namespace Environment;

  namespace Numerics {

    namespace RemeshMeandros {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<TAdptMaestro,Maestro, RemeshMeandrModule, 1> tadaptMaestroProvider("TAdptMaestro");

//////////////////////////////////////////////////////////////////////////////

void TAdptMaestro::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< std::vector<std::string> >("InitialFiles","Files necessary for the first iteration.");
}

//////////////////////////////////////////////////////////////////////////////

TAdptMaestro::TAdptMaestro(const std::string& name) : Maestro(name)
{
   addConfigOptionsTo(this);

   _initialFilesStr = std::vector<std::string>();
   setParameter("InitialFiles",&_initialFilesStr);
}

//////////////////////////////////////////////////////////////////////////////

TAdptMaestro::~TAdptMaestro()
{
}

//////////////////////////////////////////////////////////////////////////////

Common::Signal::return_t TAdptMaestro::control ( Common::Signal::arg_t input )
{
  /// @todo change in input to xml
  ///       now we assume name type format

  using namespace std;

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

  SimulationStatus& simStatus = SimulationStatus::getInstance();
  simStatus.resetAll();
  simStatus.setSubSystems(subsysnames);

  // copying the initial files to the working dir
  for ( CFuint iFile=0;iFile < _initialFilesStr.size(); ++iFile )
  {
      path user_file = path(_initialFilesStr[iFile]);

      boost::filesystem::path from_file = DirPaths::getInstance().getWorkingDir() / user_file;
      boost::filesystem::path to_file   = DirPaths::getInstance().getResultsDir() / user_file.filename();

      // delete destination file if exsits to avoid exception throw
      if ( boost::filesystem::exists (to_file) )
        boost::filesystem::remove (to_file);

      boost::filesystem::copy_file ( from_file, to_file );
  }

  // configuring and setting Up the subSystems
  for ( CFint it = 0; !m_stopcriteria->isSatisfied(); ++it )
  {

  CFout << " ### --- Iteration: " << SimulationStatus::getInstance().getNbIter() << "\n";

    simStatus.incrementNbIter();

    CFuint iter = 0;

    std::string msg;

    msg += subsysnames[iter] + "\n";
    msg += subsystypes[iter] + "\n";
    
    const string currSubsysName = subsysnames[iter];
    SubSystemStatusStack::setCurrentName(currSubsysName);
    cf_assert(currSubsysName == SubSystemStatusStack::getCurrentName());
    
    CFout << "#\n###### STARTING SUBSYSTEM [" << subsysnames[iter] << "] ######\n#\n";
    event_handler->call_signal (event_handler->key("", "CF_ON_MAESTRO_BUILDSUBSYSTEM"), msg );
    
    CFout << "#\n###### CONFIG PHASE #################\n#\n";
    event_handler->call_signal (event_handler->key("","CF_ON_MAESTRO_CONFIGSUBSYSTEM"), msg );

    CFout << "#\n###### SOCKETS PLUG PHASE ###########\n#\n";
    event_handler->call_signal (event_handler->key(currSubsysName, "CF_ON_MAESTRO_PLUGSOCKETS"), msg );

    // allow to restart from the previous saved iteration
    if ((simStatus.isRestart()) && (simStatus.getNbIter() > 1))
    {
        CFout << "#\n### MODIFY RESTART \n#\n";
        event_handler->call_signal (event_handler->key(currSubsysName, "CF_ON_MAESTRO_MODIFYRESTART"), msg );
    }

    CFout << "#\n###### BUILD PHASE ##################\n#\n";
    event_handler->call_signal (event_handler->key(currSubsysName, "CF_ON_MAESTRO_BUILDMESHDATA"), msg );

    CFout << "#\n###### SETUP PHASE ##################\n#\n";
    event_handler->call_signal (event_handler->key(currSubsysName, "CF_ON_MAESTRO_SETUP"), msg );

    CFout << "#\n###### RUN PHASE ####################\n#\n";
    event_handler->call_signal (event_handler->key(currSubsysName, "CF_ON_MAESTRO_RUN"), msg );

    CFout << "#\n###### UNSETUP PHASE ################\n#\n";
    event_handler->call_signal (event_handler->key(currSubsysName, "CF_ON_MAESTRO_UNSETUP"), msg );

    CFout << "#\n###### SOCKETS UNPLUG PHASE #########\n#\n";
    event_handler->call_signal (event_handler->key(currSubsysName, "CF_ON_MAESTRO_UNPLUGSOCKETS"), msg );

    CFout << "#\n###### DESTRUCTION SUBSYSTEM PHASE #########\n#\n";
    event_handler->call_signal (event_handler->key("", "CF_ON_MAESTRO_DESTROYSUBSYSTEM"), msg );
    
    CFout << " ------------ " << it << " END \n";
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace RemeshMeandros

  } // namespace Numerics

} // namespace COOLFluiD

