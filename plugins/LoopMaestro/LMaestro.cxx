#include <boost/algorithm/string.hpp>

#include "Common/OSystem.hh"
#include "Common/EventHandler.hh"

#include "Environment/ObjectProvider.hh"
#include "Environment/DirPaths.hh"

#include "Framework/SimulationStatus.hh"
#include "Framework/GlobalStopCriteria.hh"

#include "LoopMaestro/LoopMaestro.hh"
#include "LoopMaestro/LMaestro.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace boost::filesystem;

namespace COOLFluiD {

  using namespace Framework;
  using namespace Common;
  using namespace Environment;

    namespace LoopMaestro {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<LMaestro,
               Maestro,
               LoopMaestroPlugin,
               1>
loopMaestroProvider("LoopMaestro");

//////////////////////////////////////////////////////////////////////////////

void LMaestro::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< std::vector<std::string> >("InitialFiles","Files necessary for the first iteration.");
}
      
//////////////////////////////////////////////////////////////////////////////

LMaestro::LMaestro(const std::string& name) : Maestro(name)
{
  addConfigOptionsTo(this);
  
  m_init_files = std::vector<std::string>();
  setParameter("InitialFiles",&m_init_files);
}

//////////////////////////////////////////////////////////////////////////////

LMaestro::~LMaestro()
{
}

//////////////////////////////////////////////////////////////////////////////

Common::Signal::return_t LMaestro::control ( Common::Signal::arg_t input )
{
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

  SimulationStatus& simStatus = SimulationStatus::getInstance();
  simStatus.resetAll();
  simStatus.setSubSystems(subsysnames);

  // copying the initial files to the working dir
  for(CFuint iFile=0; iFile < m_init_files.size(); ++iFile)
  {
      path user_file = path(m_init_files[iFile]);

      boost::filesystem::path from_file = DirPaths::getInstance().getWorkingDir() / user_file;
      boost::filesystem::path to_file   = DirPaths::getInstance().getResultsDir() / user_file.filename();

      // delete destination file if exsits to avoid exception throw
      if ( boost::filesystem::exists (to_file) )
        boost::filesystem::remove (to_file);

      boost::filesystem::copy_file ( from_file, to_file );
  }

  // configuring and setting Up the subSystems
  for ( ; !m_stopcriteria->isSatisfied(); )
  {

    simStatus.incrementNbIter();

    for(CFuint iter=0;iter<subsysnames.size();++iter)
    {

      std::string msg;

      msg += subsysnames[iter] + "\n";
      msg += subsystypes[iter] + "\n";

      CFout << "#\n###### STARTING SUBSYSTEM : " << subsysnames[iter] << " ######\n#\n";
      event_handler->call_signal ( "CF_ON_MAESTRO_BUILDSUBSYSTEM", msg );

      CFout << "#\n###### CONFIG PHASE #################\n#\n";
      event_handler->call_signal ( "CF_ON_MAESTRO_CONFIGSUBSYSTEM", msg );

      CFout << "#\n###### SOCKETS PLUG PHASE ###########\n#\n";
      event_handler->call_signal ( "CF_ON_MAESTRO_PLUGSOCKETS", msg );

      // allow to restart from the previous saved iteration
      if ((simStatus.isRestart()) && (simStatus.getNbIter() > 1))
      {
        CFout << "#\n### MODIFY RESTART \n#\n";
        event_handler->call_signal ( "CF_ON_MAESTRO_MODIFYRESTART", msg );
      }

      CFout << "#\n###### BUILD PHASE ##################\n#\n";
      event_handler->call_signal ( "CF_ON_MAESTRO_BUILDMESHDATA", msg );

      CFout << "#\n###### SETUP PHASE ##################\n#\n";
      event_handler->call_signal ( "CF_ON_MAESTRO_SETUP", msg );

      CFout << "#\n###### RUN PHASE ####################\n#\n";
      event_handler->call_signal ( "CF_ON_MAESTRO_RUN", msg );

      CFout << "#\n###### UNSETUP PHASE ################\n#\n";
      event_handler->call_signal ( "CF_ON_MAESTRO_UNSETUP", msg );

      CFout << "#\n###### SOCKETS UNPLUG PHASE #########\n#\n";
      event_handler->call_signal ( "CF_ON_MAESTRO_UNPLUGSOCKETS", msg );

      CFout << "#\n###### DESTRUCTION SUBSYSTEM PHASE #########\n#\n";
      event_handler->call_signal ( "CF_ON_MAESTRO_DESTROYSUBSYSTEM", msg );
    }
  }

  return Common::Signal::return_t();
}

//////////////////////////////////////////////////////////////////////////////

} // namespace LoopMaestro
} // namespace COOLFluiD

