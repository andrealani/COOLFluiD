#include <iostream>

#include "Environment/CFEnv.hh"
#include "MutationI/MutationLibrary.hh"
#include "Environment/Factory.hh"
#include "Environment/DirPaths.hh"
#include "ComputeWithMutation.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Config;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::Mutation;
using namespace COOLFluiD::MutationUsage;

//////////////////////////////////////////////////////////////////////////////

class AppOptions : public ConfigObject {
public:

    static void defineConfigOptions(Config::OptionList& options)
    {
      options.addConfigOption< CFuint >("log","Logging level for output");
      options.addConfigOption< std::string >("bdir", "The base dir to preced all file and dir paths.");
      options.addConfigOption< bool >("help", "Show this help");
      options.addConfigOption< std::string >("action", "Action to perform.");
      options.addConfigOption< std::string >("mix", "Gas mixture to use");
    }

  AppOptions() : ConfigObject("MutationUsage")
  {
    addConfigOptionsTo(this);

    logLevel = INFO;
    setParameter("log", &logLevel);

    showUsage = false;
    setParameter("help", &showUsage);

    baseDir = "";
    setParameter("bdir", &baseDir);

    action = "";
    setParameter("action", &action);

    mix = "air5";
    setParameter("mix", &mix);
  }

public:
  CFuint   logLevel;
  bool   showUsage;
  std::string baseDir;
  std::string action;
  std::string mix;
};

void parseOptions(AppOptions& options,
		  int argc, char** argv)
{
  options.getOptionList().setStrictArgs(true);

  // unused arguments
  ConfigArgs unusedArgs;

  // process the command line
  try {
    options.getOptionList().processCommandLine(argc, argv);
  }
  // Something bad happened. Dump the usage statement:
  catch (const Common::Exception& e) {
    cerr << options.writeUsage();
    exit(1);
  }

  // showUsage and exit if any argument is unprocessed
  /// @todo when processing command line detect which options where not used
//   unusedArgs = options.getOptionList().getUnprocessedArgs();
//   if (unusedArgs.size()!=0) {
//     ConfigArgs::iterator itr = unusedArgs.begin();
//     for(; itr != unusedArgs.end(); ++itr) {
//       clog << "Unused argument: " << itr->first << " " << itr->second << endl;
//     }
//     cout << options.writeUsage();
//     exit(1);
//   }

  // showUsage if --help
  if (options.showUsage) {
    cout << options.writeUsage();
    exit(0);
  }
}

// to run the executable do:
// ./MutationUsage --bdir /data/andrea/COOLFLUID/ --action xLTE --mix air5

int main(int argc, char** argv)
{
  // Initialize Argument parser with strictArgs and helpOnError.
  AppOptions options;
  parseOptions(options, argc, argv);

  CFLogger::getInstance().setMainLoggerLevel(options.logLevel);
  Environment::CFEnv::getInstance().initiate(argc, argv);
  Environment::DirPaths::getInstance().setBaseDir(options.baseDir);

  MutationLibrary mutLib("Mutation");
  
  std::string baseDir = options.baseDir;
  const std::string path = baseDir + "/plugins/MutationI/data/mutation/";
  mutLib.setLibPathName(path);
  mutLib.setMixtureName(options.mix);
  mutLib.setup();

  SelfRegistPtr<ComputeWithMutation> action =
    Environment::Factory<ComputeWithMutation>::getInstance().getProvider(options.action)->create(&mutLib);

  action->compute();

  Environment::CFEnv::getInstance().terminate();

  return 0;
}
