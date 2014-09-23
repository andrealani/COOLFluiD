// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Environment/CFEnv.hh"
#include "Common/CFLog.hh"
#include "Environment/Factory.hh"
#include "Config/ConfigObject.hh"
#include "Environment/DirPaths.hh"
#include "ConvertStructMesh.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Config;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::CFmeshTools;

//////////////////////////////////////////////////////////////////////////////

class AppOptions : public ConfigObject {
public:

    static void defineConfigOptions(Config::OptionList& options)
    {
      options.addConfigOption< CFuint >("log","Logging level for output");
      options.addConfigOption< std::string >("bdir", "The base dir to preced all file and dir paths.");
      options.addConfigOption< std::vector<std::string> > ("ldir", "The list of dir paths where to search for module plug-in libraries. Paths are fully qualified and not preceeded with base dir.");
      options.addConfigOption< bool >("help", "Show this help");
      options.addConfigOption< std::string >("converter", "Converter name.");
      options.addConfigOption< std::string >("meshFile", "Mesh file");
      options.addConfigOption< std::string >("meshType", "Mesh type (FEM or FVM");
    }

  AppOptions() : ConfigObject("COOLFluiD converter from block structured mesh to unstructured")
  {
    addConfigOptionsTo(this);

    logLevel = INFO;
    setParameter("log", &logLevel);

    showUsage = false;
    setParameter("help", &showUsage);

    baseDir = "";
    setParameter("bdir", &baseDir);

    libDir = vector<std::string>();
    setParameter("ldir", &libDir);

    converterName = "";
    setParameter("converter", &converterName);

    meshFile = "";
    setParameter("meshFile", &meshFile);

    meshType = "";
    setParameter("meshType", &meshType);
  }

public:
  CFuint   logLevel;
  bool     showUsage;
  std::string baseDir;
  std::vector<std::string> libDir;
  std::string converterName;
  std::string meshFile;
  std::string meshType;
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
// export LD_LIBRARY_PATH="/vki_phd/lani/COOLFLUID/debug/dso/:$LD_LIBRARY_PATH"
// ./convertStructMesh --bdir /vki_phd/lani/COOLFLUID/ --meshFile cone256 --converter Quad --meshType FVM

int main(int argc, char** argv)
{
  // Initialize Argument parser with strictArgs and helpOnError.
  AppOptions options;
  parseOptions(options, argc, argv);

  CFLogger::getInstance().setMainLoggerLevel(options.logLevel);
  Environment::CFEnv::getInstance().initiate(argc, argv);
  Environment::DirPaths::getInstance().setBaseDir(options.baseDir);
  Environment::DirPaths::getInstance().addModuleDirs(options.libDir);

  SelfRegistPtr<ConvertStructMesh> obj =
    Environment::Factory<ConvertStructMesh>::getInstance().getProvider
    (options.converterName)->create(options.meshFile);
  obj->convert(options.meshType);

  Environment::CFEnv::getInstance().terminate();

  return 0;
}


