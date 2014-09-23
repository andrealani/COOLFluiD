#include <iostream>

#include "Environment/CFEnv.hh"
#include "Environment/Factory.hh"
#include "Environment/DirPaths.hh"

#include "AutoTemplateLoader/convert_to_string.hh"
#include "AutoTemplateLoader/TemplateClassLoader.hh"

#include "AutoTemplateLoader/BaseRunLib.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Environment;
using namespace COOLFluiD::Config;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

class AppOptions : public ConfigObject {
public:

    static void defineConfigOptions(Config::OptionList& options)
    {
      options.addConfigOption< CFuint >("log","Logging level for output");
      options.addConfigOption< std::string >("bdir", "The base dir to preced all file and dir paths.");
      options.addConfigOption< bool >("help", "Show this help");
      options.addConfigOption< bool >("build", "Build the benchmark before running");
      options.addConfigOption< std::string >("bench", "The name of the benchmark");
      options.addConfigOption< unsigned int >("size", "The size of the vectors to use");
      options.addConfigOption< unsigned int >("nbelems", "The number of elements in the benchmark");
    }

  AppOptions() : ConfigObject("BenchmarkLoader")
  {
    addConfigOptionsTo(this);

    logLevel = INFO;
    setParameter("log", &logLevel);

    showUsage = false;
    setParameter("help", &showUsage);

    baseDir = "";
    setParameter("bdir", &baseDir);

    build_bench = false;
    setParameter("build", &build_bench);

    benchmark_class = "";
    setParameter("bench", &benchmark_class);

    size = 3;
    setParameter("size", &size);
    nbelems = 1;
    setParameter("nbelems", &nbelems);
  }

public:
  CFuint       logLevel;
  bool         showUsage;
  std::string     baseDir;
  bool         build_bench;
  std::string     benchmark_class;
  unsigned int size;
  unsigned int nbelems;
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
//   unusedArgs = options.getOptionList().getUnprocessedArgs();
//   if (unusedArgs.size())
//   {
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

int main(int argc, char** argv)
{
  // Initialize Argument parser with strictArgs and helpOnError.
  AppOptions options;
  parseOptions(options, argc, argv);

  CFLogger::getInstance().setMainLoggerLevel(options.logLevel);
  Environment::CFEnv::getInstance().initiate(argc, argv);
  Environment::DirPaths::getInstance().setBaseDir(options.baseDir);

  // implementatios goes here

  if (options.benchmark_class.empty())
  {
    options.showUsage = true;
    clog << "Warning: no benchmark defined!" << endl;
  }

  // showUsage if --help
  if (options.showUsage)
  {
    cout << options.writeUsage();
    exit(0);
  }

  CFLog(VERBOSE,"Starting ClassLoader\n");

   string noheader = "";
   double time;
   int return_value = 0;
   try
   {
      COOLFluiD::Common::SelfRegistPtr<BaseRunLib> benchmark;

      if (options.build_bench)
      {
        TemplateClassLoader<BaseRunLib,LAMGCC> benchmark_loader;
        benchmark_loader.Compiler().addToIncludeFlags(" -O2 -fomit-frame-pointer -DNDEBUG -I/home/tlmq/workspace/COOLFluiD/DEV/FRAMEWORK/src "
                                                      " -I/home/tlmq/workspace/COOLFluiD/DEV/FRAMEWORK/plugins "
                                                      " -I/home/tlmq/workspace/COOLFluiD/DEV/FRAMEWORK/plugins/AutoTemplateLoader "
                                                      " -I/home/tlmq/workspace/Develop/ForeignCodes/FLENS-2007-07-13 "
                                                      " -I/home/tlmq/workspace/Develop/ForeignCodes/blitz-0.9 "
                                                      " -I/home/tlmq/local/arch/local/tvmet/include "
                                                      " -DCF_HAVE_CONFIG_H "
                                                      " -I/home/tlmq/workspace/COOLFluiD/DEV/FRAMEWORK "
                                                      " -I/home/tlmq/workspace/COOLFluiD/DEV/Plugins "
                                                      " -I/home/tlmq/local/arch/include"
                                                      " -I/home/tlmq/workspace/COOLFluiD/DEV/FRAMEWORK/debug -I/home/tlmq/workspace/COOLFluiD/DEV/FRAMEWORK/src "
                                                      " -I/home/tlmq/workspace/COOLFluiD/DEV/FRAMEWORK/plugins "
                                                      " -I/home/tlmq/local/arch/local/blitz/include "
                                                      " -I/usr/include/atlas ");
        benchmark_loader.Compiler().addToLinkerFlags( " -L/home/tlmq/workspace/Develop/ForeignCodes/FLENS-2007-07-13 "
                                                      " -lflens "
                                                      " -Wl,-rpath -Wl,/home/tlmq/workspace/Develop/ForeignCodes/FLENS-2007-07-13");

        benchmark_loader.Compiler().addToIncludeFlags(std::string("-I") + DirPaths::getInstance().getBaseDir().string());

        // test benchmark on libraries
        TemplateClassLoader<BaseRunLib>::ParamsType benchmark_parms;
        benchmark_parms.push_back(make_pair("double",noheader));
        benchmark_parms.push_back(make_pair(convert_to_string(options.size),noheader));

        // blitz benchmark
        benchmark_loader.loadClass(options.benchmark_class,benchmark_parms);
        benchmark = benchmark_loader.createObject(options.benchmark_class,benchmark_parms);
      }
      else // just allocate the benchmark form the factory
      {
        COOLFluiD::Common::SafePtr<BaseRunLib::PROVIDER> prov =
          COOLFluiD::Environment::Factory<BaseRunLib>::getInstance().getProvider(options.benchmark_class);
        benchmark = prov->create();
      }

      time = benchmark->test(options.nbelems);
      cout << "Benchmark to library '" << options.benchmark_class << "' took " << time << " seconds" << endl;

      // deallocate the banchmark
      benchmark.release();
    }
    catch (std::exception& e)
    {
      cerr << e.what() << endl;
      cerr << "Aborting ... " << endl;
      return_value = 1;
    }
    catch (...)
    {
      cerr << "Unknown exception thrown and not caught !!!" << endl;
      cerr << "Aborting ... " << endl;
      return_value = 1;
    }

  Environment::CFEnv::getInstance().terminate();

  return return_value;
}
