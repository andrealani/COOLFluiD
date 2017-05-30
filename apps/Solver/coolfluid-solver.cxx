// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include <iostream>
#include <fstream>

#include "boost/filesystem/operations.hpp" // includes boost/filesystem/path.hpp
#include "boost/filesystem/fstream.hpp"    
#include "boost/regex.hpp"                 

#include "Common/Exception.hh"
#include "Common/OSystem.hh"
#include "Common/ProcessInfo.hh"
#include "Common/PE.hh"
#include "Common/StringOps.hh"
#include "Common/BadValueException.hh"

#include "Config/ConfigObject.hh"
#include "Config/ConfigFileReader.hh"

#include "MathTools/MathConsts.hh"

#include "Environment/CFEnv.hh"
#include "Environment/CFEnvVars.hh"
#include "Environment/DirPaths.hh"
#include "Environment/Factory.hh"

#include "Framework/SimulationStatus.hh"
#include "Framework/Simulator.hh"
#include "Framework/Maestro.hh"
#include "Framework/SubSystemStatus.hh"

#ifdef CF_HAVE_CUDA
#include "Framework/CudaDeviceManager.hh"
#endif

#ifdef CF_HAVE_SINGLE_EXEC
#include "PluginsRegister.hh"
#endif

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Environment;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Config;

//////////////////////////////////////////////////////////////////////////////

class AppOptions : public Config::ConfigObject {
  public: // functions

    /// Defines the command line options
    static void defineConfigOptions(Config::OptionList& options)
    {
      options.addConfigOption< bool >     ("wait",    "Wait before continuing with execution");
      options.addConfigOption< bool >     ("waitend",    "Wait before ending execution");
      options.addConfigOption< bool >     ("help",    "Show this help");
      options.addConfigOption< bool >     ("testEnv",    "test only the CF environment");
      options.addConfigOption< CFuint >   ("log",     "Logging level for output");
      options.addConfigOption< CFreal >   ("residual",  "The simulation should finish with this residual. The executable will return an erro if not.");
      options.addConfigOption< CFreal >   ("tolerance",   "Acceptable percentage in the testcase residual.");
      options.addConfigOption< std::string > ("scase",   "Path to simulation CFcase file");
      options.addConfigOption< std::string > ("conf",    "Path to solver environment conf file");
      options.addConfigOption< std::string > ("bdir",    "Path of base dir to preced all file and dir paths.");
      options.addConfigOption< std::vector<std::string> >("ldir", "The list of dir paths where to search for module plug-in libraries. Paths are fully qualified and not preceeded with base dir.");
    }

    /// Constructor
    AppOptions() : ConfigObject("Solver"),
      wait(false),
      waitend(false),
      show_help(false),
      test_env(false),
      log_level(INFO),
      residual(MathConsts::CFrealMax()),
      tolerance(3.0),
      scase_file(""),
      conf_file("coolfluid-solver.xml"),
      base_dir(""),
      libDir(0)
    {
      addConfigOptionsTo(this);
      setParameter( "wait",   &wait        );
      setParameter( "waitend",   &waitend        );
      setParameter( "help",   &show_help    );
      setParameter( "testEnv",   &test_env    );
      setParameter( "log",    &log_level    );
      setParameter( "tolerance", &tolerance    );
      setParameter( "residual", &residual  );
      setParameter( "scase",  &scase_file );
      setParameter( "conf",   &conf_file );
      setParameter( "bdir",   &base_dir     );
      setParameter( "ldir",   &libDir)     ;
    }

    /// Checks the validity of the command line options
    void check_options ()
    {
        if ( solver_conf.size() != 0 ) // show help if any argument is unprocessed
        {
          ConfigArgs::iterator itr = solver_conf.begin();
          for(; itr != solver_conf.end(); ++itr)
          {
            clog << "Unused argument: " << itr->first << " " << itr->second << endl;
          }
          CFLog(INFO, writeUsage());
          if (waitend) cin.get();
          exit(1);
        }

        if (scase_file == "") // no case selected
        {
          show_help = true;
          clog << "Warning: Simulation CFcase case file is not defined!" << endl;
        }

        if (show_help) // show_help if --help
        {
          CFLog(INFO, writeUsage());
          if (waitend) cin.get();
          exit(0);
        }

    }

  public: // data
      std::string mCFcaseFile;
      bool   wait;
      bool   waitend;
      bool   show_help;
      bool   test_env;
      CFuint  log_level;
      CFreal  residual;
      CFreal  tolerance;
      std::string scase_file;
      std::string conf_file;
      std::string base_dir;
      std::vector<std::string> libDir;
      /// where to put the configuration options
      ConfigArgs solver_conf;
};

//////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv)
{
  using namespace boost;
  
#ifdef CF_HAVE_SINGLE_EXEC
  Common::SafePtr<FactoryRegistry> fRegistry = 
    Environment::CFEnv::getInstance().getFactoryRegistry();
  cf_assert(fRegistry.isNotNull());
  
  PluginsRegister pr;
  pr.registerAll(fRegistry);
#else
  Common::SafePtr<FactoryRegistry> fRegistry(CFNULL);
#endif
  
  // process the command line
  AppOptions options;
  options.getOptionList().setStrictArgs(true);
  try
  {
    options.getOptionList().processCommandLine(argc, argv);
    options.check_options ();

    if ( filesystem::exists( options.conf_file ) )
      options.setConfigFile( options.conf_file );

    options.configure ( options.solver_conf );
  }
  // Something bad happened. Dump the usage statement:
  catch (const Exception& e)
  {
    cerr << e.what() << "\n";
    cerr << options.writeUsage();
    exit(1);
  }

  int return_value = 0;
  try
  {
    // set the command line log level
    // takes precedence over the config
    CFLogger::getInstance().setMainLoggerLevel(options.log_level);
    
    // build the environment
    CFEnv& cf_env = CFEnv::getInstance();
    
    // initiate the environemnt
    cf_env.initiate(argc, argv);

    CFLog(VERBOSE, "--- coolfluid-solver ----------------------------------------\n\n");
 
    // print out starting directory
    CFLog(VERBOSE, "starting in directory [" << filesystem::current_path().string() << "]\n\n");
    
    // print out input parameters 
    CFLog(VERBOSE, "called with arguments:\n");
    for ( int iarg = 0; iarg < argc; ++iarg ) {
      CFLog(VERBOSE, "arg [" << iarg << "] : [" << argv[iarg] << "]\n");
     }
    CFLog(VERBOSE, "\n-------------------------------------------------------------\n");

    CFLog(VERBOSE, "places to search for libraries ...\n[");
    std::vector< string >::const_iterator itr = options.libDir.begin();
    for ( ; itr != options.libDir.end() ; ++itr )
    {
       CFLog(VERBOSE, "\n\tlibpath [" << *itr << "]");
     }   
     CFLog(VERBOSE, "\n]");
     CFLog(VERBOSE, "\n-------------------------------------------------------------\n\n");
   
    //////////////////////////////////////////////////////////////////
    // removes the config-p*.log
    // only processor 0 does it while the others wait
    PE::GetPE().setBarrier("Default");
    if ( cf_env.getCPURank() == 0 )
    {
      CFLog(VERBOSE, "Removing previous config logs\n");
      
      filesystem::directory_iterator ditr(filesystem::current_path()), dir_end;
      for( ; ditr != dir_end; ++ditr )
      {

#ifdef CF_HAVE_BOOST_1_42
	std::string filename = ditr->filename();
#else
	std::string filename = (ditr->path()).filename().native();
#endif

	bool is_config = StringOps::startsWith(filename, "config") && StringOps::endsWith(filename, ".log"); 
	bool is_output = StringOps::endsWith(filename, "-output.log"); 
	
	if ( !is_directory( *ditr ) && ( is_config || is_output ) )
	{
	  CFLog(VERBOSE, "removing file: " << filename << "\n");
	  try { filesystem::remove(*ditr); } catch (...) {};
	}
      }
    }
    PE::GetPE().setBarrier("Default");
    
    //////////////////////////////////////////////////////////////////
    
    // set essential paths
    DirPaths::getInstance().setBaseDir(options.base_dir);
    DirPaths::getInstance().addModuleDirs(options.libDir);
    
    // configure the environment
    filesystem::path casefile;
    if ( !Common::StringOps::startsWith( options.scase_file,".") && !Common::StringOps::startsWith(options.scase_file,"/") )
      casefile  = Environment::DirPaths::getInstance().getBaseDir() / filesystem::path(options.scase_file);
    else
      casefile  =  filesystem::path(options.scase_file);
    
    // parse configuration options for environment
    ConfigFileReader cfile_reader;
    ConfigArgs config_args;
    cfile_reader.parse ( casefile, config_args);
    // configure the environemt
    cf_env.configure ( config_args );
    // setup the runtime environment
    cf_env.setup();
    
    if (!options.test_env) {
      
#ifdef CF_HAVE_CUDA
    CudaEnv::CudaDeviceManager& cudaDev = CudaEnv::CudaDeviceManager::getInstance();
    cudaDev.configure(config_args);
#endif
    
    if (options.wait) { // wait to attach debugger
      CFuint pid = OSystem::getInstance().getProcessInfo()->getPID();
      CFLog(INFO, "Stopping to attach debugger ...\n");
      CFLog(INFO, "Current PID is [" << pid << "].\n");
      char ans;
      bool first_pass = true;
      do {
	if (!first_pass) { CFLog(INFO, "Please type a 'y' or an 'n'.\n"); }
	CFout << "Continue (y/n) ?\n";
	std::cin >> ans;
	first_pass = false;
      }
      while((ans !='Y')&&(ans !='N')&&(ans !='y')&&(ans !='n'));
    }
    
    CFLog(VERBOSE,"-------------------------------------------------------------\n");
    
    // create the maestro
    
    Common::SelfRegistPtr<Maestro> maestro;
    std::string maestro_str = "SimpleMaestro";
    if ( config_args.find ("Maestro") != config_args.end() )
      maestro_str = config_args["Maestro"];
    
    CFLog(VERBOSE,"Creating Simulation Maestro\n");
    Common::SafePtr<Maestro::PROVIDER> prov = 
      FACTORY_GET_PROVIDER(fRegistry, Maestro, maestro_str);
    CFLog(VERBOSE,"Creating Simulation Maestro 1\n");
    cf_assert(prov.isNotNull());
    CFLog(VERBOSE,"Creating Simulation Maestro 2\n");
    maestro.reset(prov->create(prov->getName()));
    maestro->setFactoryRegistry(fRegistry);
    CFLog(VERBOSE,"Creating Simulation Maestro 4\n");
    maestro->configure ( config_args );
    CFLog(VERBOSE,"Creating Simulation Maestro 5\n");
    
    CFLog(INFO,"-------------------------------------------------------------\n");
    
    // create the simulation
    
    CFLog(INFO,"Creating Simulation\n");
    // create the simulator
    SharedPtr < Simulator > sim ( new Simulator("Simulator") );
    sim->setFactoryRegistry(fRegistry);
    // give him the file to configure from
    sim->openCaseFile(casefile.string());
    
    CFLog(INFO,"-------------------------------------------------------------\n");
    
//////////////////////////////////////////////////////////////////
// dumps the tree.xml
   /* std::string sim_xmltree = sim->getTreeXML();
    ofstream xml;
    xml.open ("tree.xml");
    xml << sim_xmltree << std::endl;
    xml.close();*/
//////////////////////////////////////////////////////////////////

    // maestro will manage this simulation
    maestro->manage ( sim );

    // which subsystems will be controlled by meastro
    std::string msg;
    std::vector< string > subsysnames = sim->getSubSystemNames();
    std::vector< string > subsystypes = sim->getSubSystemTypes();
    cf_assert ( subsysnames.size() == subsystypes.size() );
    for ( CFuint i = 0; i < subsysnames.size(); ++i )
    {
      msg += subsysnames[i] + " ";
      msg += subsystypes[i] + " ";
    }

    // maestro takes control of simulation
    maestro->call_signal ( "control", msg );
    
    // if a target residual has been set in the CFcase file, run the following test
    if ( options.residual != MathConsts::CFrealMax()) {
      const CFreal totalResidual = SimulationStatus::getInstance().getLastResidual();
      cf_assert(std::abs(totalResidual) > 0.);
      const CFreal tolerance = std::abs(options.residual - totalResidual) * 100.0 / std::abs(options.residual);
      
      CFLog(INFO, "\n"
	    << "Target   residual [" << options.residual << "]\n"
	    << "Achieved residual [" << totalResidual << "]\n"
	    << "Target  diff      [" << options.tolerance << "%]\n"
	    << "Percent diff      [" << tolerance << "%]\n\n");
      
      // test succeed if the tolerance is smaller than the threshold one
      return_value = (tolerance > options.tolerance) ? 1 : 0;
      
      if (return_value == 1) {
	const string msg = "Test failed with a " + StringOps::to_str(tolerance) + 
	  "% error (max allowed is " + StringOps::to_str(options.tolerance) + "%)";
	throw Common::BadValueException(FromHere(), msg);
      }
    }
    
    // deallocate Simulator
    sim.release();
    // deallocate Maestro
    maestro.release();
    }
    
    // unsetup the runtime environment
    cf_env.unsetup();
    // terminate the runtime environment
    cf_env.terminate();
  }
  catch ( std::exception& e ) {
    cerr << e.what() << endl;
    cerr << "Exception thrown: Aborting ..." << endl;
    return_value = 1;
  }
  catch (...) {
    cerr << "Unknown exception thrown and not caught !!!" << endl;
    cerr << "Aborting ... " << endl;
    return_value = 1;
  }
  
  CFLog(VERBOSE, "Exit value " << return_value << "\n");
  if (options.waitend) cin.get();
  
  return return_value;
}
