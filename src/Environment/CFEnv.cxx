// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include <fstream>

#ifdef CF_HAVE_CONFIG_H
  #include "coolfluid_config.h"
#endif

#include "logcpp/FileAppender.hh"
#include "logcpp/OstreamAppender.hh"
#include "logcpp/PatternLayout.hh"

#include "coolfluid_svnversion.hh"

#include "Common/PE.hh"
#include "Common/EventHandler.hh"
#include "Common/CFLog.hh"
#include "Common/SignalHandler.hh"
#include "Common/OSystem.hh"
#include "Common/FactoryRegistry.hh"

#include "Environment/SingleBehaviorFactory.hh"
#include "Environment/DirPaths.hh"
#include "Environment/FileHandlerInput.hh"
#include "Environment/FileHandlerOutput.hh"
#include "Environment/ModuleRegistry.hh"
#include "Environment/CFEnv.hh"
#include "Environment/ModuleRegisterBase.hh"
#include "Environment/CFEnvVars.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Environment {

//////////////////////////////////////////////////////////////////////////////

CFEnv& CFEnv::getInstance()
{
  static CFEnv env;
  return env;
}

//////////////////////////////////////////////////////////////////////////////

void CFEnv::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< bool >    ("OnlyCPU0Writes",    "Only CPU0 writes to stdout");
  options.addConfigOption< bool >    ("DoAssertions",      "Turn off assertions dynamically");
  options.addConfigOption< bool >    ("AssertionDumps",    "If assertions should dump backtraces");
  options.addConfigOption< bool >    ("AssertionThrows",   "If assertions should throw exceptions instead of aborting code");
  options.addConfigOption< bool >    ("ExceptionOutputs",  "If exception contructor should output");
  options.addConfigOption< bool >    ("ExceptionDumps",    "If exception contructor should dump backtrace");
  options.addConfigOption< bool >    ("ExceptionAborts",   "If exception contructor should abort execution immedietly");
  options.addConfigOption< CFuint >  ("ExceptionLogLevel", "Loglevel for exceptions");
  options.addConfigOption< bool >    ("RegistSignalHandlers", "If CPU signal handlers should be registered");
  options.addConfigOption< bool >    ("TraceToStdOut",     "If Tracing should be sent to stdout also");
  options.addConfigOption< bool >    ("TraceActive",       "If Tracing should be active");
  options.addConfigOption< bool >    ("VerboseEvents",     "If Events have verbose output");
  options.addConfigOption< bool >    ("ErrorOnUnusedConfig","Signal error when some user provided config parameters are not used");
  options.addConfigOption< std::string >("MainLoggerFileName", "Name of main log file");
  options.addConfigOption< CFuint >("NbWriters", "Number of writing processes in parallel I/O");
}
    
//////////////////////////////////////////////////////////////////////////////
    
CFEnv::CFEnv() : 
  Config::ConfigObject("CFEnv"),
  m_eventHandler(new Common::EventHandler()),
  m_moduleRegistry(new Environment::ModuleRegistry()),
  m_factoryRegistry(new Common::FactoryRegistry()),
  m_env_vars (new CFEnvVars())
{
  addConfigOptionsTo(this);

  setParameter("DoAssertions",          &(AssertionManager::getInstance().DoAssertions));
  setParameter("AssertionDumps",        &(AssertionManager::getInstance().AssertionDumps));
  setParameter("AssertionThrows",       &(AssertionManager::getInstance().AssertionThrows));

  setParameter("ExceptionOutputs",      &(ExceptionManager::getInstance().ExceptionOutputs));
  setParameter("ExceptionDumps",        &(ExceptionManager::getInstance().ExceptionDumps));
  setParameter("ExceptionAborts",       &(ExceptionManager::getInstance().ExceptionAborts));

  setParameter("OnlyCPU0Writes",        &(m_env_vars->OnlyCPU0Writes));
  setParameter("RegistSignalHandlers",  &(m_env_vars->RegistSignalHandlers));
  setParameter("VerboseEvents",         &(m_env_vars->VerboseEvents));
  setParameter("ErrorOnUnusedConfig",   &(m_env_vars->ErrorOnUnusedConfig));
  setParameter("TraceToStdOut",         &(m_env_vars->TraceToStdOut));
  setParameter("TraceActive",           &(m_env_vars->TraceActive));
  setParameter("MainLoggerFileName",    &(m_env_vars->MainLoggerFileName));
  setParameter("ExceptionLogLevel",     &(m_env_vars->ExceptionLogLevel));
  setParameter("NbWriters",     &(m_env_vars->NbWriters));
}

//////////////////////////////////////////////////////////////////////////////

void CFEnv::configure ( Config::ConfigArgs& args )
{
  // set the number of writers for parallel I/O
  m_env_vars->NbWriters = PE::GetPE().GetProcessorCount("Default");
  
  ConfigObject::configure(args);
  
  CFLog(VERBOSE, "Configuring OSystem signal handlers ... \n");
  if ( m_env_vars->RegistSignalHandlers )
  {
    OSystem::getInstance().getSignalHandler()->registSignalHandlers();
    CFLog(VERBOSE, "OK\n");
  }
  else {
   CFLog(VERBOSE, "skipping\n");
  }

  CFLog(VERBOSE, "Configuring Logging ... \n");
  initLoggers();
  CFLog(VERBOSE, "OK\n");

  // clean the config.log file
 /* boost::filesystem::path fileconfig =
    Environment::DirPaths::getInstance().getResultsDir() / boost::filesystem::path("config.log");

  SelfRegistPtr<Environment::FileHandlerOutput> fhandle = Environment::SingleBehaviorFactory<Environment::FileHandlerOutput>::getInstance().create();
  std::ofstream& fout = fhandle->open(fileconfig,std::ios_base::trunc);

  fout << "CONFIG LOG:" << "\n";
  fhandle->close();*/
}

//////////////////////////////////////////////////////////////////////////////

void CFEnv::setup()
{
  CFLog(VERBOSE, "CFEnv::setup() => start\n");
  
  SetupObject::setup();
  
  // these are the default values
#ifdef CF_HAVE_CURL  
  SingleBehaviorFactory<Environment::FileHandlerInput>::getInstance().setDefaultBehavior("CurlAccessRepository");
#else
  SingleBehaviorFactory<Environment::FileHandlerInput>::getInstance().setDefaultBehavior("DirectFileAccess");
#endif
  
  SingleBehaviorFactory<Environment::FileHandlerOutput>::getInstance().setDefaultBehavior("DirectFileWrite");
  
  CFLog(VERBOSE, "CFEnv::setup() => end\n");
}

//////////////////////////////////////////////////////////////////////////////

void CFEnv::unsetup()
{
  SetupObject::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

std::string CFEnv::getVersionString () const
{
  std::string ret;
  ret += getCFVersion() + " ";
  ret += "Kernel " + getKernelVersion() + " ";
  ret += "( r" + getSvnVersion() + ", " + Common::PE::GetPE().GetName() + ", " + getBuildType() + " )";
  return ret;
}

//////////////////////////////////////////////////////////////////////////////

std::string CFEnv::getBuildType () const
{
  return CF_BUILD_TYPE;
}

//////////////////////////////////////////////////////////////////////////////

std::string CFEnv::getSvnVersion () const
{
  return COOLFLUID_SVNVERSION;
}

//////////////////////////////////////////////////////////////////////////////

std::string CFEnv::getCFVersion () const
{
  return COOLFLUID_VERSION_STR;
}

//////////////////////////////////////////////////////////////////////////////

std::string CFEnv::getKernelVersion () const
{
  return CF_KERNEL_VERSION_STR;
}

//////////////////////////////////////////////////////////////////////////////

std::string CFEnv::getBuildProcessor () const
{
  return CF_BUILD_PROCESSOR;
}

//////////////////////////////////////////////////////////////////////////////

std::string CFEnv::getBuildSystem () const
{
  std::string ret;
#ifdef CF_CMAKE_VERSION
  ret += "CMake ";
  ret += CF_CMAKE_VERSION;
#else
  ret += "Unknown";
#endif
  return ret;
}

//////////////////////////////////////////////////////////////////////////////

std::string CFEnv::getSystemName() const
{
  return CF_OS_NAME;
}

//////////////////////////////////////////////////////////////////////////////

std::string CFEnv::getLongSystemName() const
{
  return CF_OS_LONGNAME;
}

//////////////////////////////////////////////////////////////////////////////

std::string CFEnv::getSystemVersion() const
{
  return CF_OS_VERSION;
}

//////////////////////////////////////////////////////////////////////////////

std::string CFEnv::getSystemBits() const
{
  return CF_OS_BITS;
}
//////////////////////////////////////////////////////////////////////////////

CFEnv::~CFEnv()
{
  deletePtr ( m_env_vars );

  delete m_moduleRegistry;     m_moduleRegistry = CFNULL;
  delete m_factoryRegistry;    m_factoryRegistry = CFNULL;

/// @todo should be done like this and these classes probably
///       should be nested and private inside this class
//   deletePtr(m_moduleRegistry);
//   deletePtr(m_factoryRegistry);
}

//////////////////////////////////////////////////////////////////////////////

void CFEnv::initLoggers()
{
  CFLogger::getInstance().getTraceLogger().removeAllAppenders();
  CFLogger::getInstance().getMainLogger().removeAllAppenders();

  const bool append = false;
  std::ostringstream stout_filename;
  std::ostringstream trace_filename;
  stout_filename << "P" << getCPURank() << "-" << m_env_vars->MainLoggerFileName;
  trace_filename << "P" << getCPURank() << "-TRACE-" << m_env_vars->MainLoggerFileName;


  // OUTPUT TO SCREEN --------------------------------------------------------

  // create and set the standard output appender
  // only rank 0 CPU will write to the standard output
#ifndef CF_ENABLE_PARALLEL_VERBOSE
  if (getCPURank() == 0)
#endif
  {
    logcpp::PatternLayout* stdout_layout = new logcpp::PatternLayout();
    std::string stdout_format("%m");
    stdout_layout->setConversionPattern(stdout_format);

    logcpp::Appender* stdout_appender = new logcpp::OstreamAppender("StdoutAppender", &std::cout);
    stdout_appender->setLayout(stdout_layout);

    CFLogger::getInstance().getMainLogger().addAppender(stdout_appender);
  }

  // OUTPUT TO FILES --------------------------------------------------------

  // create and set the file output appender
  // always to cpu P0 and the others optionally
  if ( !m_env_vars->OnlyCPU0Writes || getCPURank() == 0 )
  {
    logcpp::PatternLayout* f_layout = new logcpp::PatternLayout();
//     std::string f_format("[%p] (%d{%H:%M:%S}) %m");
    std::string f_format("%m");
    f_layout->setConversionPattern(f_format);

    logcpp::Appender* f_appender = new logcpp::FileAppender("FileAppender",stout_filename.str(),append);
    f_appender->setLayout(f_layout);

    CFLogger::getInstance().getMainLogger().addAppender(f_appender);

    // main logger writes also to the trace file
    if ( m_env_vars->TraceActive )
    {
      logcpp::PatternLayout* t_layout = new logcpp::PatternLayout();
      std::string t_format("%m");
      t_layout->setConversionPattern(t_format);
      logcpp::Appender* t_appender = new logcpp::FileAppender("FileAppenderTrace",trace_filename.str(),append);
      t_appender->setLayout(t_layout);
      CFLogger::getInstance().getMainLogger().addAppender(t_appender);
    }
  }
  
  // TRACE TO FILES --------------------------------------------------------
  
  // create and set the trace log appender
  // always to cpu P0 and the others optionally
  // and only if tracing has been asked
  if ( m_env_vars->TraceActive )
  {
    if ( !m_env_vars->OnlyCPU0Writes || getCPURank() == 0 )
    {
      logcpp::PatternLayout* t_layout = new logcpp::PatternLayout();
//    std::string t_format("%d{%H:%M:%S %d/%m/%Y} %p %m");
      std::string t_format("[TRACE] %m");
      t_layout->setConversionPattern(t_format);

      logcpp::Appender* t_appender = new logcpp::FileAppender("TraceAppender",trace_filename.str(),append);
      t_appender->setLayout(t_layout);

      CFLogger::getInstance().getTraceLogger().addAppender(t_appender);
    }
    // sets default tracer priority
    CFLogger::getInstance().getTraceLogger().setPriority(DEBUG_MAX);
  }
}

//////////////////////////////////////////////////////////////////////////////

void CFEnv::initiate ( int argc, char** argv )
{
  // Initiate the Parallel environment
  // This modifies argc and argv! 
  
  COOLFluiD::Common::PE::InitPE(&argc, &argv);
  
  CFLog(VERBOSE, "CFEnv::initiate() => start\n");
  if (Common::PE::GetPE().GetRank("Default") == 0) {
    CFLog(INFO, "-------------------------------------------------------------\n");
    CFLog(INFO, "COOLFluiD Environment\n");
    CFLog(INFO, "-------------------------------------------------------------\n");
    CFLog(INFO, "COOLFluiD version    : " << getVersionString () << "\n");
    CFLog(INFO, "Parallel Environment : " << Common::PE::GetPE().GetName() << "\n");
    CFLog(INFO, "Build system         : " << getBuildSystem() << "\n");
    CFLog(INFO, "Build OS             : " << getLongSystemName() << " [" << getSystemBits() << "bits]\n");
    CFLog(INFO, "Build processor      : " << getBuildProcessor() << "\n");
  }
  
  m_env_vars->InitArgs.first  = argc;
  m_env_vars->InitArgs.second = argv;
  
  if (Common::PE::GetPE().GetRank("Default") == 0) {
    CFLog(VERBOSE, "Initializing Hook Modules ...\n");
  }
  
  initiateModules();
  
  if (Common::PE::GetPE().GetRank("Default") == 0) {
    CFLog(INFO, "-------------------------------------------------------------\n");
  }
  
  CFLog(VERBOSE, "CFEnv::initiate() => end\n");
}
    
//////////////////////////////////////////////////////////////////////////////

void CFEnv::initiateModules()
{
  std::vector< SafePtr<ModuleRegisterBase> > mod = m_moduleRegistry->getAllModules();
  std::for_each(mod.begin(),
                mod.end(),
                Common::safeptr_mem_fun(&ModuleRegisterBase::initiate));
}

//////////////////////////////////////////////////////////////////////////////

void CFEnv::terminateModules()
{
  std::vector< SafePtr<ModuleRegisterBase> > mod = m_moduleRegistry->getAllModules();
  std::for_each(mod.begin(),
                mod.end(),
                Common::safeptr_mem_fun(&ModuleRegisterBase::terminate));
}

//////////////////////////////////////////////////////////////////////////////

void CFEnv::terminate()
{
  CFLog(VERBOSE, "-------------------------------------------------------------\n");
  CFLog(VERBOSE, "COOLFluiD Environment Terminating\n");
  
  CFLog(VERBOSE, "Terminating Hook Modules ...\n");
  terminateModules();
  
  CFLog(VERBOSE, "Terminating Parallel Environment : Model " << Common::PE::GetPE().GetName() << "\n");
  Common::PE::DonePE ();
 
#if !defined(CF_HAVE_CRAYSTATIC) && defined(CF_HAVE_LOG4CPP) 
  CFLog(NOTICE, "-------------------------------------------------------------\n");
  CFLog(NOTICE, "COOLFluiD Environment Terminated\n");
  CFLog(NOTICE, "-------------------------------------------------------------\n");
#endif
}

//////////////////////////////////////////////////////////////////////////////

CFuint CFEnv::getCPURank()
{
  if (Common::PE::IsInitialised()) {
    return Common::PE::GetPE().GetRank("Default");
  }
  else {
    return 0;
  }
}

//////////////////////////////////////////////////////////////////////////////

Common::SafePtr<Environment::CFEnvVars> CFEnv::getVars()
{
  cf_assert(m_env_vars != CFNULL);
  return m_env_vars;
}

//////////////////////////////////////////////////////////////////////////////

Common::SafePtr<Environment::ModuleRegistry> CFEnv::getModuleRegistry()
{
  cf_assert(m_moduleRegistry != CFNULL);
  return m_moduleRegistry;
}

//////////////////////////////////////////////////////////////////////////////

Common::SafePtr<Common::FactoryRegistry> CFEnv::getFactoryRegistry()
{
  cf_assert(m_factoryRegistry != CFNULL);
  return m_factoryRegistry;
}

//////////////////////////////////////////////////////////////////////////////

Common::SafePtr<Common::EventHandler> CFEnv::getEventHandler()
{
  cf_assert(m_eventHandler.isNotNull());
  return m_eventHandler.getPtr();
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Environment

} // namespace COOLFluiD
