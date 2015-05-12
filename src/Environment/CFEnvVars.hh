// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOFluiD_Environment_CFEnvVars_hh
#define COOFluiD_Environment_CFEnvVars_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/COOLFluiD.hh"
#include "Common/NonCopyable.hh"
#include "Common/CFLog.hh"
#include "Environment/EnvironmentAPI.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Environment {

//////////////////////////////////////////////////////////////////////////////

/// Stores COOLFluiD Runtime environment variables
/// @author Tiago Quintino
class Environment_API CFEnvVars : public Common::NonCopyable<CFEnvVars>
{
  public: // functions

    /// Constructor
    CFEnvVars();

  public: // data

    /// only processor P0 outputs the log info to files
    bool OnlyCPU0Writes;
    /// assertions throw exceptions instead of aborting
    bool AssertionThrows;
    /// regist signal handlers
    bool RegistSignalHandlers;
    /// activate trace
    bool TraceActive;
    /// tracing also sento to StdOut to be put into CFEnv
    bool TraceToStdOut;
    /// If Events have verbose output
    bool VerboseEvents;
    /// Signal error when some user provided config parameters are not used
    bool ErrorOnUnusedConfig;
    /// the name of the file in which to put the logging messages
    std::string MainLoggerFileName;
    /// the loglevel for exceptions
    CFuint ExceptionLogLevel;
    /// the initial arguments with which the environment was started
    std::pair<int,char**> InitArgs;
    /// number of writing processes in parallel I/O
    CFuint NbWriters;
    
}; // end class CFEnvVars

//////////////////////////////////////////////////////////////////////////////

  } // namespace Environment

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOFluiD_Environment_CFEnvVars_hh
