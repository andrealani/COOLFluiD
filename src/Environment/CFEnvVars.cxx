// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Environment/CFEnvVars.hh"
#include "Environment/ModuleRegisterBase.hh"
#include "Common/PE.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Common;

namespace COOLFluiD {

  namespace Environment {

//////////////////////////////////////////////////////////////////////////////

CFEnvVars::CFEnvVars() :
  OnlyCPU0Writes       ( true  ),
  RegistSignalHandlers ( false ),
  TraceActive          ( false ),
  TraceToStdOut        ( false ),
  VerboseEvents        ( false ),
  ErrorOnUnusedConfig  ( false ),
  MainLoggerFileName("output.log"),
  SyncAlgo("AllToAll"),
  ExceptionLogLevel( (CFuint) VERBOSE),
  InitArgs()
{
  InitArgs.first  = 0;
  InitArgs.second = CFNULL;
  NbWriters = 1;
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Environment

} // namespace COOLFluiD
