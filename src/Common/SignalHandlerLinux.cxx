// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include <cstdio>     // for printf()
#include <cstdlib>    // for free() and abort()
#include <csignal>    // POSIX signal(), SIGFPE and SIGSEGV
#include <fenv.h>     // floating Common access
#include <sstream>    // streamstring

#include "Common/Common.hh"
#include "Common/FloatingPointException.hh"
#include "Common/ProcessInfoLinux.hh"
#include "Common/SignalHandlerLinux.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;

namespace COOLFluiD {

  namespace Common {

//////////////////////////////////////////////////////////////////////////////

SignalHandlerLinux::SignalHandlerLinux()
{
}

//////////////////////////////////////////////////////////////////////////////

SignalHandlerLinux::~SignalHandlerLinux()
{
}

//////////////////////////////////////////////////////////////////////////////

void SignalHandlerLinux::registSignalHandlers()
{
  // register handler functions for the signals
  signal(SIGFPE,    (sighandler_t) SignalHandlerLinux::handleSIGFPE);
  signal(SIGSEGV,   (sighandler_t) SignalHandlerLinux::handleSIGSEGV);

  // enable the exceptions that will raise the SIGFPE signal
  feenableexcept ( FE_DIVBYZERO );
  feenableexcept ( FE_INVALID   );
  feenableexcept ( FE_OVERFLOW  );
  feenableexcept ( FE_UNDERFLOW );
}

//////////////////////////////////////////////////////////////////////////////

int SignalHandlerLinux::handleSIGFPE (int signal)
{
  printf("\nreceived signal SIGFPE [%d] - 'Floating Point Exception'\n",signal);
  static std::string dump = ProcessInfoLinux::dumpBacktrace();
  printf( "%s\n", dump.c_str() );
  throw Common::FloatingPointException (FromHere(), "Some floating point operation has given an invalid result");
}

//////////////////////////////////////////////////////////////////////////////

int SignalHandlerLinux::handleSIGSEGV(int signal)
{
  printf("\nreceived signal SIGSEGV [%d] - 'Segmentation violation'\n",signal);
  static std::string dump = ProcessInfoLinux::dumpBacktrace();
  printf( "%s\n", dump.c_str() );
  abort();
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Common

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
