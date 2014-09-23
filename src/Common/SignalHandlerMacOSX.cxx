// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Common/FloatingPointException.hh"
#include "Common/Common.hh"
#include "Common/SignalHandlerMacOSX.hh"
#include "Common/ProcessInfoMacOSX.hh"

#include <cstdio>     // for printf()
#include <cstdlib>    // for free() and abort()
#include <csignal>    // POSIX signal(), SIGFPE and SIGSEGV
#include <fenv.h>     // floating Common access
#include <sstream>    // streamstring

//////////////////////////////////////////////////////////////////////////////

using namespace std;

namespace COOLFluiD {

  namespace Common {
    
//////////////////////////////////////////////////////////////////////////////

/*
 * Following functions are required since they are not available for Mac OSX
 * This only works for intel architecture
 * http://www-personal.umich.edu/~williams/archive/computation/fe-handling-example.c
 */ 

static int
fegetexcept (void)
{
  static fenv_t fenv;

  return fegetenv (&fenv) ? -1 : (fenv.__control & FE_ALL_EXCEPT);
}

static int
feenableexcept (unsigned int excepts)
{
  static fenv_t fenv;
  unsigned int new_excepts = excepts & FE_ALL_EXCEPT,
               old_excepts;  // previous masks

  if ( fegetenv (&fenv) ) return -1;
  old_excepts = fenv.__control & FE_ALL_EXCEPT;

  // unmask
  fenv.__control &= ~new_excepts;
  fenv.__mxcsr   &= ~(new_excepts << 7);

  return ( fesetenv (&fenv) ? -1 : old_excepts );
}

static int
fedisableexcept (unsigned int excepts)
{
  static fenv_t fenv;
  unsigned int new_excepts = excepts & FE_ALL_EXCEPT,
               old_excepts;  // all previous masks

  if ( fegetenv (&fenv) ) return -1;
  old_excepts = fenv.__control & FE_ALL_EXCEPT;

  // mask
  fenv.__control |= new_excepts;
  fenv.__mxcsr   |= new_excepts << 7;

  return ( fesetenv (&fenv) ? -1 : old_excepts );
}

//////////////////////////////////////////////////////////////////////////////

SignalHandlerMacOSX::SignalHandlerMacOSX()
{
}

//////////////////////////////////////////////////////////////////////////////

SignalHandlerMacOSX::~SignalHandlerMacOSX()
{
}

//////////////////////////////////////////////////////////////////////////////

void SignalHandlerMacOSX::registSignalHandlers()
{
  // register handler functions for the signals
  signal(SIGFPE,    (sighandler_t) SignalHandlerMacOSX::handleSIGFPE);
  signal(SIGSEGV,   (sighandler_t) SignalHandlerMacOSX::handleSIGSEGV);

  // enable the exceptions that will raise the SIGFPE signal
  feenableexcept ( FE_DIVBYZERO );
  feenableexcept ( FE_INVALID   );
  feenableexcept ( FE_OVERFLOW  );
  feenableexcept ( FE_UNDERFLOW );
}

//////////////////////////////////////////////////////////////////////////////

int SignalHandlerMacOSX::handleSIGFPE (int signal)
{
  printf("\nreceived signal SIGFPE [%d] - 'Floating Point Exception'\n",signal);
  static std::string dump = ProcessInfoMacOSX::dumpBackTrace();
  printf( "%s\n", dump.c_str() );
  throw Common::FloatingPointException (FromHere(), "Some floating point operation has given an invalid result");
}

//////////////////////////////////////////////////////////////////////////////

int SignalHandlerMacOSX::handleSIGSEGV(int signal)
{
  printf("\nreceived signal SIGSEGV [%d] - 'Segmentation violation'\n",signal);
  static std::string dump = ProcessInfoMacOSX::dumpBackTrace();
  printf( "%s\n", dump.c_str() );
  abort();
}


//////////////////////////////////////////////////////////////////////////////

  } // namespace Common

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
