// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include <execinfo.h>    // for backtrace() from glibc

#include <unistd.h>      // for getting the PID of the process
#include <sys/types.h>   // for getting the PID of the process

#include <mach/mach_types.h> 
#include <mach/mach_init.h>
#include <mach/task.h>

#include "Common/ProcessInfoMacOSX.hh"
#include "Common/Common.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;

namespace COOLFluiD {

  namespace Common {

//////////////////////////////////////////////////////////////////////////////

ProcessInfoMacOSX::ProcessInfoMacOSX()
{
}

//////////////////////////////////////////////////////////////////////////////

ProcessInfoMacOSX::~ProcessInfoMacOSX()
{
}

//////////////////////////////////////////////////////////////////////////////

std::string ProcessInfoMacOSX::getBackTrace () const
{
  return dumpBackTrace ();
}

//////////////////////////////////////////////////////////////////////////////

std::string ProcessInfoMacOSX::dumpBackTrace ()
{
  #define BUFFER_SIZE 500

  printf ("\n\ndumping backtrace ...\n");

  std::ostringstream oss;
  int j, nptrs;
  void *buffer[BUFFER_SIZE];
  char **strings;

  nptrs = backtrace(buffer, BUFFER_SIZE);
  oss << "\nbacktrace() returned " << nptrs << " addresses\n";

  strings = backtrace_symbols(buffer, nptrs);
  if (strings == NULL)
    oss << "\nno backtrace_symbols found\n";
  for (j = 0; j < nptrs; j++)
    oss << strings[j] << "\n";
  free(strings);

  #undef BUFFER_SIZE

  return oss.str();
}

//////////////////////////////////////////////////////////////////////////////

CFuint ProcessInfoMacOSX::getPID() const
{
  pid_t pid = getpid();
  return static_cast<CFuint> ( pid );
}

//////////////////////////////////////////////////////////////////////////////

CFdouble ProcessInfoMacOSX::memoryUsageBytes() const
{
  
  struct task_basic_info t_info;
  mach_msg_type_number_t t_info_count = TASK_BASIC_INFO_COUNT;

  if (KERN_SUCCESS != task_info(mach_task_self(),
                                TASK_BASIC_INFO, (task_info_t)&t_info, 
                                &t_info_count))
  {
      return -1;
  }
  // resident size is in t_info.resident_size;
  // virtual size is in t_info.virtual_size;
  
  return static_cast<CFdouble>(t_info.resident_size);
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Common

} // namespace COOLFluiD

