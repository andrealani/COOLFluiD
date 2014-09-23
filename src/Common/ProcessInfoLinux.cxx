// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include <execinfo.h>    // for backtrace() from glibc
#include <sys/types.h>   // for getting the PID of the process

#ifdef CF_HAVE_UNISTD_H
  #include <unistd.h>
#endif

#include <malloc.h>      //  for mallinfo

#include <cstdio>
#include <sstream>       // streamstring
#include <cstdlib>       // for free() and abort()

#include "Common/ProcessInfoLinux.hh"
#include "Common/Common.hh"
#include "Common/CFLog.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;

namespace COOLFluiD {

  namespace Common {

//////////////////////////////////////////////////////////////////////////////

ProcessInfoLinux::ProcessInfoLinux()
{
}

//////////////////////////////////////////////////////////////////////////////

ProcessInfoLinux::~ProcessInfoLinux()
{
}

//////////////////////////////////////////////////////////////////////////////

std::string ProcessInfoLinux::getBackTrace () const
{
  return dumpBacktrace ();
}

//////////////////////////////////////////////////////////////////////////////

std::string ProcessInfoLinux::dumpBacktrace ()
{
#define CF_BUFFER_SIZE 256

  printf ("\n\ndumping backtrace ...\n");

  std::ostringstream oss;
  int j, nptrs;
  void *buffer[CF_BUFFER_SIZE];
  char **strings;

  nptrs = backtrace(buffer, CF_BUFFER_SIZE);
  oss << "\nbacktrace() returned " << nptrs << " addresses\n";

  strings = backtrace_symbols(buffer, nptrs);
  if (strings == NULL)
    oss << "\nno backtrace_symbols found\n";
  for (j = 0; j < nptrs; j++)
    oss << strings[j] << "\n";
  free(strings);

#undef CF_BUFFER_SIZE

  return oss.str();
}

//////////////////////////////////////////////////////////////////////////////

CFuint ProcessInfoLinux::getPID() const
{
  pid_t pid = getpid();
  return static_cast<CFuint> ( pid );
}

//////////////////////////////////////////////////////////////////////////////

CFdouble ProcessInfoLinux::memoryUsageBytes() const
{
  struct mallinfo info;

  // get current memory
  info = mallinfo();
#if 0
  std::cout
  << "### MALLINFO ###"  << std::endl
  << "  number of free chunks " <<  info.ordblks << std::endl
  << "  number of fastbin blocks " << info.smblks << std::endl
  << "  number of mmapped regions " << info.hblks << std::endl
  << "  non-mmapped space allocated from system " << info.arena /1024 << "KB" << std::endl
  << "  space in mmapped regions " << info.hblkhd /1024 << "KB" << std::endl
  << "  maximum total allocated space " << info.usmblks/1024 << "KB"<< std::endl
  << "  space available in freed fastbin blocks " << info.fsmblks/1024 << "KB"<< std::endl
  << "  total allocated space " << info.uordblks/1024 << "KB"<< std::endl
  << "  total free space " << info.fordblks/1024 << "KB" << std::endl
  << "  top-most, releasable (via malloc_trim) space " << info.keepcost/1024 << "KB"<< std::endl
  << "################" << std::endl
  << "SUM BEFORE : " << ( info.arena + info.hblkhd ) / 1024  << "KB" << std::endl
  << "SUM ALL    : " << ( info.arena + info.hblkhd + info.fsmblks + info.uordblks ) / 1024 << "KB" << std::endl
  << "################" << std::endl;

  char buf[30];
  snprintf(buf, 30, "/proc/%u/statm", (unsigned)getpid());
        FILE* pf = fopen(buf, "r");
        if (pf) {
            unsigned size; //       total program size
            unsigned resident;//   resident set size
            unsigned share;//      shared pages
            unsigned text;//       text (code)
            unsigned lib;//        library
            unsigned data;//       data/stack
            unsigned dt;//         dirty pages (unused in Linux 2.6)
            fscanf(pf, "%u %u %u %u %u %u", &size, &resident, &share, &text, &lib, &data);
            CFLog(INFO, "NEW METHOD size" << size / (1024.0) << " MB mem used\n");
            CFLog(INFO, "NEW METHOD resident" << resident / (1024.0) << " MB mem used\n");
            CFLog(INFO, "NEW METHOD share" << share / (1024.0) << " MB mem used\n");
            CFLog(INFO, "NEW METHOD text" << text / (1024.0) << " MB mem used\n");
            CFLog(INFO, "NEW METHOD lib" << lib / (1024.0) << " MB mem used\n");
            CFLog(INFO, "NEW METHOD data" << data / (1024.0) << " MB mem used\n");
        }
#endif

  return static_cast<CFdouble>(info.arena) +
         static_cast<CFdouble>(info.hblkhd);
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Common

} // namespace COOLFluiD

