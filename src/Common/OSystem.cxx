// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include <cstdlib>  // provides system call

#include "Common/ProcessInfo.hh"
#include "Common/SignalHandler.hh"
#include "Common/OSystem.hh"
#include "Common/StringOps.hh"

#ifndef CF_HAVE_ALLSTATIC
#ifdef CF_HAVE_DLOPEN
  #include "Common/PosixDlopenLibLoader.hh"
#endif
#endif

#ifdef CF_OS_LINUX
  #include "Common/ProcessInfoLinux.hh"
  #include "Common/SignalHandlerLinux.hh"
#endif

#ifdef CF_OS_MACOSX
  #include "Common/ProcessInfoMacOSX.hh"
  #include "Common/SignalHandlerMacOSX.hh"
#endif

#ifdef CF_OS_WINDOWS
  #include "Common/ProcessInfoWin32.hh"
  #include "Common/SignalHandlerWin32.hh"
  #include "Common/Win32LibLoader.hh"
#endif


//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;

namespace COOLFluiD {

  namespace Common {

//////////////////////////////////////////////////////////////////////////////

OSystem::OSystem() :
  m_process_info (CFNULL),
  m_sig_handler(CFNULL), 
  m_lib_loader(CFNULL)
{
#ifndef CF_HAVE_ALLSTATIC
#ifdef CF_HAVE_DLOPEN
    if ( m_lib_loader == CFNULL )   m_lib_loader = new PosixDlopenLibLoader();
#endif
#endif

#ifdef CF_OS_LINUX
    if ( m_process_info == CFNULL ) m_process_info = new ProcessInfoLinux();
    if ( m_sig_handler == CFNULL )  m_sig_handler = new SignalHandlerLinux();
#else
#ifdef CF_OS_MACOSX
    if ( m_process_info == CFNULL ) m_process_info = new ProcessInfoMacOSX();
    if ( m_sig_handler == CFNULL )  m_sig_handler = new SignalHandlerMacOSX();
#else
#ifdef CF_OS_WINDOWS
    if ( m_process_info == CFNULL ) m_process_info = new ProcessInfoWin32();
    if ( m_sig_handler == CFNULL )  m_sig_handler = new SignalHandlerWin32();
    if ( m_lib_loader == CFNULL )   m_lib_loader = new Win32LibLoader();
#else
  #error "Unkown operating system: not Windows, MacOSX or Linux"
#endif
#endif
#endif
}

//////////////////////////////////////////////////////////////////////////////

OSystem::~OSystem()
{
  deletePtr(m_process_info);
  deletePtr(m_sig_handler);
#ifndef CF_HAVE_ALLSTATIC
  deletePtr(m_lib_loader);
#endif
}

//////////////////////////////////////////////////////////////////////////////

OSystem& OSystem::getInstance()
{
  static OSystem osystem;
  return osystem;
}

//////////////////////////////////////////////////////////////////////////////

SafePtr<ProcessInfo> OSystem::getProcessInfo()
{
  return m_process_info;
}

//////////////////////////////////////////////////////////////////////////////

SafePtr<SignalHandler> OSystem::getSignalHandler()
{
  return m_sig_handler;
}

//////////////////////////////////////////////////////////////////////////////

SafePtr<LibLoader> OSystem::getLibLoader()
{
  cf_assert(m_lib_loader != CFNULL);
  return m_lib_loader;
}

//////////////////////////////////////////////////////////////////////////////

void OSystem::executeCommand(const std::string& call)
{
  int return_value = system ( call.c_str() );
  if ( return_value == -1)
  {
    std::string msg;
    msg += "Command \'";
    msg += call;
    msg += "\' return error code";
    throw OSystemException ( FromHere(), msg );
  }
}

//////////////////////////////////////////////////////////////////////////////

void OSystem::sleep (const CFuint& seconds)
{
  std::string callSleep = "sleep " + Common::StringOps::to_str(seconds);
  executeCommand (callSleep);
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Common

} // namespace COOLFluiD
