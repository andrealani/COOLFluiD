// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include <cstdlib>

#include "Common/CFLog.hh" // CF_HAVE_WINDOWSH is defined via this header

// windows header ( maybe this define should be in the coolfluid_config.h )
#define _WIN32_WINNT 0x0502 // minimum requirement is ( "Windows Server 2003 with SP1" ) or ( "Windows XP with SP2" )
#include <windows.h>

#include "Common/Win32LibLoader.hh"
#include "Common/Common.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Common {

//////////////////////////////////////////////////////////////////////////////

Win32LibLoader::Win32LibLoader()
{
}

//////////////////////////////////////////////////////////////////////////////

Win32LibLoader::~Win32LibLoader()
{
}

//////////////////////////////////////////////////////////////////////////////

void Win32LibLoader::set_search_paths(std::vector< boost::filesystem::path >& paths)
{
  m_search_paths = paths;
}

//////////////////////////////////////////////////////////////////////////////

void Win32LibLoader::load_library(const std::string& lib)
{
  CFLog( VERBOSE , "Win32LibLoader: Attempting to load '" << lib << "'\n" );
  using boost::filesystem::path;
  std::string libname = lib + ".dll";

  // attempt to load a module
  HINSTANCE hdl = CFNULL;
  std::vector<path>::const_iterator itr = m_search_paths.begin();
  for (; itr != m_search_paths.end() ; ++itr)
  {
    // set the current search directory
    CFLog ( WARN, "searching in dir: '" <<  (*itr).string() << "'\n" );
	  SetDllDirectory((*itr).string().c_str());

	  // try to load the library
    hdl = LoadLibrary(TEXT(libname.c_str()));
    if( hdl != CFNULL ) break;
    CFLog ( WARN, "didnt find '" <<  lib << "' in dir '" <<  (*itr).string() << "'\n" );
  }

  // searhc in the global path
  hdl = LoadLibrary(TEXT(libname.c_str()));
  if( hdl == CFNULL )
  {
    char * pPath;
    pPath = getenv ("PATH");
    if (pPath != CFNULL)
      CFLog ( WARN, "didnt find '" <<  lib << "' in global path variable PATH='" << pPath << "'\n" );
  }

  // check for success
  if(hdl != NULL)
  {
    CFLog( WARN, "Win32LibLoader: Loaded '" << libname << "'\n" );
  }
  else
  {
    CFLog( WARN, "LoadLibrary() failed to load module : '" << libname << "'\n" );
    throw LibLoaderException (FromHere(), "Module failed to load '" + libname + "'"  );
  }
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Common

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

