// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include <sstream>
#include <iostream>

#include "Common/Common.hh"
#include "Common/CFAssert.hh"
#include "Common/FailedAssertionException.hh"

#include "Common/OSystem.hh"
#include "Common/ProcessInfo.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;

namespace COOLFluiD {

//////////////////////////////////////////////////////////////////////////////

AssertionManager::AssertionManager() :
  DoAssertions    ( true ),
  AssertionDumps  ( true ),
  AssertionThrows ( true ) {}

//////////////////////////////////////////////////////////////////////////////

AssertionManager& AssertionManager::getInstance()
{
  static AssertionManager assertion_manager;
  return assertion_manager;
}

//////////////////////////////////////////////////////////////////////////////

void AssertionManager::do_assert ( bool condition,
                                   const char * cond_str,
                                   const char * file,
                                   int line,
                                   const char * func,
                                   const char * desc )
{
  if ( (!condition) && AssertionManager::getInstance().DoAssertions )
  {
    std::ostringstream out;
    out << "Assertion failed: [" << cond_str << "] ";

    if (desc)
      out << "'" << desc << "' ";

    out << "in " << file << ":" << line;

    if (func)
      out << " [function " << func << "]";

    if ( AssertionManager::getInstance().AssertionDumps )
      out << "\n" << OSystem::getInstance().getProcessInfo()->getBackTrace();

    if ( AssertionManager::getInstance().AssertionThrows )
    {
      throw FailedAssertionException (FromHere(),out.str());
    }
    else
    {
      std::cerr << out.str() << std::endl;
      cerr.flush ();
      abort ();
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

} // namespace COOLFluiD
