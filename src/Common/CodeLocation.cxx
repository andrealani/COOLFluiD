// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include <cstdio>

#include "Common/CodeLocation.hh"

//////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Common {

//////////////////////////////////////////////////////////////////////

CodeLocation::CodeLocation(const char * file, int line, const char * function)
 : m_file(file), m_function(function), m_line (line)
{
}

//////////////////////////////////////////////////////////////////////

std::string CodeLocation::str () const
{
  char line [50];
  sprintf (line, "%d", m_line);
  std::string place (m_file);
  place += ":";
  place += line;
  if (!m_function.empty()) // skip if compiler doees not set function
  {
    place += ":";
    place += m_function;
  }
  return place;
}

//////////////////////////////////////////////////////////////////////

  } // namespace Common

} // namespace COOLFluiD

