// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Common_CodeLocation_hh
#define COOLFluiD_Common_CodeLocation_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/COOLFluiD.hh"

namespace COOLFluiD {
  namespace Common {

//////////////////////////////////////////////////////////////////////////////

/// This class stores the information about a location in the source code
/// @author Tiago Quintino
class Common_API CodeLocation
{
public:

  /// constructor of the code location
  explicit CodeLocation (const char * file, int line, const char * function);

  /// @returns a string where the location is
  std::string str () const;

private:
  /// from which file the exception was thrown
  std::string m_file;
  /// from which function the exception was thrown
  /// @note will be empty if the compiler does not support it
  std::string m_function;
  /// from which line the exception was thrown
  int      m_line;

}; // end of class CodeLocation

//////////////////////////////////////////////////////////////////////////////

#define FromHere() COOLFluiD::Common::CodeLocation( __FILE__ , __LINE__ , __FUNCTION__ )

//////////////////////////////////////////////////////////////////////////////

  } // namespace Common
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Common_CodeLocation_hh
