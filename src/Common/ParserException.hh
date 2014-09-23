// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Common_ParserException_hh
#define COOLFluiD_Common_ParserException_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/Exception.hh"
#include "Common/Common.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Common {

//////////////////////////////////////////////////////////////////////////////

/// Exception thrown when user provides a bad value input
/// @author Tiago Quintino
class Common_API ParserException : public Common::Exception {
public:

  /// Constructor
  ParserException ( const Common::CodeLocation& where, const std::string& what);

  /// Copy constructor
  ParserException ( const ParserException& e) throw ();

}; // end of class ParserException

//////////////////////////////////////////////////////////////////////////////

  } // namespace Common

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Common_ParserException_hh
