// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Common_BadValueException_hh
#define COOLFluiD_Common_BadValueException_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/Exception.hh"
#include "Common/Common.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Common {

//////////////////////////////////////////////////////////////////////////////

/// Exception thrown when user provides a bad value input
/// @author Tiago Quintino
class Common_API BadValueException : public Common::Exception {
public:

  /// Constructor
  BadValueException ( const Common::CodeLocation& where, const std::string& what) :
    Common::Exception(where, what, "BadValueException") {}

  /// Copy constructor
  BadValueException ( const BadValueException& e) throw () : Exception(e) {}

}; // end of class BadValueException

//////////////////////////////////////////////////////////////////////////////

  } // namespace Common

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Common_BadValueException_hh
