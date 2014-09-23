// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Common_FailedCastException_hh
#define COOLFluiD_Common_FailedCastException_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/Exception.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Common {

//////////////////////////////////////////////////////////////////////////////

/// This class represents an Exception thrown
/// when a dynamic cast of a pointer fails.
/// @author Tiago Quintino
class FailedCastException : public Common::Exception {
public:

  /// Constructor
  FailedCastException ( const Common::CodeLocation& where, const std::string& what) :
    Common::Exception(where, what, "FailedCastException") {}

  /// Copy constructor
  FailedCastException ( const FailedCastException& e) throw () : Exception(e) {}

}; // end of class FailedCastException

//////////////////////////////////////////////////////////////////////////////

  } // namespace Common

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Common_FailedCastException_hh
