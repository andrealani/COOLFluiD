// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Common_FailedAssertionException_hh
#define COOLFluiD_Common_FailedAssertionException_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/Exception.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Common {

//////////////////////////////////////////////////////////////////////////////

/// This class represents an Exception thrown  when an assertion fails but
/// the code is configured to throw an exception rather than crash.
/// @author Tiago Quintino
class Common_API FailedAssertionException : public Common::Exception {
public:

  /// Constructor
  FailedAssertionException (const Common::CodeLocation& where, const std::string& what) :
    Common::Exception(where, what, "FailedAssertionException") {}

  /// Copy constructor
  FailedAssertionException(const FailedAssertionException& e) throw  () : Exception(e) {}

}; // end of class FailedAssertionException

//////////////////////////////////////////////////////////////////////////////

  } // namespace Common

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Common_FailedAssertionException_hh
