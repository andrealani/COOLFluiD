// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Common_URLException_hh
#define COOLFluiD_Common_URLException_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/Exception.hh"

#include "Common/Common.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Common {

//////////////////////////////////////////////////////////////////////////////

/// This class represents the exception thrown
/// if an error occurs when accessing a network URL.
/// @author Tiago Quintino
class Common_API URLException : public Common::Exception {
public:

  /// Constructor
  URLException (const Common::CodeLocation& where, const std::string& what) :
    Exception(where, what, "URLException") {}

  /// Copy constructor
  URLException(const URLException& e) throw () : Exception(e) {}

}; // end of class URLException

//////////////////////////////////////////////////////////////////////////////

  } // namespace Common

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Common_URLException_hh
