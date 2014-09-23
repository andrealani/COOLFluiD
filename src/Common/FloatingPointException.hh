// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Common_FloatingPointException_hh
#define COOLFluiD_Common_FloatingPointException_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/Exception.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Common {

//////////////////////////////////////////////////////////////////////////////

/// This class represents an Exception thrown when
/// a floating point error happens.
/// @author Tiago Quintino
class Common_API FloatingPointException : public Common::Exception {
public:

  /// Constructor
  FloatingPointException ( const Common::CodeLocation& where, const std::string& what) :
    Common::Exception(where, what, "FloatingPointException") {}

  /// Copy constructor
  FloatingPointException ( const FloatingPointException& e) throw () : Exception(e) {}

}; // end of class FloatingPointException

//////////////////////////////////////////////////////////////////////////////

  } // namespace Common

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Common_FloatingPointException_hh
