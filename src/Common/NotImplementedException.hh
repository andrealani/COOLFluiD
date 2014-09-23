// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Common_NotImplementedException_hh
#define COOLFluiD_Common_NotImplementedException_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/Exception.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Common {

//////////////////////////////////////////////////////////////////////////////

/// This class represents an Exception throwna certain functionality is not implemented
/// @author Andrea Lani
/// @author Tiago Quintino
class Common_API  NotImplementedException : public Common::Exception {
public:

  /// Constructor
  /// @see COOLFluiD::Exception()
  NotImplementedException(const Common::CodeLocation& where, const std::string& what) :
    Exception(where, what,"NotImplementedException") {}

  /// A copy constructor is necessary for exceptions, for the C++
  /// exception mechanism to work.
  NotImplementedException(const NotImplementedException& e) throw() : Exception(e) {}

}; // end of class NotImplementedException

//////////////////////////////////////////////////////////////////////////////

  } // namespace Common

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Common_NotImplementedException_hh
