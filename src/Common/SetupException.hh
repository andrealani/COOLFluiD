// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Common_SetupException_hh
#define COOLFluiD_Common_SetupException_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/Exception.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Common {

//////////////////////////////////////////////////////////////////////////////

/// This class represents the exception thrown if an error occurs when setting up an object.
/// @author Tiago Quintino
class SetupException : public Common::Exception {
public:

  /// Constructor
  SetupException ( const Common::CodeLocation& where, const std::string& what) :
    Common::Exception(where, what, "SetupException") {}

  /// Copy constructor
  SetupException ( const SetupException& e) throw () : Exception(e) {}

}; // end of class SetupException

//////////////////////////////////////////////////////////////////////////////

  } // namespace Common

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Common_SetupException_hh
