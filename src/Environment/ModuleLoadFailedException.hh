// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Environment_ModuleLoadFailedException_hh
#define COOLFluiD_Environment_ModuleLoadFailedException_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/Exception.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Environment {

//////////////////////////////////////////////////////////////////////////////

/// This class represents an Exception thrown when a certain
/// module was attempted to load but by some reason failed.
/// @author Tiago Quintino
class Environment_API ModuleLoadFailedException :
  public Common::Exception {
public:

  /// Constructor
  ModuleLoadFailedException ( const Common::CodeLocation& where, const std::string& what) :
    Common::Exception(where, what, "ModuleLoadFailedException") {}

  /// Copy constructor
  ModuleLoadFailedException ( const ModuleLoadFailedException& e) throw () : Exception(e) {}

}; // end of class ModuleLoadFailedException

//////////////////////////////////////////////////////////////////////////////

  } // namespace Environment

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Environment_ModuleLoadFailedException_hh
