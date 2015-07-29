// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_NullPhysicalPropertyLibrary_hh
#define COOLFluiD_Framework_NullPhysicalPropertyLibrary_hh

//////////////////////////////////////////////////////////////////////////////

#include "PhysicalPropertyLibrary.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// Provides an abstract interface for libraries that compute the physical
/// properties.
/// @author Andrea Lani
class Framework_API NullPhysicalPropertyLibrary : public PhysicalPropertyLibrary {
public:

  /// Constructor
  NullPhysicalPropertyLibrary(const std::string& name);

  /// Default destructor
  ~NullPhysicalPropertyLibrary();

  /// Setups the data of the library
  void setup();

  /// Unsetups the data of the library
  void unsetup();
  
  /// Checks if this object is a Null object.
  /// By default is always false, except for
  /// the concrete Null types that should implement
  /// it as returning true.
  /// @return true if Null and false otherwise
  virtual bool isNull() const {return true;}
  
}; // end of class NullPhysicalPropertyLibrary

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_NullPhysicalPropertyLibrary_hh
