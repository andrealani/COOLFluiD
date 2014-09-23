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

}; // end of class NullPhysicalPropertyLibrary

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_NullPhysicalPropertyLibrary_hh
