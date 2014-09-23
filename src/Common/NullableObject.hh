// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Common_NullableObject_hh
#define COOLFluiD_Common_NullableObject_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/COOLFluiD.hh"
#include "Common/Common.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Common {

//////////////////////////////////////////////////////////////////////////////

/// This class represents an interface of objects
/// that have a Null counterpart.
/// @author Tiago Quintino
class Common_API NullableObject {
public:

  /// Checks if this object is a Null object.
  /// By default is always false, except for
  /// the concrete Null types that should implement
  /// it as returning true.
  /// @return true if Null and false otherwise
  virtual bool isNull() const;

  /// Checks if this object is a not Null object.
  /// @return true if not Null and false otherwise
  bool isNotNull() const {  return !isNull();  }

protected:

  /// Constructor.
  /// Protected constructor enforces this class to be non instantiable by clients.
  NullableObject();

  /// Default destructor
  /// Protected destructor enforces this class to be non deletable by clients.
  virtual ~NullableObject();

}; // end of class NullableObject

//////////////////////////////////////////////////////////////////////////////

  } // namespace Common

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

// #include "NullableObject.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Common_NullableObject_hh
