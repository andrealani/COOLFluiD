// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_NullStateInterpolator_hh
#define COOLFluiD_Framework_NullStateInterpolator_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/StateInterpolator.hh"
#include "Common/OldLookupTable.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This class represents an interpolator based on lookup
/// @author Andrea Lani
class Framework_API NullStateInterpolator : public StateInterpolator {

public: // functions

  /// Default constructor without arguments
  NullStateInterpolator(const std::string& name);
  
  /// Default destructor
  virtual ~NullStateInterpolator();
  
  /// Interpolate from given input scalar value
  /// @param varID  ID of the  output variable to interpolate
  /// @param in     input variable
  /// @param out    output variable
  virtual void interpolate(const CFuint varID, const CFreal& in, CFreal& out);
  
  /// Interpolate from given input state vector
  /// @param in     input state vector
  /// @param out    output state vector
  virtual void interpolate(const RealVector& in, RealVector& out);
  
  /// Setup the object
  virtual void setup();
  
  /// Unsetup the object
  virtual void unsetup();
  
  /// Checks if this object is a Null object.
  /// By default is always false, except for
  /// the concrete Null types that should implement
  /// it as returning true.
  /// @return true if Null and false otherwise
  virtual bool isNull() const {return true;}
  
}; // end of class NullStateInterpolator

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

Framework_Factory(NullStateInterpolator) // define the factory

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_NullStateInterpolator_hh
