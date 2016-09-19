// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_StateInterpolator_hh
#define COOLFluiD_Framework_StateInterpolator_hh

//////////////////////////////////////////////////////////////////////////////

#include "MathTools/RealVector.hh"
#include "Environment/ConcreteProvider.hh"
#include "Framework/NumericalStrategy.hh"
#include "Framework/Storage.hh"
#include "Common/NullableObject.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a generic interface for interpolating state vectors
/// @author Andrea Lani
class Framework_API StateInterpolator : public Common::NullableObject, 
					public NumericalStrategy {

public: // typedefs

  typedef Environment::ConcreteProvider<StateInterpolator,1> PROVIDER;
  typedef const std::string& ARG1;

public: // functions

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Default constructor without arguments
  StateInterpolator(const std::string& name);

  /// Default destructor
  virtual ~StateInterpolator();
  
  /// Interpolate from given input scalar value
  /// @param varID  ID of the  output variable to interpolate
  /// @param in     input variable
  /// @param out    output variable
  virtual void interpolate(const CFuint varID, const CFreal& in, CFreal& out) = 0;
  
  /// Interpolate from given input state vector
  /// @param in     input state vector
  /// @param out    output state vector
  virtual void interpolate(const RealVector& in, RealVector& out) = 0;
  
  /// Setup the object
  virtual void setup();
  
  /// Unsetup the object
  virtual void unsetup();
  
  /// Returns the DataSocket's that this numerical strategy needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();
  
  /// Gets the Class name
  static std::string getClassName()
  {
    return "StateInterpolator";
  }

  /// Gets the polymorphic type name
  virtual std::string getPolymorphicTypeName() {return getClassName();}

protected:
  
}; // end of class StateInterpolator

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

Framework_Factory(StateInterpolator) // define the factory

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_StateInterpolator_hh
