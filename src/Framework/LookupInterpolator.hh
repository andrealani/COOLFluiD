// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_LookupInterpolator_hh
#define COOLFluiD_Framework_LookupInterpolator_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/StateInterpolator.hh"
#include "Common/OldLookupTable.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This class represents an interpolator based on lookup
/// @author Andrea Lani
class Framework_API LookupInterpolator : public StateInterpolator {

public: // functions

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);
  
  /// Default constructor without arguments
  LookupInterpolator(const std::string& name);
  
  /// Default destructor
  virtual ~LookupInterpolator();
  
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
  
 protected: // functions

  /// Fill in the lookup table
  virtual void fillTable();
  
 protected:  // data
 
  /// look up table for u(y)
  std::vector<Common::LookUpTable<CFreal,CFreal>*> m_lookupState;
  
  /// input data file name
  std::string m_infile;
  
}; // end of class LookupInterpolator

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

Framework_Factory(LookupInterpolator) // define the factory

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_LookupInterpolator_hh
