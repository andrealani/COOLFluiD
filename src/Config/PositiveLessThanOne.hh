// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Config_PositiveLessThanOne_hh
#define COOLFluiD_Config_PositiveLessThanOne_hh

//////////////////////////////////////////////////////////////////////////////

#include "Config/OptionValidation.hh"
#include "Config/Option.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Config {

//////////////////////////////////////////////////////////////////////////////

/// This class validates the option value to be bigger than zero
/// @author Tiago Quintino
class Config_API PositiveLessThanOne : public OptionValidation {

public: // functions

  /// Default constructor without arguments.
  explicit PositiveLessThanOne(Option * opt);

  /// Virtual destructor.
  ~PositiveLessThanOne();

  /// Validate the option
  virtual bool isValid ();

  /// Give the reason why validation fails
  /// @return string with the failing reason
  virtual std::string getReason ();

}; // end of class PositiveLessThanOne

//////////////////////////////////////////////////////////////////////////////

  } // namespace Config

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Config_PositiveLessThanOne_hh
