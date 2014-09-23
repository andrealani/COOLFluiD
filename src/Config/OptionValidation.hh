// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Config_OptionValidation_hh
#define COOLFluiD_Config_OptionValidation_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/COOLFluiD.hh"
#include "Config/Config.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Config {

  class Option;

//////////////////////////////////////////////////////////////////////////////

/// This class is attached to a given option to validate its value according to
/// some restriction defined in the derived classes.
/// This class is an abstract interface to be implemented by derived classes.
/// @author Tiago Quintino
class Config_API OptionValidation {

public: // functions

  /// Default constructor without arguments.
  OptionValidation(Option * m_opt);

  /// Virtual destructor.
  virtual ~OptionValidation();

  /// Validate the option
  /// @throw OptionValidationException when value for option is not valid
  virtual bool isValid () = 0;

  /// Give the reason why validation fails
  /// @return string with the failing reason
  virtual std::string getReason () = 0;

protected: // data

  /// pointer to the concrete option to validate
  Option * m_opt;

}; // end of class OptionValidation

//////////////////////////////////////////////////////////////////////////////

  } // namespace Config

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Config_OptionValidation_hh
