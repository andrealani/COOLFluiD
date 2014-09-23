// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Config_OptionValidationException_hh
#define COOLFluiD_Config_OptionValidationException_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/Exception.hh"
#include "Config/Config.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Config {

//////////////////////////////////////////////////////////////////////////////

/// Exception raised when an option validatio failed
/// @author Tiago Quintino
class Config_API OptionValidationException : public Common::Exception {
public:

  /// Constructor
  OptionValidationException ( const Common::CodeLocation& where, const std::string& what) :
    Common::Exception(where, what, "OptionValidationException") {}

  /// Copy constructor
  OptionValidationException ( const OptionValidationException& e) throw () : Exception(e) {}

}; // end of class OptionValidationException

//////////////////////////////////////////////////////////////////////////////

  } // namespace Config

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Config_OptionValidationException_hh
