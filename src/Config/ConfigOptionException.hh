// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Config_ConfigOptionException_hh
#define COOLFluiD_Config_ConfigOptionException_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/Exception.hh"
#include "Config/Config.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Config {

//////////////////////////////////////////////////////////////////////////////

/// Exception raised when config options are set incorrectly.
/// @author Tiago Quintino
class Config_API ConfigOptionException : public Common::Exception {
public:

  /// Constructor
  ConfigOptionException ( const Common::CodeLocation& where, const std::string& what) :
    Common::Exception(where, what, "ConfigOptionException") {}

  /// Copy constructor
  ConfigOptionException ( const ConfigOptionException& e) throw () : Exception(e) {}

}; // end of class ConfigOptionException

//////////////////////////////////////////////////////////////////////////////

  } // namespace Config

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Config_ConfigOptionException_hh
