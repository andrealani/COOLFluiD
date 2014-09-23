// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Config_NoEnvVarException_hh
#define COOLFluiD_Config_NoEnvVarException_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/Exception.hh"
#include "Config/Config.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Config {

//////////////////////////////////////////////////////////////////////////////

/// This class represents exception that arise
/// when the searched Enviromental Variable does not exist.
/// @author Tiago Quintino
class Config_API NoEnvVarException : public Common::Exception {
public:

  /// Constructor
  NoEnvVarException ( const Common::CodeLocation& where, const std::string& what) :
    Common::Exception(where, what, "NoEnvVarException") {}

  /// Copy constructor
  NoEnvVarException ( const NoEnvVarException& e) throw () : Exception(e) {}

}; // end of class NoEnvVarException

//////////////////////////////////////////////////////////////////////////////

  } // namespace Config
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Config_NoEnvVarException_hh
