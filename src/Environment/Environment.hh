// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Environment_hh
#define COOLFluiD_Environment_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/ModuleRegister.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  /// Basic classes that form the COOLFluiD Environment that bind together the external components
  namespace Environment {

//////////////////////////////////////////////////////////////////////////////

/// This class defines the Module Environment
class Environment_API EnvironmentModule :
    public Environment::ModuleRegister<EnvironmentModule> {

public:

  /// Static function that returns the module name.
  /// Must be implemented for the ModuleRegister template
  /// @return name of the module
  static std::string getModuleName()  {  return "Environment";  }

  /// Static function that returns the description of the module.
  /// Must be implemented for the ModuleRegister template
  /// @return descripton of the module
  static std::string getModuleDescription() {  return "This module implementes the environment that binds together the external components.";  }

}; // end EnvironmentModule

//////////////////////////////////////////////////////////////////////////////

  } // namespace Environment

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Environment_hh
