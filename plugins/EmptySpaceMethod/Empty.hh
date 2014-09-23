// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_EmptySpaceMethod_EmptyModule_hh
#define COOLFluiD_EmptySpaceMethod_EmptyModule_hh

#include "Environment/ModuleRegister.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace EmptySpaceMethod {

//////////////////////////////////////////////////////////////////////////////

/// This class defines an empty module
class EmptyModule : public Environment::ModuleRegister< EmptyModule > {
public:

  /// Static function that returns the module name. Must be implemented for the
  /// ModuleRegister template
  /// @return name of the module
  static std::string getModuleName() {
    return "Empty";
  }

  /// Static function that returns the description of the module. Must be
  /// implemented for the ModuleRegister template
  /// @return descripton of the module
  static std::string getModuleDescription() {
    return "This module interfaces and empty solver.";
  }

}; // end EmptyModule

//////////////////////////////////////////////////////////////////////////////

  }  // namespace EmptySpaceMethod
}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_EmptySpcaeMethod_EmptyModule_hh

