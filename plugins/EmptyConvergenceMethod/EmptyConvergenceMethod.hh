// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_EmptyConvergenceMethod_hh
#define COOLFluiD_EmptyConvergenceMethod_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/ModuleRegister.hh"
#include "EmptyConvergenceMethod/EmptyConvergenceMethodAPI.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  /// The classes that implement EmptyConvergenceMethod time stepper.
  namespace EmptyConvergenceMethod {

//////////////////////////////////////////////////////////////////////////////

/// This class defines the Module EmptyConvergenceMethod
class  EmptyConvergenceMethod_API EmptyConvergenceMethodLib : 
      public Environment::ModuleRegister<EmptyConvergenceMethodLib> {

public:
  
  /// Static function that returns the module name.
  /// Must be implemented for the ModuleRegister template
  /// @return name of the module
  static std::string getModuleName() { return "EmptyConvergenceMethod"; }
  
  /// Static function that returns the description of the module.
  /// Must be implemented for the ModuleRegister template
  /// @return descripton of the module
  static std::string getModuleDescription()
  {
    return "This module implements a empty iterator method.";
  }
  
}; // end EmptyConvergenceMethodLib

//////////////////////////////////////////////////////////////////////////////

  } // namespace EmptyConvergenceMethod

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFLUID_EmptyConvergenceMethod_hh

