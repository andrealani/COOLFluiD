// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Numerics_ForwardEuler_hh
#define COOLFluiD_Numerics_ForwardEuler_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/ModuleRegister.hh"
#include "ForwardEuler/ForwardEulerAPI.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



  /// The classes that implement ForwardEuler time stepper.
  namespace ForwardEuler {

//////////////////////////////////////////////////////////////////////////////

/// This class defines the Module ForwardEuler
class  ForwardEuler_API ForwardEulerLib : public Environment::ModuleRegister<ForwardEulerLib> {
public:

  /// Static function that returns the module name.
  /// Must be implemented for the ModuleRegister template
  /// @return name of the module
  static std::string getModuleName() { return "ForwardEuler"; }

  /// Static function that returns the description of the module.
  /// Must be implemented for the ModuleRegister template
  /// @return descripton of the module
  static std::string getModuleDescription()
  {
    return "This module implements a forward Euler time stepper.";
  }

}; // end ForwardEulerLib

//////////////////////////////////////////////////////////////////////////////

  } // namespace ForwardEuler

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFLUID_Numerics_ForwardEuler_hh

