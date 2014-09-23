// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Physics_LinearAdv_hh
#define COOLFluiD_Physics_LinearAdv_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/ModuleRegister.hh"
#include "Common/ExportAPI.hh"

//////////////////////////////////////////////////////////////////////////////

/// Define the macro LinearAdv_API
/// @note build system defines LinearAdv_EXPORTS when compiling FluctSplit files
#ifdef LinearAdv_EXPORTS
#   define LinearAdv_API CF_EXPORT_API
#else
#   define LinearAdv_API CF_IMPORT_API
#endif

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {
  /// The classes that implement Linear Advection physical model.
  namespace LinearAdv {

//////////////////////////////////////////////////////////////////////////////

/// This class defines the Module LinearAdv
class LinearAdv_API LinearAdvModule : public Environment::ModuleRegister<LinearAdvModule> {
public:

  /// Static function that returns the module name.
  /// Must be implemented for the ModuleRegister template
  /// @return name of the module
  static std::string getModuleName() { return "LinearAdv"; }

  /// Static function that returns the description of the module.
  /// Must be implemented for the ModuleRegister template
  /// @return descripton of the module
  static std::string getModuleDescription()
  {
    return "This module implements a Linear Advection physical model.";
  }

}; // end LinearAdvModule

//////////////////////////////////////////////////////////////////////////////

} // namespace LinearAdv
  }
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFLUID_Physics_LinearAdv_hh
