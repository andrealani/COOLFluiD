// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_hh
#define COOLFluiD_Framework_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/ExportAPI.hh"
#include "Environment/ModuleRegister.hh"

//////////////////////////////////////////////////////////////////////////////

/// Define the macro Framework_API
/// @note build system defines Framework_EXPORTS when compiling Environment files
#ifdef Framework_EXPORTS
#   define Framework_API      CF_EXPORT_API
#   define Framework_TEMPLATE
#else
#   define Framework_API      CF_IMPORT_API
#   define Framework_TEMPLATE CF_TEMPLATE_EXTERN
#endif

/// Macro defining the factories in Environment
#ifdef  CF_HAVE_CXX_EXPLICIT_TEMPLATES
#   define Framework_Factory(__fac__) namespace COOLFluiD { namespace Environment { CF_TEMPLATE_EXTERN template class Framework_API COOLFluiD::Environment::Factory<COOLFluiD::Framework::__fac__>; } }
#else
#   define Framework_Factory(__fac__)
#endif

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  /// Basic classes that form the Framework which is a basis for
  /// the implementation of external Components
  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This class defines the Module Framework
class Framework_API FrameworkLib : public Environment::ModuleRegister<FrameworkLib> {
public:

  /// Static function that returns the module name.
  /// Must be implemented for the ModuleRegister template
  /// @return name of the module
  static std::string getModuleName() {  return "Framework";  }

  /// Static function that returns the description of the module.
  /// Must be implemented for the ModuleRegister template
  /// @return descripton of the module
  static std::string getModuleDescription()
  {
    return "This module implementes the kernel functionality for the COOLFluiD framework.";
  }

}; // end FrameworkLib

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_hh
