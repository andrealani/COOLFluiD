// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Environment_EnvironmentAPI_hh
#define COOLFluiD_Environment_EnvironmentAPI_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/ExportAPI.hh"

//////////////////////////////////////////////////////////////////////////////

/// Define the macro Environment_API
/// @note build system defines Environment_EXPORTS when compiling Environment files
#ifdef Environment_EXPORTS
#   define Environment_API      CF_EXPORT_API
#   define Environment_TEMPLATE
#else
#   define Environment_API      CF_IMPORT_API
#   define Environment_TEMPLATE CF_TEMPLATE_EXTERN
#endif

/// Macro defining the factories in Environment
#ifdef  CF_HAVE_CXX_EXPLICIT_TEMPLATES
#   define Environment_Factory(__fac__) namespace COOLFluiD { namespace Environment { CF_TEMPLATE_EXTERN template class Environment_API COOLFluiD::Environment::Factory<COOLFluiD::Environment::__fac__>; } }
#else
#   define Environment_Factory(__fac__)
#endif

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Environment_EnvironmentAPI_hh
