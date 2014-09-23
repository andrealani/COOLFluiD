// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_ShapeFunctions_ShapeFunctionsAPI_hh
#define COOLFluiD_ShapeFunctions_ShapeFunctionsAPI_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/ExportAPI.hh"

//////////////////////////////////////////////////////////////////////////////

/// Define the macro ShapeFunctions_API
/// @note build system defines ShapeFunctions_EXPORTS when compiling ShapeFunctions files
#ifdef ShapeFunctions_EXPORTS
#   define ShapeFunctions_API CF_EXPORT_API
#else
#   define ShapeFunctions_API CF_IMPORT_API
#endif

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_ShapeFunctions_ShapeFunctionsAPI_hh
