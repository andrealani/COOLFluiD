// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_ForwardEuler_ForwardEulerAPI_hh
#define COOLFluiD_ForwardEuler_ForwardEulerAPI_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/ExportAPI.hh"

//////////////////////////////////////////////////////////////////////////////

/// Define the macro ForwardEuler_API
/// @note build system defines ForwardEuler_EXPORTS when compiling ForwardEuler files
#ifdef ForwardEuler_EXPORTS
#   define ForwardEuler_API CF_EXPORT_API
#else
#   define ForwardEuler_API CF_IMPORT_API
#endif

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_ForwardEuler_ForwardEulerAPI_hh
