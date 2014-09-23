// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_MarcoTest_MarcoTestAPI_hh
#define COOLFluiD_MarcoTest_MarcoTestAPI_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/ExportAPI.hh"

//////////////////////////////////////////////////////////////////////////////

/// Define the macro MarcoTest_API
/// @note build system defines MarcoTest_EXPORTS when compiling MarcoTest files
#ifdef MarcoTest_EXPORTS
#   define MarcoTest_API CF_EXPORT_API
#else
#   define MarcoTest_API CF_IMPORT_API
#endif

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_MarcoTest_MarcoTestAPI_hh
