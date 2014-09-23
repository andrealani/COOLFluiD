// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Common_CommonAPI_hh
#define COOLFluiD_Common_CommonAPI_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/ExportAPI.hh"

//////////////////////////////////////////////////////////////////////////////

/// Define the macro Common_API
/// @note build system defines Common_EXPORTS when compiling Common files
#ifdef Common_EXPORTS
#   define Common_API CF_EXPORT_API
#else
#   define Common_API CF_IMPORT_API
#endif

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Common_CommonAPI_hh
