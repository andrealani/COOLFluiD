// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_TecplotWriter_TecplotWriterAPI_hh
#define COOLFluiD_TecplotWriter_TecplotWriterAPI_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/ExportAPI.hh"

//////////////////////////////////////////////////////////////////////////////

/// Define the macro TecplotWriter_API
/// @note build system defines TecplotWriter_EXPORTS when compiling TecplotWriter files
#ifdef TecplotWriter_EXPORTS
#   define TecplotWriter_API CF_EXPORT_API
#else
#   define TecplotWriter_API CF_IMPORT_API
#endif

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_TecplotWriter_TecplotWriterAPI_hh
