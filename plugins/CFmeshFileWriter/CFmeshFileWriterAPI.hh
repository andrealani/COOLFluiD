// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_CFmeshFileWriter_CFmeshFileWriterAPI_hh
#define COOLFluiD_CFmeshFileWriter_CFmeshFileWriterAPI_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/ExportAPI.hh"

//////////////////////////////////////////////////////////////////////////////

/// Define the macro CFmeshFileWriter_API
/// @note build system defines CFmeshFileWriter_EXPORTS when compiling CFmeshFileWriter files
#ifdef CFmeshFileWriter_EXPORTS
#   define CFmeshFileWriter_API CF_EXPORT_API
#else
#   define CFmeshFileWriter_API CF_IMPORT_API
#endif

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_CFmeshFileWriter_CFmeshFileWriterAPI_hh
