// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_CFmeshFileReader_CFmeshFileReaderAPI_hh
#define COOLFluiD_CFmeshFileReader_CFmeshFileReaderAPI_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/ExportAPI.hh"

//////////////////////////////////////////////////////////////////////////////

/// Define the macro CFmeshFileReader_API
/// @note build system defines CFmeshFileReader_EXPORTS when compiling CFmeshFileReader files
#ifdef CFmeshFileReader_EXPORTS
#   define CFmeshFileReader_API CF_EXPORT_API
#else
#   define CFmeshFileReader_API CF_IMPORT_API
#endif

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_CFmeshFileReader_CFmeshFileReaderAPI_hh
