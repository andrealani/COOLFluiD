// Copyright (C) 2016 KU Leuven, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFLUID_CGNSWriter_CGNSWriterAPI_hh
#define COOLFLUID_CGNSWriter_CGNSWriterAPI_hh

#include "Common/ExportAPI.hh"

// Define the macro CGNSWriter_API
// @note build system defines CGNSWriter_EXPORTS when compiling CGNSWriter files
#ifdef CGNSWriter_EXPORTS
#   define CGNSWriter_API CF_EXPORT_API
#else
#   define CGNSWriter_API CF_IMPORT_API
#endif

#endif // COOLFLUID_CGNSWriter_CGNSWriterAPI_hh
