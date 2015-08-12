// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_EmptyConvergenceMethod_EmptyConvergenceMethodAPI_hh
#define COOLFluiD_EmptyConvergenceMethod_EmptyConvergenceMethodAPI_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/ExportAPI.hh"

//////////////////////////////////////////////////////////////////////////////

/// Define the macro EmptyConvergenceMethod_API
/// @note build system defines EmptyConvergenceMethod_EXPORTS when compiling EmptyConvergenceMethod files
#ifdef EmptyConvergenceMethod_EXPORTS
#   define EmptyConvergenceMethod_API CF_EXPORT_API
#else
#   define EmptyConvergenceMethod_API CF_IMPORT_API
#endif

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_EmptyConvergenceMethod_EmptyConvergenceMethodAPI_hh
