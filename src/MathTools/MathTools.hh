// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFLUID_MATHTOOLSAPI_HH
#define COOLFLUID_MATHTOOLSAPI_HH

//////////////////////////////////////////////////////////////////////////////

#include "Common/ExportAPI.hh"

//////////////////////////////////////////////////////////////////////////////

/// Define the macro MathTools_API
/// @note build system defines MathTools_EXPORTS when compiling MathTools files
#ifdef MathTools_EXPORTS
#   define MathTools_API      CF_EXPORT_API
#   define MathTools_TEMPLATE
#else
#   define MathTools_API      CF_IMPORT_API
#   define MathTools_TEMPLATE CF_TEMPLATE_EXTERN
#endif

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  /// Basic Classes for Mathematical applications used by the COOLFluiD
  namespace MathTools {}

} // namespace MathTools

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFLUID_MATHTOOLSAPI_HH
