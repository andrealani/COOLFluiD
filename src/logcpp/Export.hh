// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef CF_logcpp_API_HH
#define CF_logcpp_API_HH

//////////////////////////////////////////////////////////////////////////////

#include "Common/ExportAPI.hh"

//////////////////////////////////////////////////////////////////////////////

/// Define the macro logcpp_API
/// @note build system defines logcpp_EXPORTS when compiling logcpp files
#ifdef logcpp_EXPORTS
#   define logcpp_API CF_EXPORT_API
#else
#   define logcpp_API CF_IMPORT_API
#endif

//////////////////////////////////////////////////////////////////////////////

#endif // CF_logcpp_API_HH

