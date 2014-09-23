// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Common_ExportAPI_hh
#define COOLFluiD_Common_ExportAPI_hh

//////////////////////////////////////////////////////////////////////////////

#ifdef CF_HAVE_CONFIG_H
#  include "coolfluid_config.h"
#endif // CF_HAVE_CONFIG_H

//////////////////////////////////////////////////////////////////////////////

/// Define the macros:
///    CF_LOCAL_API
///    CF_IMPORT_API
///    CF_EXPORT_API

// Define or not the extern keyword for templates
#ifdef CF_HAVE_EXTERN_TEMPLATES
#   define CF_TEMPLATE_EXTERN extern
#else
#   define CF_TEMPLATE_EXTERN
#endif


// Visual Studio :
//    all symbols are local (invisible) by default
#  if defined (_MSC_VER)
#    define CF_LOCAL_API
#    define CF_IMPORT_API    __declspec(dllimport)
#    define CF_EXPORT_API    __declspec(dllexport)
#  endif

// GNU Compiler :  all symbols are global (visible) by default
#  if defined (__GNUC__) // && defined (__unix__)
#    define CF_LOCAL_API    __attribute__ ((visibility("hidden")))
#    define CF_EXPORT_API   __attribute__ ((visibility("default")))
#    define CF_IMPORT_API   __attribute__ ((visibility("default")))
#  endif

// For unrecognized compilers
#  ifndef CF_LOCAL_API
#    if defined ( CF_ACCEPT_ANY_COMPILER )
#      define CF_LOCAL_API
#      define CF_EXPORT_API
#      define CF_IMPORT_API
#    else
#      error "Unrecognised compiler and / or platform"
#    endif
#  endif

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Common_ExportAPI_hh
