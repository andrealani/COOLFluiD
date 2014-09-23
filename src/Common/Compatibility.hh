// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Common_Compatibility_hh
#define COOLFluiD_Common_Compatibility_hh

/// @note This header should be included by including COOLFluiD.hh instead.

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

//////////////////////////////////////////////////////////////////////////////

/// Macro necessary to cope with operator<< in combination with templates.
/// Compilers that have to have this set: cxx ( SGI )
/// Compilers that should not have this set: gcc, CC, icc (8.0)
#ifdef CXX_NEEDS_FRIEND_TMPL_DECL
#  define LTGT <>
#else
#  define LTGT
#endif

/// Macro necessary for compilers, that do not define the __FUNCTION__ variable
#ifndef CF_HAVE_FUNCTION_DEF
#  define __FUNCTION__ ""
#endif

/// Macro necessary for compiling functions to be used on device and/or host
#ifndef CF_HAVE_CUDA
#define HOST_DEVICE
#else
#define HOST_DEVICE __host__ __device__
#endif 
  
//////////////////////////////////////////////////////////////////////////////

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Common_Compatibility_hh
