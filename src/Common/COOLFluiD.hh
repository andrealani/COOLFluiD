// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_COOLFluiD_hh
#define COOLFluiD_COOLFluiD_hh

//////////////////////////////////////////////////////////////////////////////

#ifdef CF_HAVE_CONFIG_H
#  include "coolfluid_config.h"
#endif // CF_HAVE_CONFIG_H

//////////////////////////////////////////////////////////////////////////////

#include "Common/Compatibility.hh"
#include "Common/StlHeaders.hh"
#include "Common/Common.hh"
#include "Common/CFAssert.hh"
#include "Common/DemangledTypeID.hh"
#include "Common/PtrAlloc.hh"

//////////////////////////////////////////////////////////////////////////////

/// Definition of COOLFluiD namespace.
/// @author Tiago Quintino
/// @author Andrea Lani
namespace COOLFluiD {

//////////////////////////////////////////////////////////////////////////////

  /// Definition of the basic types for possible portability conflicts

  /// typedef for float
  typedef float              CFfloat;
  /// typedef for double
  typedef double             CFdouble;
  /// typedef for long double
  typedef long double        CFldouble;
  
#ifdef CF_HAVE_LONG
  /// typedef for int
  typedef long long int      CFint;
  /// typedef for unsigned int
  typedef long long int      CFuint;
  //typedef long unsigned int  CFuint;
#else
  /// typedef for int
  typedef int                CFint;
  /// typedef for unsigned int
  typedef unsigned int       CFuint;
#endif
  
  /// typedef for char
  typedef char               CFchar;

  /// Enumeration of the dimensions
  enum CFDim         {DIM_0D, DIM_1D, DIM_2D, DIM_3D};

  /// Enumeration of the coordinates indexes
  enum CoordXYZ       {XX, YY, ZZ};
  
  /// Enumeration of the reference coordinates indexes
  enum CoordRefXiEtaZeta    {KSI, ETA, ZTA};

  /// Enumeration of the device types
  enum DeviceType {CPU=0, GPU=1};

  /// class to be used to define a default type
  class NOTYPE {};
  
//////////////////////////////////////////////////////////////////////////////

/// Definition of the default precision
#ifdef CF_PRECISION_LONG_DOUBLE
  typedef CFldouble CFreal;
#else
  #ifdef CF_PRECISION_DOUBLE
    typedef CFdouble CFreal;
  #else
    #ifdef CF_PRECISION_SINGLE
      typedef CFfloat CFreal;
    #endif
  #endif
#endif
// if nothing defined, use double
#if !defined CF_PRECISION_DOUBLE && !defined CF_PRECISION_SINGLE && !defined CF_PRECISION_LONG_DOUBLE
  typedef CFdouble CFreal;
#endif

typedef std::complex<CFreal>  CFcomplex;

//////////////////////////////////////////////////////////////////////////////

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_COOLFLUID_hh
