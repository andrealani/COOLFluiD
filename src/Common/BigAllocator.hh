// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Common_BIGALLOCATOR_HH
#define COOLFluiD_Common_BIGALLOCATOR_HH

#ifdef CF_HAVE_CONFIG_H
#  include "coolfluid_config.h"
#endif // CF_HAVE_CONFIG_H

#ifdef CF_HAVE_ALLOC_MMAP
#  include "Common/MemoryAllocatorMMap.hh"
#else
#  include "Common/MemoryAllocatorNormal.hh"
#endif // CF_HAVE_ALLOC_MMAP

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
    namespace Common  {

//////////////////////////////////////////////////////////////////////////////

#ifdef CF_HAVE_ALLOC_MMAP
  typedef MemoryAllocatorMMap   BigAllocator;
#else
  typedef MemoryAllocatorNormal BigAllocator;
#endif

//////////////////////////////////////////////////////////////////////////////

    } // Common
} // COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Common_BIGALLOCATOR_HH
