// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef MEM_ALLOC_NORMAL_HH
#define MEM_ALLOC_NORMAL_HH

//////////////////////////////////////////////////////////////////////////////

#include "Common/MemoryAllocator.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace Common {

//////////////////////////////////////////////////////////////////////////////

/// class MemoryAllocatorNormal
///    This class uses virtual functions, because the overhead is neglectable
///    here. (And this makes it easy to switch between versions)
///    This versions uses malloc / realloc
/// @todo
///   - Integrate with autoconf:
///        * Use zero-copy version when mremap is present
///        * Fall back to malloc/realloc
///   - Exceptions
///  @author Dries Kimpe
class  Common_API  MemoryAllocatorNormal : public MemoryAllocator {
private:

    MA_Ptr _ptr;
    MA_Size _size;

public:

  MemoryAllocatorNormal (MA_Size _initsize = 0);

  virtual ~MemoryAllocatorNormal ();

  /// Return current size (bytes)
  virtual MA_Size GetSize () const
  {
    return _size;
  }

  /// Granularity
  /// We assume libC's malloc returns only aligned memory references
  virtual MA_Size GetGranularity () const
  {
    return sizeof(int);
  }

  /// Resizes memory
  /// Returns new size, will be at least NewSize:
  /// [WARNING: may change ptr!! ]
  virtual MA_Size Resize (MA_Size NewSize);

  /// Is state valid (= memory allocated?)
  virtual bool IsValid () const
  {
    return (_ptr != 0);
  }

  /// Return pointer to memory
  virtual MA_Ptr GetPtr () const
  {
    return _ptr;
  }

  /// Is resize zero copy?
  virtual bool IsZeroCopy () const
  {
    return false;
  }

}; // end class MemoryAllocatorNormal

//////////////////////////////////////////////////////////////////////////////

  } // namespace Common

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // MEM_ALLOC_NORMAL_HH
