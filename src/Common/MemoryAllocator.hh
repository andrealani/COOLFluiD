// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef MEM_ALLOC_HH
#define MEM_ALLOC_HH

//////////////////////////////////////////////////////////////////////////////

#include "Common/Exception.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace Common {

//////////////////////////////////////////////////////////////////////////////

/// Exception used in MemoryAllocator
/// @author Dries Kimpe
class  Common_API  MemoryAllocatorException : public Common::Exception {
public:

  /// Constructor
  MemoryAllocatorException (const Common::CodeLocation& where, const std::string& what = "Memory exhausted") :
    Common::Exception(where, what, "MemoryAllocatorException") {}

}; // end class MemoryAllocatorException

//////////////////////////////////////////////////////////////////////////////

///  This class uses virtual functions, because the overhead is neglectable
///  here. (And this makes it easy to switch between versions)
///  Using this interface allows a (linux specific) zero-copy implementation
/// @todo Probably this shouldn't use virtual functions, as we decide on
///        the implementation (mmap or normal) on compile time.
///  @author Dries Kimpe
class Common_API MemoryAllocator {
public:

  typedef void * MA_Ptr;
  typedef size_t MA_Size;

  /// Return current size (bytes)
  virtual MA_Size GetSize () const = 0;

  /// Granularity
  virtual MA_Size GetGranularity () const = 0;

  /// Resizes memory
  /// Returns new size, will be at least NewSize:
  /// [WARNING: may change ptr!! ]
  virtual MA_Size Resize (MA_Size NewSize) = 0;

  /// Is state valid (= memory allocated?)
  virtual bool IsValid () const = 0;

  /// Return pointer to memory
  virtual MA_Ptr GetPtr () const = 0;

  /// Is resize zero copy?
  virtual bool IsZeroCopy () const = 0;

  ///  Not so sure this is a good idea
  operator MA_Ptr () const
  {
    return GetPtr();
  };

  /// To make testing easy
  operator bool () const
  {
    return IsValid ();
  };

  /// Virtual destructor
  virtual ~MemoryAllocator ()
  {
  }

}; // end class MemoryAllocator

//////////////////////////////////////////////////////////////////////////////

  } // end Namespace Common

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // MEM_ALLOC_HH
