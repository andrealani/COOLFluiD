// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Common/COOLFluiD.hh"
#include "Common/MemoryAllocatorMMap.hh"

#include <unistd.h>
#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace Common {

//////////////////////////////////////////////////////////////////////////////

MemoryAllocatorMMap::~MemoryAllocatorMMap ()
{
  Free ();
}

//////////////////////////////////////////////////////////////////////////////

MemoryAllocatorMMap::MemoryAllocatorMMap (MA_Size InitialSize)
    : DataPtr(0),CurrentSize(0),FileDesc(-1)
{
    Alloc(InitialSize);
}

//////////////////////////////////////////////////////////////////////////////

MemoryAllocatorMMap::MA_Size MemoryAllocatorMMap::GetGranularity () const
{
  // alternate implementation:
  cf_assert ( int(sysconf(_SC_PAGESIZE)) == getpagesize () );
  return (int) sysconf(_SC_PAGESIZE);
}

//////////////////////////////////////////////////////////////////////////////

void MemoryAllocatorMMap::Alloc (MA_Size size)
{
  cf_assert (DataPtr==0);
  cf_assert (FileDesc == -1);

  // Minimum size
  if (size==0)
      size=1;

  /* map /dev/zero */
  FileDesc = open ("/dev/zero", O_RDWR);
  if (FileDesc < 0)
    throw MemoryAllocatorException (FromHere());

  /* do memory mapping */
  DataPtr = mmap (0, size, PROT_READ|PROT_WRITE, MAP_PRIVATE, FileDesc, 0);
  if (DataPtr == MAP_FAILED)
  {
    DataPtr=0;
    throw MemoryAllocatorException (FromHere());
  }
  cf_assert (DataPtr!=0);

  CurrentSize = size;
}

//////////////////////////////////////////////////////////////////////////////

void MemoryAllocatorMMap::Free ()
{
  cf_assert (DataPtr != 0);
  cf_assert (FileDesc != -1);

  int Ret = munmap (DataPtr, CurrentSize);
  if (Ret < 0)
    throw MemoryAllocatorException (FromHere());
  DataPtr = 0;

  Ret = close (FileDesc);
  cf_assert (Ret>= 0);
  FileDesc = -1;
  CurrentSize = 0;
}

//////////////////////////////////////////////////////////////////////////////

MemoryAllocatorMMap::MA_Size MemoryAllocatorMMap::Resize (MA_Size NewSize)
{
  if (NewSize==0)
      NewSize=1;

  if (CurrentSize==NewSize)
      return CurrentSize;

  // Resize current allocation
  cf_assert (DataPtr != 0);
  cf_assert (FileDesc != -1);

  MA_Ptr NewData = mremap (DataPtr, CurrentSize, NewSize, MREMAP_MAYMOVE);
  if (NewData == MA_Ptr(-1))
    throw MemoryAllocatorException (FromHere());

  CurrentSize = NewSize;
  DataPtr = NewData;
  return CurrentSize;
}

//////////////////////////////////////////////////////////////////////////////

bool MemoryAllocatorMMap::IsValid () const
{
  return (FileDesc != -1);
}

//////////////////////////////////////////////////////////////////////////////

} // End namespace Common

} // End namespace COOLFluiD
