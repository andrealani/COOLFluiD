// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Common_PtrAlloc_hh
#define COOLFluiD_Common_PtrAlloc_hh

/// @note This header should be included by including COOLFluiD.hh instead.

//////////////////////////////////////////////////////////////////////////////

#include "Common/CFAssert.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

//////////////////////////////////////////////////////////////////////////////

///  Definition of CFNULL
#define CFNULL 0

//////////////////////////////////////////////////////////////////////////////

  /// @brief Deletes a pointer and makes sure it is set to CFNULL afterwards
  /// It would not have to check for CFNULL before deletion, as
  /// deleting a null is explicitely allowed by the standard.
  /// Nevertheless it does check, to avoid problems with not so compliant compilers.
  /// Do not use this function with data allocate with new [].
  /// @author Tiago Quintino
  /// @pre ptr has been allocated with operator new
  /// @param ptr pointer to be deleted
  /// @post  ptr equals CFNULL
  template <class TYPE>
  void deletePtr(TYPE*& ptr)
  {
    if (ptr != CFNULL)
    {
      delete ptr;
      ptr = CFNULL;
    }
    cf_assert(ptr == CFNULL);
  }

  /// @brief Deletes a pointer and makes sure it is set to CFNULL afterwards
  /// It would not have to check for CFNULL before deletion, as
  /// deleting a null is explicitely allowed by the standard.
  /// Nevertheless it does check, to avoid problems with not so compliant compilers.
  /// Do not use this function with data allocate with new.
  /// @author Tiago Quintino
  /// @pre ptr has been allocated with operator new []
  /// @param ptr pointer to be deleted
  /// @post  ptr equals CFNULL
  template <class TYPE>
  void deletePtrArray(TYPE*& ptr)
  {
    if (ptr != CFNULL)
    {
      delete [] ptr;
      ptr = CFNULL;
    }
    cf_assert(ptr == CFNULL);
  }

//////////////////////////////////////////////////////////////////////////////

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Common_PtrAlloc_hh
