// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef UTILS_ARRAY2D_HH
#define UTILS_ARRAY2D_HH

#include "Common/COOLFluiD.hh"

namespace COOLFluiD {
    namespace Common {

//////////////////////////////////////////////////////////////////////////////

/// Fast, memory efficient 2D array (first dimension is
/// extendable)
/// @author Dries Kimpe
///  @TODO: iterator support, more STL like functions
///  some more template programming to help [] operators & constructors
template <typename T>
class Array2D
{
public:

    Array2D ()
: ElementSize_(0), Size_(0)
    {
    }

    Array2D (unsigned int Size1, unsigned int Size2)
: ElementSize_(Size2)
    {
Resize(Size1);
    }


    void SetElementSize (unsigned int S)
    {
cf_assert (S);
const unsigned int OldSize = GetSize();
ElementSize_ = S;
Resize (OldSize);
    }

    unsigned int GetSize () const
    {
return Size_;
    };

    // Compatibility
    unsigned int size () const
    {
return GetSize();
    }

    void Resize (unsigned int NewSize)
    {
cf_assert (ElementSize_ > 0);
Data_.resize(NewSize*ElementSize_);
Size_ = NewSize;
    }

    void Reserve (unsigned int NewSize)
    {
Data_.reserve (NewSize*ElementSize_);
    }

    T & Get (unsigned int T1, unsigned int T2)
    {
cf_assert (T1 < Size_);
cf_assert (T2 < ElementSize_);
return Data_[T1*ElementSize_+T2];
    }

    const T & Get (unsigned int T1, unsigned int T2) const
    {
cf_assert (T1 < Size_);
cf_assert (T2 < ElementSize_);
return Data_[T1*ElementSize_+T2];
    }

    T & PushBack ()
    {
Resize (Size_+1);
return Get(Size_-1, 0);
    }

    void clear ()
    {
      Resize (0);
    }

    unsigned int GetElementSize () const
    {
return ElementSize_;
    }

protected:

    /// Element size
    unsigned int ElementSize_;

    unsigned int Size_;

    /// Raw data
    std::vector<T> Data_;
};

//////////////////////////////////////////////////////////////////////////////


    }
}


#endif
