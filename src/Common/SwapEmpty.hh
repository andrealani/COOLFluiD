// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef UTILS_SWAPEMPTY_HH
#define UTILS_SWAPEMPTY_HH

namespace COOLFluiD {
    namespace Common  {
//////////////////////////////////////////////////////////////////////////////

/// Really free storage associated with a vector
/// Can also be used for other things
template <typename T>
inline void SwapEmpty (T & E)
{
    T().swap (E);
//    std::swap(E, T());
}

//////////////////////////////////////////////////////////////////////////////
    }
}

#endif
