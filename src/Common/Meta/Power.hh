// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Common_Meta_Power_hh
#define COOLFluiD_Common_Meta_Power_hh

///////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace Common {

//////////////////////////////////////////////////////////////////////////////

/// Template meta-programming function for power function
template < int N , int P >
struct Power
{
    enum { value = N * Power<N,P-1>::value };
};

template <int N>
struct Power<N,0>
{
    enum { value = 1 };
};

//////////////////////////////////////////////////////////////////////////////

  } // Utils
} // COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Common_Meta_Power_hh


