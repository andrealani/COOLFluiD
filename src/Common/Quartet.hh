// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Common_Quartet_hh
#define COOLFluiD_Common_Quartet_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/COOLFluiD.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Common {

//////////////////////////////////////////////////////////////////////////////

template <typename T1,
          typename T2,
          typename T3,
          typename T4> class Quartet;

template <typename T1,
          typename T2,
          typename T3,
          typename T4>
class Quartet {
public:

  typedef T1 first_type;
  typedef T2 second_type;
  typedef T3 third_type;
  typedef T4 fourth_type;

public:

  T1 first;
  T2 second;
  T3 third;
  T4 fourth;

public:

  ///  Empty Constructor
  Quartet() : first(T1()), second(T2()), third(T3()), fourth(T4()) {}

  ///  Constructor
  Quartet(const T1& a, const T2& b, const T3& c, const T4& d) : first(a), second(b), third(c), fourth(d) {}

  ///  Copy Constructor with conversion of types
  template <typename U1, typename U2, typename U3, typename U4>
  Quartet(const Quartet<U1, U2, U3, U4>& p) : first(p.first), second(p.second), third(p.third), fourth(p.fourth) {}

}; // end class Quartet

template <class T1, class T2, class T3, class T4>
inline bool operator> (const Quartet<T1, T2, T3, T4>& p1, const Quartet<T1, T2, T3, T4>& p2)
{
  bool result;

  if (p1.first == p2.first)
    if (p1.second == p2.second)
      if (p1.third == p2.third)
        result = (p1.fourth > p2.fourth);
      else
        result = (p1.third > p2.third);
    else
      result = (p1.second > p2.second);
  else
    result = (p1.first > p2.first);

  return result;
}

template <class T1, class T2, class T3, class T4>
inline bool operator< (const Quartet<T1, T2, T3, T4>& p1, const Quartet<T1, T2, T3, T4>& p2)
{
  bool result;

  if (p1.first == p2.first)
    if (p1.second == p2.second)
      if (p1.third == p2.third)
        result = (p1.fourth < p2.fourth);
      else
        result = (p1.third < p2.third);
    else
      result = (p1.second < p2.second);
  else
    result = (p1.first < p2.first);

  return result;
}
template <class T1, class T2, class T3, class T4>
inline bool operator==(const Quartet<T1, T2, T3, T4>& x, const Quartet<T1, T2, T3, T4>& y)
{
  return x.first == y.first && x.second == y.second  && x.third == y.third && x.fourth == y.fourth;
}

template <class T1, class T2, class T3, class T4>
inline bool operator!=(const Quartet<T1, T2, T3, T4>& x, const Quartet<T1, T2, T3, T4>& y) {
  return !(x == y);
}

template <class T1, class T2, class T3, class T4>
inline Quartet<T1, T2, T3, T4> make_Quartet(const T1& x1, const T2& x2, const T3& x3, const T4& x4)
{
  return Quartet<T1, T2, T3, T4>(x1, x2, x3, x4);
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Common

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // end COOLFluiD_Common_Quartet_hh
