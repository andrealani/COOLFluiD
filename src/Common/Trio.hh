// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Common_Trio_hh
#define COOLFluiD_Common_Trio_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/COOLFluiD.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Common {

//////////////////////////////////////////////////////////////////////////////

template <typename T1,
          typename T2,
          typename T3> class Trio;

template <typename T1,
          typename T2,
          typename T3> std::istream& operator>> (std::istream& in, Trio<T1, T2, T3>&);

template <typename T1,
          typename T2,
          typename T3> std::ostream& operator<< (std::ostream& in, Trio<T1, T2, T3>&);

template <typename T1,
          typename T2,
          typename T3>
class Trio {
public:

  typedef T1 first_type;
  typedef T2 second_type;
  typedef T3 third_type;

public:

  T1 first;
  T2 second;
  T3 third;

public:

  ///  Empty Constructor
  Trio() : first(T1()), second(T2()), third(T3())
  {
  }

  ///  Constructor
  Trio(const T1& a, const T2& b, const T3& c) : first(a), second(b), third(c)
  {
  }

  ///  Copy Constructor with conversion of types
  template <typename U1, typename U2, typename U3>
  Trio(const Trio<U1, U2, U3>& p) : first(p.first), second(p.second), third(p.third)
  {
  }

  friend std::istream& operator>> (std::istream& in, Trio<T1, T2, T3>& q)
  {
    in >> q.first  >> q.second >> q.third;
    return in;
  }

  friend std::ostream& operator<< (std::ostream& out, Trio<T1, T2, T3>& q)
  {
    out << q.first << " " << q.second << " " << q.third;
    return out;
  }

}; // end class Trio

template <class T1, class T2, class T3>
inline bool operator==(const Trio<T1, T2, T3>& x, const Trio<T1, T2, T3>& y)
{
  return x.first == y.first && x.second == y.second && x.third == y.third;
}

template <class T1, class T2, class T3>
inline bool operator!=(const Trio<T1, T2, T3>& x, const Trio<T1, T2, T3>& y) {
  return !(x == y);
}

template <class T1, class T2, class T3>
inline Trio<T1, T2, T3> make_Trio(const T1& x, const T2& y, const T3& z)
{
  return Trio<T1, T2, T3>(x, y, z);
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Common

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // end COOLFluiD_Common_Trio_hh
