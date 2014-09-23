// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Common_CppTypeInfo_hh
#define COOLFluiD_Common_CppTypeInfo_hh

namespace COOLFluiD {
    namespace Common {

  /// void type information
  template < typename T >  struct is_void        { enum { type = 0 }; };
  template < >             struct is_void <void> { enum { type = 1 }; };

  /// integer type information
  template < typename T >  struct is_integer                  { enum { type = 0 }; };
  template <>              struct is_integer<bool>            { enum { type = 1 }; };
  template <>              struct is_integer<char>            { enum { type = 1 }; };
  template <>              struct is_integer<signed char>     { enum { type = 1 }; };
  template <>              struct is_integer<unsigned char>   { enum { type = 1 }; };
  template <>              struct is_integer<wchar_t>         { enum { type = 1 }; };
  template <>              struct is_integer<short>           { enum { type = 1 }; };
  template <>              struct is_integer<unsigned short>  { enum { type = 1 }; };
  template <>              struct is_integer<int>             { enum { type = 1 }; };
  template <>              struct is_integer<unsigned int>    { enum { type = 1 }; };
  template <>              struct is_integer<long>            { enum { type = 1 }; };
  template <>              struct is_integer<unsigned long>   { enum { type = 1 }; };

  /// floating point type information
  template < typename T >  struct is_floating                  { enum { type = 0 }; };
  template <>              struct is_floating<float>           { enum { type = 1 }; };
  template <>              struct is_floating<double>          { enum { type = 1 }; };
  template <>              struct is_floating<long double>     { enum { type = 1 }; };

  // a type is arithmetic if it is floating or integer
  template < typename T >  struct is_arithmetic                { enum { type = is_integer<T>::type || is_floating<T>::type }; };

  // a fundamental type is either void or arithmetic type
  template < typename T > struct  is_fundamental               { enum { type = is_void<T>::type || is_arithmetic<T>::type  }; };

  } // namespace Common
} // namespace COOLFluiD

#endif // CppTypeInfo
