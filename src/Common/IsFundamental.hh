// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Common_IsFundamental_hh
#define COOLFluiD_Common_IsFundamental_hh

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Common {

//////////////////////////////////////////////////////////////////////////////

/// This series of classes extrapolated from the GNU C++ system files (version 4.3)
/// provide utility function to do some type checking based on metaprogramming
/// techniques
/// @author Gabriel Dos Reis <dosreis@cmla.ens-cachan.fr>
    struct TrueType { };
    struct FalseType { };

    template<bool>
    struct TruthType
    { typedef FalseType TTYPE; };

    template<>
    struct TruthType<true>
    { typedef TrueType TTYPE; };

    // N.B. The conversions to bool are needed due to the issue
    // explained in c++/19404.
    template<class SP, class TP>
    struct TraitorClass
    {
      enum { VALUE = bool(SP::VALUE) || bool(TP::VALUE) };
      typedef typename TruthType<VALUE>::TTYPE TTYPE;
    };

    // Compare for equality of types.
    template<typename, typename>
    struct AreEqual
    {
      enum { VALUE = 0 };
      typedef FalseType TTYPE;
    };

    template<typename TP>
    struct AreEqual<TP, TP>
    {
      enum { VALUE = 1 };
      typedef TrueType TTYPE;
    };

    // Holds if the template-argument is a void type.
    template<typename TP>
    struct IsVoid
    {
      enum { VALUE = 0 };
      typedef FalseType TTYPE;
    };

    template<>
    struct IsVoid<void>
    {
      enum { VALUE = 1 };
      typedef TrueType TTYPE;
    };

    //
    // Integer types
    //
    template<typename TP>
    struct IsInteger
    {
      enum { VALUE = 0 };
      typedef FalseType TTYPE;
    };

    // Thirteen specializations (yes there are eleven standard integer
    // types; 'long long' and 'unsigned long long' are supported as
    // extensions)
    // template<>
//     struct IsInteger<bool>
//     {
//       enum { VALUE = 1 };
//       typedef TrueType TTYPE;
//     };

//     template<>
//     struct IsInteger<char>
//     {
//       enum { VALUE = 1 };
//       typedef TrueType TTYPE;
//     };

//     template<>
//     struct IsInteger<signed char>
//     {
//       enum { VALUE = 1 };
//       typedef TrueType TTYPE;
//     };

//     template<>
//     struct IsInteger<unsigned char>
//     {
//       enum { VALUE = 1 };
//       typedef TrueType TTYPE;
//     };

    template<>
    struct IsInteger<short>
    {
      enum { VALUE = 1 };
      typedef TrueType TTYPE;
    };

    template<>
    struct IsInteger<unsigned short>
    {
      enum { VALUE = 1 };
      typedef TrueType TTYPE;
    };

    template<>
    struct IsInteger<int>
    {
      enum { VALUE = 1 };
      typedef TrueType TTYPE;
    };

    template<>
    struct IsInteger<unsigned int>
    {
      enum { VALUE = 1 };
      typedef TrueType TTYPE;
    };

    template<>
    struct IsInteger<long>
    {
      enum { VALUE = 1 };
      typedef TrueType TTYPE;
    };

    template<>
    struct IsInteger<unsigned long>
    {
      enum { VALUE = 1 };
      typedef TrueType TTYPE;
    };

    template<>
    struct IsInteger<long long>
    {
      enum { VALUE = 1 };
      typedef TrueType TTYPE;
    };

    template<>
    struct IsInteger<unsigned long long>
    {
      enum { VALUE = 1 };
      typedef TrueType TTYPE;
    };

    //
    // Floating point types
    //
    template<typename TP>
    struct IsFloat
    {
      enum { VALUE = 0 };
      typedef FalseType TTYPE;
    };

    // three specializations (float, double and 'long double')
    template<>
    struct IsFloat<float>
    {
      enum { VALUE = 1 };
      typedef TrueType TTYPE;
    };

    template<>
    struct IsFloat<double>
    {
      enum { VALUE = 1 };
      typedef TrueType TTYPE;
    };

    template<>
    struct IsFloat<long double>
    {
      enum { VALUE = 1 };
      typedef TrueType TTYPE;
    };

    //
    // Pointer types
    //
    template<typename TP>
    struct IsPointer
    {
      enum { VALUE = 0 };
      typedef FalseType TTYPE;
    };

    template<typename TP>
    struct IsPointer<TP*>
    {
      enum { VALUE = 1 };
      typedef TrueType TTYPE;
    };

    //
    // An arithmetic type is an integer type or a floating point type
    //
    template<typename TP>
    struct IsArithmetic
      : public TraitorClass<IsInteger<TP>, IsFloat<TP> >
    { };

    //
    // A fundamental type is `void' or and arithmetic type
    //
    template<typename TP>
    struct IsFundamental
      : public TraitorClass<IsVoid<TP>, IsArithmetic<TP> >
    { };

    //
    // A scalar type is an arithmetic type or a pointer type
    //
    template<typename TP>
    struct IsScalar
    : public TraitorClass<IsArithmetic<TP>, IsPointer<TP> >
    { };

    //
    // For use in std::copy and std::find overloads for streambuf iterators.
    //
    template<typename TP>
    struct IsChar
    {
      enum { VALUE = 0 };
      typedef FalseType TTYPE;
    };

    template<>
    struct  IsChar<char>
    {
      enum { VALUE = 1 };
      typedef TrueType TTYPE;
    };

    template<typename TP>
    struct IsByte
    {
      enum { VALUE = 0 };
      typedef FalseType TTYPE;
    };

    template<>
    struct IsByte<char>
    {
      enum { VALUE = 1 };
      typedef TrueType TTYPE;
    };

    template<>
    struct IsByte<signed char>
    {
      enum { VALUE = 1 };
      typedef TrueType TTYPE;
    };

    template<>
    struct IsByte<unsigned char>
    {
      enum { VALUE = 1 };
      typedef TrueType TTYPE;
    };

//////////////////////////////////////////////////////////////////////////////

  }
}

//////////////////////////////////////////////////////////////////////////////

#endif //  COOLFluiD_Common_IsFundamental_hh
