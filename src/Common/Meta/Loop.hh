// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Common_Meta_Loop_hh
#define COOLFluiD_Common_Meta_Loop_hh

///////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace Common {

//////////////////////////////////////////////////////////////////////////////

/// Template meta-programming version of a Loop structure
/// FUNC trait class should implement a static function exec
/// taking no parameters
template< typename FUNC, unsigned int I >
class Loop
{ 
  public:
    static inline void run ( )
    {
      Loop < FUNC, I-1 >::run();
      FUNC::template exec<I-1>();
    }
};

/// Specialization provides end case of recursion
template< typename FUNC >
class Loop < FUNC, 0 >
{
  public:
    static inline void run () { }
};

//////////////////////////////////////////////////////////////////////////////

/// Template meta-programming version of a Loop structure
/// FUNC trait class should implement a static function exec
/// taking 1 parameter parameters
template< typename FUNC, unsigned int I >
class Loop1
{ 
  public:
    template < typename P1 >
    static inline void run ( P1& p1 )
    {
      Loop1 < FUNC, I-1 >::run(p1);
      FUNC::template exec<I-1>(p1);
    }
};

/// Specialization provides end case of recursion
template< typename FUNC >
class Loop1 < FUNC, 0 >
{
  public:
    template < typename P1 >
    static inline void run ( P1& p1 ) { }
};

//////////////////////////////////////////////////////////////////////////////

/// Template meta-programming version of a Loop structure
/// FUNC trait class should implement a static function exec
/// taking 2 parameter parameters
template< typename FUNC, unsigned int I >
class Loop2
{ 
  public:
    template < typename P1, typename P2 >
    static inline void run ( P1& p1, P2& p2 )
    {
      Loop2 < FUNC, I-1 >::run(p1,p2);
      FUNC::template exec<I-1>(p1,p2);
    }
};

/// Specialization provides end case of recursion
template< typename FUNC >
class Loop2 < FUNC, 0 >
{
  public:
    template < typename P1, typename P2 >
    static inline void run ( P1& p1, P2& p2 ) { }
};

//////////////////////////////////////////////////////////////////////////////

/// Template meta-programming version of a Loop structure
/// FUNC trait class should implement a static function exec
/// taking 3 parameters
template< typename FUNC, unsigned int I >
class Loop3
{ 
  public:
    template < typename P1, typename P2 , typename P3 >
    static inline void run ( P1& p1, P2& p2, P3& p3 )
    {
      Loop3 < FUNC, I-1 >::run(p1,p2,p3);
      FUNC::template exec<I-1>(p1,p2,p3);
    }
};

/// Specialization provides end case of recursion
template< typename FUNC >
class Loop3 < FUNC, 0 >
{
  public:
    template < typename P1, typename P2, typename P3 >
    static inline void run ( P1& p1, P2& p2, P3& p3 ) { }
};

//////////////////////////////////////////////////////////////////////////////

/// Template meta-programming version of a Loop structure
/// FUNC trait class should implement a static function exec
/// taking 4 parameters
template< typename FUNC, unsigned int I >
class Loop4
{
  public:
    template < typename P1, typename P2, typename P3, typename P4 >
    static inline void run ( P1& p1, P2& p2, P3& p3, P4& p4 )
    {
      Loop4 < FUNC, I-1 >::run(p1,p2,p3,p4);
      FUNC::template exec<I-1>(p1,p2,p3,p4);
    }
};

/// Specialization provides end case of recursion
template< typename FUNC >
class Loop4 < FUNC, 0 >
{
  public:
    template < typename P1, typename P2, typename P3, typename P4 >
    static inline void run ( P1& p1, P2& p2, P3& p3, P4& p4 ) { }
};

//////////////////////////////////////////////////////////////////////////////

  } // Utils
} // COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Common_Meta_Loop_hh


