// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Common_MemFunArg_hh
#define COOLFluiD_Common_MemFunArg_hh

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Common {

//////////////////////////////////////////////////////////////////////////////

/// Class to allow STL algorithms to call member functions on this pointer
/// This is similar to mem_fun_t from the STL standard
template < typename TYPE, typename RET , typename ARG >
class mem_fun_arg_t
{
public: // typedefs

  /// @c caller_type is the type of the caller object
  typedef TYPE caller_type;
  /// @c argument_type is the type of the argument
  typedef ARG argument_type;
  /// @c result_type is the return type
  typedef RET result_type;

public: // functions

  /// constructor
  explicit mem_fun_arg_t( TYPE& c, RET ( TYPE::*pf )( ARG ) ) : m_c(c), m_f(pf) {}

  /// this is the operator we need to add in order for STL algorithms to function
  RET operator()(ARG a) const  { return (m_c.*m_f)(a); }

private: // data

  /// reference to the caller object
  TYPE& m_c;
  /// storage of the pointer to the member function to be called
  RET (TYPE::*m_f)(ARG);

}; // end class mem_fun_arg_t

//////////////////////////////////////////////////////////////////////////////

/// Class to allow STL algorithms to call member functions on this pointer
/// This is similar to mem_fun_t from the STL standard
template < typename TYPE, typename RET , typename ARG >
class const_mem_fun_arg_t
{
public: // typedefs

  /// @c caller_type is the type of the caller object
  typedef TYPE caller_type;
  /// @c argument_type is the type of the argument
  typedef ARG argument_type;
  /// @c result_type is the return type
  typedef RET result_type;

public: // functions

  /// constructor
  explicit const_mem_fun_arg_t( const TYPE& c, RET ( TYPE::*pf )( ARG ) const ) : m_c(c), m_f(pf) {}

  /// this is the operator we need to add in order for STL algorithms to function
  RET operator()(ARG a) const  { return (m_c.*m_f)(a); }

private: // data

  /// reference to the caller object
  const TYPE& m_c;
  /// storage of the pointer to the member function to be called
  RET (TYPE::*m_f)(ARG) const;

}; // end class mem_fun_arg_t

//////////////////////////////////////////////////////////////////////////////

/// Helper function to allow easy use of the mem_func_t class
/// This is similar to mem_fun from the STL standard
template < typename TYPE, typename RET , typename ARG >
  inline mem_fun_arg_t<TYPE,RET,ARG>
  mem_fun_arg( TYPE& c, RET (TYPE::*f)(ARG) )
  {
    return mem_fun_arg_t<TYPE,RET,ARG>(c,f);
  }

//////////////////////////////////////////////////////////////////////////////

/// Helper function to allow easy use of the mem_func_t class
/// This is similar to mem_fun from the STL standard
template < typename TYPE, typename RET , typename ARG >
  inline const_mem_fun_arg_t<TYPE,RET,ARG>
  mem_fun_arg( const TYPE& c, RET (TYPE::*f)(ARG) const )
  {
    return const_mem_fun_arg_t<TYPE,RET,ARG>(c,f);
  }

//////////////////////////////////////////////////////////////////////////////

  } // namespace Common

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

// #include "Common/MemFunArg.ci"

//////////////////////////////////////////////////////////////////////////////

#endif //COOLFluiD_Common_MemFunArg_hh
