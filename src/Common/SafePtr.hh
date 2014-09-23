// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Common_SafePtr_hh
#define COOLFluiD_Common_SafePtr_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/COOLFluiD.hh"
#include "Common/FailedCastException.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
namespace Common {

  template < class TYPE > class SharedPtr;

//////////////////////////////////////////////////////////////////////////////

/// Class that defines a pointer that allows pass onto
/// plugin's pointers that cannot be used for deletion.
/// A safe pointer can be a null pointer, so it can bei unsafe
/// to dereference from. If it would be fully safe,
/// it would be just like a reference hence not needed.
/// @author Tiago Quintino
template < typename TYPE >
class SafePtr {

  friend class SharedPtr<TYPE>;

public: // functions

  /// Constructor.
  SafePtr();

  /// Constructor.
  /// Allows implicit conversions.
  /// @param ptr naked pointer
  SafePtr(TYPE* ptr);

  /// Copy constructor.
  /// @param other the pointer from which to copy
  SafePtr(const SafePtr& other);

  /// Destructor.
  ~SafePtr();

  /// Reset the ptr
  void reset(TYPE* other);

  /// Reset the ptr
  /// @post m_ptr = other
  void reset(const SafePtr& other);

  /// Dynamic casts the m_ptr to the TYPE
  /// @throw Common::FailedCastException
  template < typename DTYPE > SafePtr<DTYPE> d_castTo() const;

  /// Check if it is CFNULL
  /// @return true if m_ptr == CFNULL
  inline bool isNull() const;

  /// Check if it is Not CFNULL
  /// @return true if m_ptr != CFNULL
  inline bool isNotNull() const;

  /// Overloading of "=" with SafePtr
  /// @param other pointer from which to copy
  /// @return a reference to this pointer
  inline const SafePtr& operator= (const SafePtr& other);

  /// Overloading of "=" with naked pointers
  /// @param other pointer from which to copy
  /// @return a reference to this pointer
  inline const SafePtr& operator= (TYPE* other);

  /// Overloading of "=="
  /// @param other pointer to which to compare
  /// @return true if the wrapped pointer is the same
  inline bool operator== (const SafePtr& other) const;

  /// Overloading of "!="
  /// @param other pointer to which to compare
  /// @return true if the wrapped pointer is the same
  inline bool operator!= (const SafePtr& other) const;

  /// Overloading of "<"
  /// @param other pointer to which to compare
  /// @return true if the wrapped pointer is the same
  inline bool operator< (const SafePtr& other) const;

  /// Overloading of ">"
  /// @param other pointer to which to compare
  /// @return true if the wrapped pointer is the same
  inline bool operator> (const SafePtr& other) const;

  /// Overloading of "*"
  /// @return reference to the object pointed
  TYPE& operator*() const {return *m_ptr;}

  /// Overloading of "->"
  /// Notice that this is safe because the following code would not compile:
  /// \code
  ///    SafePtr<Obj> p = getObjPtr();
  ///    delete p->();
  /// \endcode
  TYPE* operator->() const {return m_ptr;}

private: // data

  /// Raw pointer to the object.
  TYPE* m_ptr;

public: // nested classes

  /// Class to allow STL algorithms to call member functions on this pointer
  /// This is similar to mem_fun_t from the STL standard
  template < typename Ret >
    class mem_fun_t
    {
    public: // typedefs
      /// @c argument_type is the type of the argument
      typedef TYPE argument_type;
      /// @c result_type is the return type
      typedef Ret result_type;

    public: // functions
      /// constructor
      explicit mem_fun_t( Ret ( TYPE::*pf )() ) : m_f(pf) {}

      /// this is the operator we need to add in order for STL algorithms to function
      Ret operator()(const SafePtr<TYPE>& p) const  { return (p.m_ptr->*m_f)(); }

    private: // data
      /// storage of the pointer to the member function to be called
      Ret (TYPE::*m_f)();
    };

}; // end class SafePtr

//////////////////////////////////////////////////////////////////////////////

/// Helper function to allow easy use of the mem_func_t class
/// This is similar to mem_fun from the STL standard
template < typename TYPE , typename Ret >
  inline typename Common::SafePtr<TYPE>::template mem_fun_t<Ret>
  safeptr_mem_fun( Ret (TYPE::*f)() )
  {
    return typename Common::SafePtr<TYPE>::template mem_fun_t<Ret> (f);
  }

//////////////////////////////////////////////////////////////////////////////

template< typename TYPE>
SafePtr<TYPE>::SafePtr() : m_ptr(CFNULL) {}

//////////////////////////////////////////////////////////////////////////////

template< typename TYPE>
SafePtr<TYPE>::SafePtr(TYPE* ptr) {  m_ptr = ptr; }

//////////////////////////////////////////////////////////////////////////////

template< typename TYPE>
SafePtr<TYPE>::~SafePtr() {}

//////////////////////////////////////////////////////////////////////////////

template< typename TYPE>
SafePtr<TYPE>::SafePtr(const SafePtr& other) { m_ptr = other.m_ptr; }

//////////////////////////////////////////////////////////////////////////////

template< typename TYPE>
void SafePtr<TYPE>::reset(TYPE* other) { m_ptr = other; }

//////////////////////////////////////////////////////////////////////////////

template< typename TYPE>
void SafePtr<TYPE>::reset(const SafePtr& other) {  m_ptr = other.m_ptr; }

//////////////////////////////////////////////////////////////////////////////

template< typename TYPE>
inline bool SafePtr<TYPE>::isNull() const { return (m_ptr == CFNULL); }

//////////////////////////////////////////////////////////////////////////////

template< typename TYPE>
inline bool SafePtr<TYPE>::isNotNull() const { return !isNull(); }

//////////////////////////////////////////////////////////////////////////////

template< typename TYPE>
inline const SafePtr<TYPE>& SafePtr<TYPE>::operator= (const SafePtr& other)
{
  reset(other);
  return *this;
}

//////////////////////////////////////////////////////////////////////////////

template< typename TYPE>
inline const SafePtr<TYPE>& SafePtr<TYPE>::operator= (TYPE* other)
{
  reset(other);
  return *this;
}

//////////////////////////////////////////////////////////////////////////////

template< typename TYPE>
inline bool SafePtr<TYPE>::operator== (const SafePtr& other) const
{
  return (m_ptr == other.m_ptr);
}

//////////////////////////////////////////////////////////////////////////////

template< typename TYPE>
inline bool SafePtr<TYPE>::operator!= (const SafePtr& other) const
{
  return (m_ptr != other.m_ptr);
}

//////////////////////////////////////////////////////////////////////////////

template< typename TYPE>
inline bool SafePtr<TYPE>::operator< (const SafePtr& other) const
{
  return (m_ptr < other.m_ptr);
}

//////////////////////////////////////////////////////////////////////////////

template< typename TYPE>
inline bool SafePtr<TYPE>::operator> (const SafePtr& other) const
{
  return (m_ptr > other.m_ptr);
}

//////////////////////////////////////////////////////////////////////////////

template< typename TYPE>
template< typename DTYPE>
SafePtr<DTYPE> SafePtr<TYPE>::d_castTo() const
{
  DTYPE* rPtr = dynamic_cast<DTYPE*>(m_ptr);
  if (rPtr == CFNULL)
  {
    std::string msg ("SafePtr failed dynamic cast from ");
    msg += DEMANGLED_TYPEID(TYPE);
    msg += " to ";
    msg += DEMANGLED_TYPEID(DTYPE);
    throw Common::FailedCastException (FromHere(), msg);
  }
  return SafePtr<DTYPE>(rPtr);
}

//////////////////////////////////////////////////////////////////////////////

} // namespace Common
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif //COOLFluiD_Common_SafePtr_hh
