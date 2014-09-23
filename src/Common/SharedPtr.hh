// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Common_SharedPtr_hh
#define COOLFluiD_Common_SharedPtr_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/SafePtr.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
namespace Common {

//////////////////////////////////////////////////////////////////////////////

/// Class that defines a smart pointer that allows to share objects (but ONLY
/// objects deriving from OwnedObject) using reference counting.
/// @see OwnedObject
/// @author Andrea Lani
/// @author Tiago Quintino
template<class T>
class SharedPtr {
public:

  /// Constructor.
  SharedPtr ();

  /// Constructor.
  /// @param ptr naked pointer
  explicit SharedPtr ( T* ptr );

  /// Copy constructor.
  /// @param other missing documentation
  SharedPtr ( const SharedPtr& other );

  /// Copy constructor taking a SafePtr
  /// @param other missing documentation
  SharedPtr ( const SafePtr<T>& other );

  /// Destructor.
  ~SharedPtr();

  /// Release ownership
  /// @post m_ptr = CFNULL
  void release();

  /// Reset the ptr
  void reset(T* other);

  /// Reset the ptr
  /// @post m_ptr = other
  void reset(const SharedPtr& other);

  /// @returns true if m_ptr == CFNULL
  inline bool isNull() const { return (m_ptr == CFNULL); }

  /// @returns true if m_ptr != CFNULL
  inline bool isNotNull() const { return (m_ptr != CFNULL); }

  /// Overloading of "=" with SharedPtr
  /// @param other missing documentation
  /// @return missing documentation
  inline const SharedPtr& operator= (const SharedPtr& other);

  /// Overloading of "=" with SharedPtr
  /// @param other missing documentation
  /// @return missing documentation
  inline const SharedPtr& operator= (T* other);

  /// Overloading of "=="
  /// @param other missing documentation
  /// @return missing documentation
  inline bool operator== (const SharedPtr& other);

  /// @returns the naked pointer to the object
  T* getPtr() const { return m_ptr; }

  /// Overloading of "->"
  /// @returns the naked pointer to the object
  T* operator->() const
  {
    cf_assert(m_ptr != CFNULL);
    return m_ptr;
  }

  /// Overloading of "*"
  /// @return missing documentation
  T& operator*() const
  {
    cf_assert(m_ptr != CFNULL);
    return *m_ptr;
  }

private:

  /// Raw pointer to the object.
  T* m_ptr;

}; // end class SharedPtr

//////////////////////////////////////////////////////////////////////////////

template<class T>
SharedPtr<T>::SharedPtr() : m_ptr(0) {}

//////////////////////////////////////////////////////////////////////////////

template<class T>
SharedPtr<T>::SharedPtr(T* ptr)
{
  m_ptr = ptr;
  if (m_ptr != CFNULL) {
    m_ptr->addOwner();
  }
}

//////////////////////////////////////////////////////////////////////////////

template<class T>
SharedPtr<T>::~SharedPtr()
{
  if(m_ptr != CFNULL) {
    m_ptr->removeOwner();
    if(m_ptr->hasNoOwner()) {
      deletePtr(m_ptr);
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

template<class T>
SharedPtr<T>::SharedPtr(const SharedPtr& other) : m_ptr(other.m_ptr)
{
  if (m_ptr != CFNULL) {
    m_ptr->addOwner();
  }
}

//////////////////////////////////////////////////////////////////////////////

template<class T>
SharedPtr<T>::SharedPtr(const SafePtr<T>& other) : m_ptr(other.m_ptr)
{
  if (m_ptr != CFNULL)
    m_ptr->addOwner();
}

//////////////////////////////////////////////////////////////////////////////

template<class T>
void SharedPtr<T>::release()
{
  if (m_ptr != CFNULL) {
    m_ptr->removeOwner();
    if(m_ptr->hasNoOwner()) {
      deletePtr(m_ptr);
    }
  }
  m_ptr = CFNULL;
}

//////////////////////////////////////////////////////////////////////////////

template<class T>
void SharedPtr<T>::reset(T* other)
{
  if (m_ptr != CFNULL) {
    m_ptr->removeOwner();
    if(m_ptr->hasNoOwner()) {
      deletePtr(m_ptr);
    }
  }
  m_ptr = other;
  if (m_ptr != CFNULL) {
    m_ptr->addOwner();
  }
}

//////////////////////////////////////////////////////////////////////////////

template<class T>
void SharedPtr<T>::reset(const SharedPtr& other)
{
  if (m_ptr != CFNULL) {
    m_ptr->removeOwner();
    if(m_ptr->hasNoOwner()) {
      deletePtr(m_ptr);
    }
  }
  m_ptr = other.m_ptr;
  if (m_ptr != CFNULL) {
    m_ptr->addOwner();
  }
}

//////////////////////////////////////////////////////////////////////////////

template<class T>
inline const SharedPtr<T>& SharedPtr<T>::operator= (const SharedPtr& other)
{
  reset(other);
  return *this;
}

//////////////////////////////////////////////////////////////////////////////

template<class T>
inline const SharedPtr<T>& SharedPtr<T>::operator= (T* other)
{
  reset(other);
  return *this;
}

//////////////////////////////////////////////////////////////////////////////

template<class T>
inline bool SharedPtr<T>::operator== (const SharedPtr& other)
{
  return (m_ptr == other.m_ptr);
}

//////////////////////////////////////////////////////////////////////////////

} // namespace Common
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif //COOLFluiD_Common_SharedPtr_hh
