// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_DataHandleInternal_hh
#define COOLFluiD_Framework_DataHandleInternal_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/SafePtr.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This class is the internal representation of Array handle
/// to work with DataStorage facility.
/// This is a template class, where the template represents the class of
/// objects that are store in the array.
/// @see DataStorage
/// @author Tiago Quintino
/// @author Dries Kimpe
/// @author Andrea Lani
template < typename TYPE , typename COMTYPE >
class DataHandleInternal {

public: // typedefs

  typedef typename COMTYPE::template StoragePolicy<TYPE>::ContainerType StorageType;
  typedef typename COMTYPE::template StoragePolicy<TYPE>::ElemType      ElemType;

public: // functions

  /// Constructor.
  explicit DataHandleInternal(StorageType *const ptr) : _ptr(ptr)  {}

  /// Copy constructor.
  DataHandleInternal (const DataHandleInternal & other) : _ptr(other._ptr) {}

  /// Overloading of the operator =
  const DataHandleInternal & operator= (const DataHandleInternal & other)
  {
    _ptr = other._ptr;
    return *this;
  }

  /// Default destructor
  ~DataHandleInternal () {}

  /// Operator"[]" for assignment.
  /// @pre idx less than size
  /// @param idx missing documentation
  /// @return missing documentation
  ElemType& operator[] (const CFuint& idx)
  {
    cf_assert(_ptr != CFNULL);
    cf_assert(idx < _ptr->size());
    return (*_ptr)[idx];
  }

  /// Operator"[]" for acessing.
  /// @pre idx less than size
  /// @param idx missing documentation
  /// @return missing documentation
  const ElemType& operator[] (const CFuint& idx) const
  {
    cf_assert(_ptr != CFNULL);
    cf_assert(idx < _ptr->size());
    return (*_ptr)[idx];
  }

  /// Overloading of the operator ()
  /// param i index row
  /// param j index column
  /// @return a reference to the element i,j in a
  ///         one dimensional array of size == sizeI*sizeJ
  /// @pre  i*sizeN+j less than size
  ElemType& operator() (const CFuint& i,
      const CFuint& j,
      const CFuint& sizeJ)
  {
    const CFuint idx = i*sizeJ+j;
  return (*this)[idx];
  }

  /// Overloading of the operator ()
  /// param i index row
  /// param j index column
  /// @return a const reference to the element i,j in a
  ///         one dimensional array of size == sizeI*sizeJ
  /// @pre  i*sizeN+j less than size
  const ElemType& operator() (const CFuint& i,
    const CFuint& j,
    const CFuint& sizeJ) const
  {
    const CFuint idx = i*sizeJ+j;
    cf_assert(_ptr != CFNULL);
    cf_assert(idx < _ptr->size());
    return (*_ptr)[idx];
  }

  /// Operator"=" for assigning to ALL members.
  /// @param value missing documentation
  /// @return missing documentation
  const DataHandleInternal & operator= (const ElemType& value)
  {
    cf_assert(_ptr != CFNULL);
    const CFuint size = _ptr->size();
    for (CFuint i = 0; i < size; i++) {
      (*_ptr)[i] = value;
    }
    return *this;
  }

  /// Missing documentation
  /// @return the allocated size of this Array
  CFuint size() const
  {
    cf_assert(_ptr != CFNULL);
    return _ptr->size();
  }

  /// Resize the allocated size of this Array
  /// @param size missing documentation
  void resize(const CFuint size) const
  {
    cf_assert(_ptr != CFNULL);
    _ptr->resize(size);
  }
  
  /// get a pointer to the global array
  Common::SafePtr<StorageType> getLocalArray() const {return _ptr;}
  
protected: // data

  /// the pointer to be handled.
  StorageType* _ptr;

}; // end of class DataHandleInternal

//////////////////////////////////////////////////////////////////////////////

} // namespace Framework
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_DataHandleInternal_hh
