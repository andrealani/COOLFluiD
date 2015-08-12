// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_DataHandle_hh
#define COOLFluiD_Framework_DataHandle_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/LocalCommTypes.hh"
#include "Framework/DataHandleInternal.hh"
#include "Common/ShouldNotBeHereException.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This is the default implementation of a DataHandle.
/// By default it uses the LOCAL communication type.
/// @author Dries Kimpe
/// @author Tiago Quintino
template < class TYPE , class COMTYPE = LOCAL >
class DataHandle : public  DataHandleInternal < TYPE , COMTYPE >
{

public:

    typedef DataHandleInternal < TYPE , COMTYPE > BaseClass;
    typedef typename BaseClass::StorageType       StorageType;
    typedef typename BaseClass::ElemType          ElemType;

public:

  /// Constructor taking the storage type
  explicit DataHandle (StorageType *const ptr) : BaseClass(ptr) {}

  /// Constructor taking the actual storage types
  DataHandle (StorageType *const ptr, StorageType *const global_ptr) : BaseClass(ptr) {}

  /// Constructor taking void*
  DataHandle (void* ptr, void* global_ptr) : BaseClass ( static_cast<StorageType*>(ptr) ) {}

  /// Copy constructor
  DataHandle (const DataHandle<TYPE, COMTYPE> & t) : BaseClass (t) {}

  /// Assignment operator
  /// This is just passed on to the underlying DataHandleInternal
  DataHandle<TYPE,COMTYPE>& operator = (const TYPE t) { BaseClass::operator = (t); return *this; }

  /// This function returns the global (cross-processes) size of
  /// the underlying parallel array
  /// @return the global size of the parallel array
  CFuint getGlobalSize() const 
  { return (this->_ptr != CFNULL) ? this->_ptr->size() : 0; }

  /// This function returns the local size of the underlying parallel array.
  /// Return the same as getGlobalSize, as this has a pure COMTYPE vector.
  /// @return the local size of the parallel array
  CFuint getLocalSize() const
  { return (this->_ptr != CFNULL) ? this->_ptr->size() : 0; }
    
  /// This does nothing on a local datahandle
  void beginSync ()  {}
  
  /// This does nothing on a local datahandle
  void endSync () {}

  /// This does nothing on a local datahandle
  void DumpContents () {}

  /// This does nothing on a local datahandle
  void DumpInfo () {}

  /// Adding ghost points on non-parallel vectors is not allowed
  unsigned int addGhostPoint (unsigned int GlobalID)
  {
    throw Common::ShouldNotBeHereException
      (FromHere(),"This COMTYPE does not know how to adding ghost points to storage");
    return 0;
  }

  /// Returns the list of ghost nodes (by processor rank) to be sent to
  /// another processor
  const std::vector< std::vector< unsigned int > >&
    getGhostSendList() const
  {
    throw Common::ShouldNotBeHereException
      (FromHere(),"This COMTYPE does not know how to send ghost points list");
  }

  /// Returns the list of ghost nodes (by processor rank) to be received
  /// from another processor
  const std::vector< std::vector< unsigned int > >&
    getGhostReceiveList() const
  {
    throw Common::ShouldNotBeHereException
      (FromHere(),"This COMTYPE does not know how to ask ghost points list");
  }
  
  /// Add local point and return the index of the newly created.
  unsigned int addLocalPoint (unsigned int GlobalID)  { return BaseClass::_ptr->increase (); }

  /// Make sure there is at least size capacity
  void reserve (unsigned int S)  {   BaseClass::_ptr->reserve (S);  }

  /// buildMap does nothing on a pure local vector
  void buildMap () {}

}; // end class DataHandle

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_DataHandle_hh
