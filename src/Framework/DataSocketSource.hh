// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_DataSocketSource_hh
#define COOLFluiD_Framework_DataSocketSource_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/NonCopyable.hh"

#include "Framework/BaseDataSocketSource.hh"
#include "Framework/DataSocketHelper.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This class is the base class for objects that need a DataSocket of a certain TYPE
/// @author Tiago Quintino
/// @author Andrea Lani
template < typename TYPE , typename STORAGETYPE = LOCAL>
class DataSocketSource : public BaseDataSocketSource {

public:

  /// Constructor
  DataSocketSource(const std::string& name);

  /// Default destructor.
  virtual ~DataSocketSource();

  /// Returns the DataHandle associated with the Socket
  DataHandle<TYPE,STORAGETYPE> getDataHandle() const { cf_assert ( isAllocated() ) ; return m_handle;  }
  
  /// Allocation of the DataHandle in this DataSocket
  void allocate (Common::SafePtr<DataStorage> storage, const std::string& nspaceName);
  
  /// Deallocation of the DataHandle in this DataSocket
  void deallocate ();
  
  /// Checks if this socket has been allocated
  bool isAllocated () const { return m_storage.isNotNull(); }
  
  /// @return the global size of the underlying data array
  virtual CFuint getGlobalSize() const {return getDataHandle().getGlobalSize();}
  
  /// @return the local size of the underlying data array
  virtual CFuint getLocalSize() const {return getDataHandle().getLocalSize();}
  
private:
  
  /// access to the datahandle on the datastorage
  DataHandle<TYPE,STORAGETYPE> m_handle;

  /// Acquaintance of the DataStorage from which the DataHandle was allocated
  Common::SafePtr<DataStorage> m_storage;

}; // end of class DataSocketSource

//////////////////////////////////////////////////////////////////////////////

template < typename TYPE, typename STORAGETYPE >
void deleteAllPtr(DataSocketSource<TYPE,STORAGETYPE>& socket)
{
}

//////////////////////////////////////////////////////////////////////////////

template < typename TYPE, typename STORAGETYPE >
void deleteAllPtr(DataSocketSource<TYPE*,STORAGETYPE>& socket)
{
  CFLogDebugMin ( "Deleting all data in socket " << socket.getDataSocketName() << "\n" );

  DataHandle<TYPE*,STORAGETYPE> handle = socket.getDataHandle();
  const size_t size = handle.size();
  for ( CFuint i = 0; i < size; ++i )
  {
    delete handle[i];
    handle[i] = CFNULL;
  }
}

//////////////////////////////////////////////////////////////////////////////

template < typename TYPE , typename STORAGETYPE >
DataSocketSource<TYPE,STORAGETYPE>::
DataSocketSource(const std::string& name) :
  BaseDataSocketSource(name,DEMANGLED_TYPEID(STORAGETYPE),DEMANGLED_TYPEID(TYPE)),
  m_handle(CFNULL),
  m_storage(CFNULL)
{
  DataBroker::getInstance().registerSource ( this );
}

//////////////////////////////////////////////////////////////////////////////GTYPE

template < typename TYPE , typename STORAGETYPE >
DataSocketSource<TYPE,STORAGETYPE>::~DataSocketSource()
{
  DataBroker::getInstance().unregisterSource ( this );
}

//////////////////////////////////////////////////////////////////////////////

template < typename TYPE , typename STORAGETYPE>
void DataSocketSource<TYPE,STORAGETYPE>::allocate(Common::SafePtr<DataStorage> storage,
						  const std::string& nspaceName)
{
  deallocate();
  m_storage = storage;

  const CFuint MIN_SIZE = 0; // 1 creates memory leaks but 0 is also problematic
  // creates an empty data storage
  CreateDataHandle<TYPE, STORAGETYPE>(m_storage, getDataSocketFullStorageName(), 
				      nspaceName, (CFuint)(MIN_SIZE), m_handle);
}

//////////////////////////////////////////////////////////////////////////////

template < typename TYPE , typename STORAGETYPE >
void DataSocketSource<TYPE,STORAGETYPE>::deallocate()
{
  if ( m_storage.isNotNull() )
    DeleteDataHandle<TYPE, STORAGETYPE>(m_storage, getDataSocketFullStorageName());
}

//////////////////////////////////////////////////////////////////////////////

} // namespace Framework
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_DataSocketSource_hh
