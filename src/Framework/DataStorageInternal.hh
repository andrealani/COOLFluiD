// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_DataStorageInternal_hh
#define COOLFluiD_Framework_DataStorageInternal_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/NoSuchStorageException.hh"
#include "Common/StorageExistsException.hh"

#include "Common/CFLog.hh"

#include "Framework/DataHandle.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// DataStorage creates, holds and manages the storage of
/// auxiliary data in the MeshData.
/// It is implemented in a general way, holding pointers to void.
/// This class is a Helper class.
/// DataStorageInternal encapsulates the functions common to all types of
/// DataStorages to avoid code duplication.
/// It is based on the fact that the current DataStorage's only differ in
/// underlying vector type.
/// @author Tiago Quintino
/// @author Dries Kimpe
template < class COMTYPE >
class DataStorageInternal {
public:

  /// Default minimal constructor.
  DataStorageInternal()
  {
  }

  /// Destructor.
  ~DataStorageInternal()
  {
    if(!_dataStorage.empty()) {

      std::vector<std::string> eraseList;

      std::map<std::string,void*>::iterator itr = _dataStorage.begin();
      for(; itr != _dataStorage.end(); ++itr) {
        eraseList.push_back(itr->first);
      }

      std::vector<std::string>::iterator jtr = eraseList.begin();
      for(; jtr != eraseList.end(); ++jtr) {
        removeDataPtr(*jtr);
      }
    }
  }

  /// Creates and Initializes a storage.
  /// @param name std::string identifier for the storage
  /// @param size size of storage
  /// @param init Object of type Type that initializes ALL entries of storage
  /// @return the pointer to the created storage
  template <class TYPE>
  DataHandleInternal< TYPE, COMTYPE >
  createData(const std::string& name, const CFuint size, const TYPE& init = TYPE ());

  template <class TYPE>
  DataHandleInternal< TYPE, COMTYPE >
  createDataDynamic (const std::string & name,    const CFuint size,
                     const size_t ElementSize, const TYPE & Init = TYPE());


  /// Deletes a storage and frees its memory.
  /// @param name std::string identifier for the storage
  /// @return the number of storages deleted (normally 0 or 1)
  template <class TYPE>
  CFuint deleteData(const std::string& name);

  /// Deletes a storage and frees its memory.
  /// @param name std::string identifier for the storage
  /// @return the number of storages deleted (normally 0 or 1)
  template <class TYPE>
  CFuint deleteDataPlusPtr(const std::string& name);

  /// Gets a storage.
  /// @param name std::string identifier for the storage
  /// @return the pointer to the storage, but casted to the correct type
  template <class TYPE>
  DataHandleInternal< TYPE, COMTYPE >
  getData(const std::string& name) ;

  /// Checks if a storage exists.
  /// @param name std::string identifier for the storage
  /// @return true if exists
  inline bool checkData(const std::string& name) ;

private: // helper methods

  /// Places a storage inside.
  /// @param name std::string identifier for the storage
  /// @param store pointer holding the storage
  void setDataPtr(const std::string& name, void* store);

  /// Gets the storage known by the name supplied.
  /// @param name is the name of the storage to get
  /// @return the pointer to the storage
  void* getDataPtr(const std::string& name) ;

  /// Removes the storage pointer by the name supplied.
  /// @param name is the name of the storage to remove
  /// @return the number of storages deleted (normally 0 or 1)
  CFuint removeDataPtr(const std::string& name);

private: // data

  /// map to store the pointers that hold the data
  std::map<std::string,void*> _dataStorage;

}; // end class DataStorageInternal

//////////////////////////////////////////////////////////////////////////////

template < class COMTYPE >
template < class TYPE >
DataHandleInternal< TYPE, COMTYPE >
DataStorageInternal<COMTYPE>::
  createDataDynamic(const std::string & name,    const CFuint size,
                    const size_t ElementSize, const TYPE & init)
{
  typedef typename COMTYPE::template StoragePolicy<TYPE>::ContainerType ContainerType;
  typedef typename COMTYPE::template StoragePolicy<TYPE>::ElemType      ElemType;
  typedef DataHandleInternal< TYPE, COMTYPE > ReturnType;

  if(!checkData(name))
  {
    ContainerType* ptr = new ContainerType(init,size,ElementSize);
    setDataPtr(name,ptr);
    CFLogDebugMin("Created Storage(dynamic): " << name << "\n");
    return ReturnType (ptr);
  }
  else
  {
    throw Common::StorageExistsException (FromHere(),name + " already exists.");
  }
}

//////////////////////////////////////////////////////////////////////////////

template < class COMTYPE >
template < class TYPE >
DataHandleInternal< TYPE, COMTYPE >
DataStorageInternal<COMTYPE>::
  createData(const std::string& name, const CFuint size, const TYPE& init)
{
  typedef typename COMTYPE::template StoragePolicy<TYPE>::ContainerType ContainerType;
  typedef typename COMTYPE::template StoragePolicy<TYPE>::ElemType      ElemType;
  typedef DataHandleInternal< TYPE, COMTYPE > ReturnType;

  if(!checkData(name))
  {
    ContainerType* ptr = new ContainerType(init,size);
    setDataPtr(name,ptr);
    CFLogDebugMin("Created Storage: " << name << "\n");
    return ReturnType(ptr);
  }
  else
  {
    throw Common::StorageExistsException (FromHere(),name + " already exists.");
  }
}

//////////////////////////////////////////////////////////////////////////////

template < class COMTYPE >
template < class TYPE >
CFuint DataStorageInternal<COMTYPE>::deleteData(const std::string& name)
{
  typedef typename COMTYPE::template StoragePolicy<TYPE>::ContainerType ContainerType;
  typedef typename COMTYPE::template StoragePolicy<TYPE>::ElemType      ElemType;

  ContainerType* ptr = static_cast<ContainerType*>(getDataPtr(name));
  if (ptr == CFNULL)
  {
    throw Common::NoSuchStorageException (FromHere()," didn't find " + name);
  }
  delete ptr;
  CFLogDebugMin("Deleted Storage: " << name << "\n");
  return removeDataPtr(name);
}

//////////////////////////////////////////////////////////////////////////////

template < class COMTYPE >
template < class TYPE >
CFuint DataStorageInternal<COMTYPE>::deleteDataPlusPtr(const std::string& name)
{
  typedef typename COMTYPE::template StoragePolicy<TYPE>::ContainerType ContainerType;
  typedef typename COMTYPE::template StoragePolicy<TYPE>::ElemType      ElemType;

  ContainerType* ptr = static_cast<ContainerType*>(getDataPtr(name));
  if (ptr == CFNULL)
  {
    throw Common::NoSuchStorageException (FromHere()," didn't find " + name);
  }
  for (CFuint i = 0; i < ptr->size(); ++i)
  {
    delete (*ptr)[i];
  }
  delete ptr;
  CFLogDebugMin("Deleted Storage: " << name << "\n");
  return removeDataPtr(name);
}

//////////////////////////////////////////////////////////////////////////////

template < class COMTYPE >
template < class TYPE >
DataHandleInternal< TYPE, COMTYPE > DataStorageInternal<COMTYPE>::getData(const std::string& name)
{
  typedef typename COMTYPE::template StoragePolicy<TYPE>::ContainerType ContainerType;
  typedef typename COMTYPE::template StoragePolicy<TYPE>::ElemType      ElemType;
  typedef DataHandleInternal< TYPE, COMTYPE > ReturnType;

  ContainerType* ptr = static_cast<ContainerType* >(getDataPtr(name));
  if (ptr == CFNULL)
  {
    throw Common::NoSuchStorageException (FromHere()," didn't find " + name);
  }
  CFLogDebugMax("Got Storage: " << name << "\n");
  return ReturnType(ptr);
}

//////////////////////////////////////////////////////////////////////////////

template < class COMTYPE >
inline bool DataStorageInternal<COMTYPE>::checkData(const std::string& name)
{
  return (_dataStorage.count(name) > 0);
}

//////////////////////////////////////////////////////////////////////////////

template < class COMTYPE >
void DataStorageInternal<COMTYPE>::setDataPtr(const std::string& name, void* store)
{
  cf_assert(_dataStorage.count(name) == 0);
  _dataStorage[name] = store;
}

//////////////////////////////////////////////////////////////////////////////

template < class COMTYPE >
void* DataStorageInternal<COMTYPE>::getDataPtr(const std::string& name)
{
  return (_dataStorage.count(name) == 0) ? CFNULL : _dataStorage[name];
}

//////////////////////////////////////////////////////////////////////////////

template < class COMTYPE >
CFuint DataStorageInternal<COMTYPE>::removeDataPtr(const std::string& name)
{
  return static_cast<CFuint>(_dataStorage.erase(name));
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_DataStorageInternal_HH
