// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_DataStorage_hh
#define COOLFluiD_Framework_DataStorage_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/DataHandle.hh"
#include "Framework/LocalCommTypes.hh"
#include "Framework/GlobalCommTypes.hh"
#include "Framework/GlobalTypeTrait.hh"

#include "Common/NoSuchStorageException.hh"
#include "Common/StorageExistsException.hh"

#include "Common/CFLog.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// Class DataStorage
/// This is class provides a default using the LOCAL communicator type
/// @author Tiago Quintino
/// @author Andrea Lani
/// @author Dries Kimpe
class Framework_API DataStorage {
public:

  /// Default minimal constructor.
  DataStorage();

  /// Destructor.
  ~DataStorage();

  /// Gets a storage.
  /// @param name std::string identifier for the storage
  /// @return the pointer to the storage, but casted to the correct type
  template <class TYPE>
  DataHandle<TYPE, LOCAL> getData(const std::string& name) ;

  /// Creates and Initializes a storage.
  /// @param name std::string identifier for the storage
  /// @param size size of storage
  /// @param init Object of type Type that initializes ALL entries of storage
  /// @return the pointer to the created storage
  template <class TYPE>
  DataHandle<TYPE, GLOBAL> getGlobalData(const std::string& name);

  /// Creates and Initializes a storage.
  /// @param name std::string identifier for the storage
  /// @param size size of storage
  /// @param init Object of type Type that initializes ALL entries of storage
  /// @return the pointer to the created storage
  template <class TYPE>
  DataHandle< TYPE, LOCAL> createData (const std::string& name, CFuint size, const TYPE& init = TYPE ());

  /// Creates and Initializes a storage.
  /// @param name std::string identifier for the storage
  /// @param size size of storage
  /// @param init Object of type Type that initializes ALL entries of storage
  /// @return the pointer to the created storage
  template <class TYPE> DataHandle<TYPE, GLOBAL> createGlobalData
    (const std::string & name, const std::string& nspaceName, CFuint size,
     const typename Framework::GlobalTypeTrait<TYPE>::GTYPE& init);
  
  /// Creates and initializes local data dynamically
  /// @param name std::string identifier for the storage
  /// @param size size of storage
  /// @param init Object of type Type that initializes ALL entries of storage
  /// @return the pointer to the created storage
  template <class TYPE>
  DataHandle< TYPE, LOCAL> createDataDynamic (const std::string & name, CFuint size, CFuint elementSize, const TYPE& init = TYPE());

  /// Creates and initializes local data dynamically
  /// @param name std::string identifier for the storage
  /// @param size size of storage
  /// @param init Object of type Type that initializes ALL entries of storage
  /// @return the pointer to the created storage
  template <class TYPE>
  DataHandle<TYPE, GLOBAL> createGlobalDataDynamic (const std::string & name, CFuint size, CFuint elementSize, const typename Framework::GlobalTypeTrait<TYPE>::GTYPE & init);

  /// Deletes a storage and frees its memory.
  /// @param name std::string identifier for the storage
  /// @return the number of storages deleted (normally 0 or 1)
  template <class TYPE>
  CFuint deleteDataPlusPtr(const std::string& name);

  /// Deletes a storage and frees its memory.
  /// @param name std::string identifier for the storage
  /// @return the number of storages deleted (normally 0 or 1)
  template <class TYPE>
  CFuint deleteData(const std::string& name);

  /// Deletes a storage and frees its memory.
  /// @param name std::string identifier for the storage
  /// @return the number of storages deleted (normally 0 or 1)
  template <class TYPE>
  CFuint deleteGlobalData(const std::string& name);

  /// Checks if a storage exists.
  /// @param name std::string identifier for the storage
  /// @return true if exists
  bool checkData(const std::string& name);
  
  /// Dumps the contents of the DataStorage to a string
  std::string dump () const;

private:

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

  typedef std::map<std::string,void*> MapType;

  /// map to store the pointers that hold the data
  MapType m_dataStorage;

}; // end class DataStorageInternal

//////////////////////////////////////////////////////////////////////////////

template <class TYPE>
DataHandle<TYPE, LOCAL> DataStorage::getData(const std::string& name)
{
  typedef typename LOCAL::template StoragePolicy<TYPE>::ContainerType ContainerType;
  // typedef typename LOCAL::template StoragePolicy<TYPE>::ElemType      ElemType;
  typedef DataHandle< TYPE, LOCAL> ReturnType;

  ContainerType* ptr = static_cast<ContainerType* >(getDataPtr(name));
  if (ptr == CFNULL)
  {
    throw Common::NoSuchStorageException (FromHere()," didn't find " + name);
  }
  CFLogDebugMax("Got Storage: " << name << "\n");
  return ReturnType(ptr);
}

//////////////////////////////////////////////////////////////////////////////

template < class TYPE >
DataHandle< TYPE, LOCAL> DataStorage::createDataDynamic
(const std::string & name, CFuint size,
 CFuint elementSize, const TYPE & init)
{
  typedef typename LOCAL::template StoragePolicy<TYPE>::ContainerType ContainerType;
  // typedef typename LOCAL::template StoragePolicy<TYPE>::ElemType      ElemType;
  typedef DataHandle< TYPE, LOCAL> ReturnType;

  if(!checkData(name))
  {
    ContainerType* ptr = new ContainerType(init,size,elementSize);
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

template <class TYPE>
DataHandle< TYPE, LOCAL> DataStorage::createData
(const std::string& name, CFuint size, const TYPE& init)
{
  typedef typename LOCAL::template StoragePolicy<TYPE>::ContainerType ContainerType;
  // typedef typename LOCAL::template StoragePolicy<TYPE>::ElemType      ElemType;
  typedef DataHandle<TYPE, LOCAL> ReturnType;

  if(!checkData(name))
  {
    ContainerType* ptr = new ContainerType(init,size); // this will not work with std::vector(size,init)
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

template < class TYPE >
CFuint DataStorage::deleteData(const std::string& name)
{
  typedef typename LOCAL::template StoragePolicy<TYPE>::ContainerType ContainerType;
  // typedef typename LOCAL::template StoragePolicy<TYPE>::ElemType      ElemType;

  ContainerType* ptr = static_cast<ContainerType*>(getDataPtr(name));
  if (ptr == CFNULL)
  {
    throw Common::NoSuchStorageException (FromHere()," didn't find " + name);
  }
  delete ptr; ptr = CFNULL;
  CFLogDebugMin("Deleted Storage: " << name << "\n");
  return removeDataPtr(name);
}

//////////////////////////////////////////////////////////////////////////////

template < class TYPE >
CFuint DataStorage::deleteDataPlusPtr(const std::string& name)
{
  typedef typename LOCAL::template StoragePolicy<TYPE>::ContainerType ContainerType;
  //  typedef typename LOCAL::template StoragePolicy<TYPE>::ElemType      ElemType;

  ContainerType* ptr = static_cast<ContainerType*>(getDataPtr(name));
  if (ptr == CFNULL)
  {
    throw Common::NoSuchStorageException (FromHere()," didn't find " + name);
  }
  for (CFuint i = 0; i < ptr->size(); ++i)
  {
    delete (*ptr)[i]; (*ptr)[i] = CFNULL;
  }
  delete ptr;
  CFLogDebugMin("Deleted Storage: " << name << "\n");
  return removeDataPtr(name);
}

//////////////////////////////////////////////////////////////////////////////

inline bool DataStorage::checkData(const std::string& name)
{
  return (m_dataStorage.count(name) > 0);
}

//////////////////////////////////////////////////////////////////////////////

inline void DataStorage::setDataPtr(const std::string& name, void* store)
{
  cf_assert(m_dataStorage.count(name) == 0);
  m_dataStorage[name] = store;
}

//////////////////////////////////////////////////////////////////////////////

inline void* DataStorage::getDataPtr(const std::string& name)
{
  return (m_dataStorage.count(name) == 0) ? CFNULL : m_dataStorage[name];
}

//////////////////////////////////////////////////////////////////////////////

inline CFuint DataStorage::removeDataPtr(const std::string& name)
{
  return static_cast<CFuint>(m_dataStorage.erase(name));
}

//////////////////////////////////////////////////////////////////////////////

template <class TYPE>
DataHandle<TYPE, GLOBAL> DataStorage::createGlobalDataDynamic (const std::string & name, CFuint size, CFuint elementSize,
const typename Framework::GlobalTypeTrait<TYPE>::GTYPE & init)
{
  typedef typename Framework::GlobalTypeTrait<TYPE>::GTYPE GTYPE;
  typedef typename GLOBAL::template StoragePolicy<GTYPE>::ContainerType GlobalVectorType;
  typedef typename LOCAL::template StoragePolicy<TYPE>::ContainerType LocalVectorType;

  MapType::const_iterator iter = m_dataStorage.find (name);
  if (iter != m_dataStorage.end ())
    throw Common::StorageExistsException
      (FromHere(),name + " already exists (createDataDynamic)!");

  LocalVectorType* vLocal = new LocalVectorType(CFNULL, size, elementSize);
  const std::string localname  = name + "_local";
  m_dataStorage[localname] = static_cast<void *>(vLocal);

  GlobalVectorType* vGlobal = new GlobalVectorType (init, size, elementSize);
  const std::string globalname = name + "_global";
  m_dataStorage[globalname] = static_cast<void *>(vGlobal);

  return DataHandle<TYPE,GLOBAL>(static_cast<void *>(vLocal),static_cast<void *>(vGlobal));
}

//////////////////////////////////////////////////////////////////////////////

template <class TYPE>
DataHandle<TYPE,GLOBAL> DataStorage::getGlobalData (const std::string & name)
{
  // local storage whose entries type are proxies for the "raw" gobal storage
  const std::string localname = name + "_local";
  MapType::const_iterator iterL = m_dataStorage.find (localname);

  if (iterL == m_dataStorage.end ()) {
    throw Common::NoSuchStorageException (FromHere(),"Storage " + localname + " not found!");
  }


  // global storage whose entries type is suitable for parallel communication
  const std::string globalname = name + "_global";

  MapType::const_iterator iterG = m_dataStorage.find (globalname);
  if (iterG == m_dataStorage.end ()) {
    throw Common::NoSuchStorageException (FromHere(),"Storage " + globalname + " not found!");
  }

  return DataHandle<TYPE,GLOBAL>(iterL->second, iterG->second);
}

//////////////////////////////////////////////////////////////////////////////

template <class TYPE>
DataHandle<TYPE,GLOBAL> DataStorage::createGlobalData
(const std::string & name, const std::string& nspaceName, 
 CFuint size, const typename Framework::GlobalTypeTrait<TYPE>::GTYPE & init)
{
  typedef typename Framework::GlobalTypeTrait<TYPE>::GTYPE GTYPE;
  typedef typename GLOBAL::template StoragePolicy<GTYPE>::ContainerType GlobalVectorType;
  typedef typename LOCAL::template StoragePolicy<TYPE>::ContainerType LocalVectorType;

  MapType::const_iterator iter = m_dataStorage.find (name);
  if (iter != m_dataStorage.end ())
    throw Common::StorageExistsException
      (FromHere(),name + " already exists (createDataDynamic)!");

  LocalVectorType* vLocal = new LocalVectorType(CFNULL, size);
  const std::string localname  = name + "_local";
  m_dataStorage[localname] = static_cast<void *>(vLocal);
  
  GlobalVectorType* vGlobal = new GlobalVectorType (nspaceName, init, size);
  const std::string globalname = name + "_global";
  m_dataStorage[globalname] = static_cast<void *>(vGlobal);
  
  return DataHandle<TYPE,GLOBAL>(static_cast<void *>(vLocal),static_cast<void *>(vGlobal));
}

//////////////////////////////////////////////////////////////////////////////

template <class TYPE>
CFuint DataStorage::deleteGlobalData (const std::string & name)
{
  typedef typename Framework::GlobalTypeTrait<TYPE>::GTYPE GTYPE;
  typedef typename GLOBAL::template StoragePolicy<GTYPE>::ContainerType GlobalVectorType;
  typedef typename LOCAL::template StoragePolicy<TYPE>::ContainerType LocalVectorType;

  // local storage whose entries type are proxies for the "raw" gobal storage
  const std::string localname = name + "_local";
  MapType::iterator iterL = m_dataStorage.find (localname);

  if (iterL == m_dataStorage.end ())
    throw Common::NoSuchStorageException (FromHere(),"Storage " + localname + " doesn't exist!");
  delete static_cast<LocalVectorType*>(iterL->second);
  iterL->second = CFNULL;
  removeDataPtr(localname);

  // global storage whose entries type is suitable for parallel communication
  const std::string globalname = name + "_global";
  MapType::iterator iterG = m_dataStorage.find (globalname);
  if (iterG == m_dataStorage.end ())
    throw Common::NoSuchStorageException (FromHere(),"Storage " + globalname + " doesn't exist!");

  delete static_cast<GlobalVectorType*>(iterG->second);
  iterG->second = CFNULL;
  return removeDataPtr(globalname);
}

//////////////////////////////////////////////////////////////////////////////

}  //  namespace Framework
}  //  namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_DataStorage_hh
