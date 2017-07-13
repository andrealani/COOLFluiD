// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Common_GeneralStorage_hh
#define COOLFluiD_Common_GeneralStorage_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/COOLFluiD.hh"
#include "Common/NoSuchStorageException.hh"
#include "Common/StorageExistsException.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace Common {

//////////////////////////////////////////////////////////////////////////////

/// GeneralStorage holds and manages the storage of pointers by name
/// @author Tiago Quintino
template <typename TYPE>
class GeneralStorage {

protected: // typedefs

  /// typedef for the inner storage
  typedef typename std::map<std::string,TYPE*> StorageType;

public: // typedefs

  /// typedef for the elements being stored
  typedef TYPE* element_type;
  /// typedef for the iterator to inner storage
  typedef typename StorageType::iterator iterator;
  /// typedef for the iterator to inner storage
  typedef typename StorageType::const_iterator const_iterator;
  /// typedef for the element type stored in the storage
  typedef typename StorageType::value_type value_type;

public: // functions

  /// Default minimal constructor.
  GeneralStorage() : m_storage() {}

  /// Destructor.
  ~GeneralStorage();

  /// @return the number of entries currently in the storage
  size_t size() const { return m_storage.size();  }
  /// @return the iterator to the beginning of the storage
  iterator begin()      { return m_storage.begin(); }
  /// @return the iterator pointing to the end of the storage
  iterator end()        { return m_storage.end();   }
  /// @return the iterator to the beginning of the storage
  const_iterator begin() const { return m_storage.begin(); }
  /// @return the iterator pointing to the end of the storage
  const_iterator end() const   { return m_storage.end();   }

  /// Adds a storage.
  /// @param name std::string identifier for the storage
  /// @param ptr pointer to the storage
  /// @throw Common::StorageExistsException when entry with the same name already exists
  void addEntry(const std::string& name, TYPE* ptr)
  {
    if (!checkEntry(name)) {
      addPtr(name,ptr);
    }
    else {
      throw Common::StorageExistsException (FromHere(),name + " already exists");
    }
  }

  /// Removes a storage.
  /// @param name std::string identifier for the storage
  /// @throw Common::NoSuchStorageException when entry with the same name doesnt exists
  void removeEntry(const std::string& name)
  {
    TYPE * ptr = getPtr(name);
    if (ptr != CFNULL)
    {
      erasePtr(name);
    }
    else {
      throw Common::NoSuchStorageException (FromHere()," didn't find " + name);
    }
  }

  /// Removes a storage and frees its memory.
  /// @param name std::string identifier for the storage
  /// @throw Common::NoSuchStorageException when entry with the same name doesnt exists
  void deleteEntry(const std::string& name)
  {
    TYPE * ptr = getPtr(name);
    if (ptr != CFNULL)
    {
      erasePtr(name);
      deletePtr(ptr);
    }
    else {
      throw Common::NoSuchStorageException (FromHere()," didn't find " + name);
    }
  }

  /// Gets a storage.
  /// @param name std::string identifier for the storage
  /// @return the pointer to the storage
  /// @throw Common::NoSuchStorageException when entry with the same name doesnt exists
  TYPE * getEntry(const std::string& name) const
  {
    TYPE * ptr = getPtr(name);
    if (ptr == CFNULL) {
      throw Common::NoSuchStorageException (FromHere()," didn't find " + name);
    }
    return ptr;
  }

  /// Removes all the entries and deletes the pointers
  /// @post the storage is empty
  void deleteAllEntries()
  {
    std::vector<std::string> deleteList;

    iterator itr = m_storage.begin();
    for(; itr != m_storage.end(); ++itr) {
      deleteList.push_back(itr->first);
    }
    std::vector<std::string>::iterator jtr = deleteList.begin();
    for(; jtr != deleteList.end(); ++jtr) {
      deleteEntry(*jtr);
    }
  }

  /// Removes all the entries
  /// @post the storage is empty
  void removeAllEntries()
  {
    std::vector<std::string> rmList;

    iterator itr = m_storage.begin();
    for(; itr != m_storage.end(); ++itr) {
      rmList.push_back(itr->first);
    }
    std::vector<std::string>::iterator jtr = rmList.begin();
    for(; jtr != rmList.end(); ++jtr) {
      removeEntry(*jtr);
    }
  }

  /// Checks if a storage exists.
  /// @param name std::string identifier for the storage
  /// @return true if exists
  inline bool checkEntry(const std::string& name) const
  {
    return (m_storage.count(name) > 0);
  }

  /// Extracts the entry from the value_type stored.
  /// Usefull for avoiding exposure of the concrete value_type to the client code.
  /// Function is static to be used with STL algorithms
  /// @param entry value_type entry in the storage to extract the actual value
  /// @return the value stored
  static TYPE* extract( const value_type& entry ) { return entry.second; }

private: // helper methods

  /// Access the pointer in storage by its name
  /// @param name std::string identifier for the storage
  /// @post return CFNULL if pointer does not exis
  TYPE* getPtr(const std::string& name) const
  {
    const_iterator citr = m_storage.find(name);
    return (citr == m_storage.end()) ? CFNULL : (*citr).second;
  }

  /// Removes the pointer from the storage
  /// @param name std::string identifier for the storage
  void erasePtr(const std::string& name)
  {
    m_storage.erase(name);
  }

  /// Adds the pointer from the storage
  /// @param name std::string identifier for the storage
  void addPtr(const std::string& name, TYPE* ptr)
  {
    m_storage[name] = ptr;
  }

private: // data

  /// map to store the pointers
  StorageType m_storage;

}; // end class DataStorage

//////////////////////////////////////////////////////////////////////////////

template <typename TYPE>
GeneralStorage<TYPE>::~GeneralStorage()
{
  if(!m_storage.empty())
  {
    std::vector<std::string> eraseList;

    iterator itr = m_storage.begin();
    for(; itr != m_storage.end(); ++itr) {
      eraseList.push_back(itr->first);
    }

    std::vector<std::string>::iterator jtr = eraseList.begin();
    for(; jtr != eraseList.end(); ++jtr) {
      erasePtr(*jtr);
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace Common

} // namespace COOLFluiD

#endif // COOLFluiD_Common_GeneralStorage_hh
