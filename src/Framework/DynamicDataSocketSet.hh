// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_DynamicDataSocketSet_hh
#define COOLFluiD_Framework_DynamicDataSocketSet_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/SafePtr.hh"
#include "Common/GeneralStorage.hh"

#include "Framework/NamespaceMember.hh"
#include "Framework/DataSocketSink.hh"
#include "Framework/DataSocketSource.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a set of sockets that aren't build at
/// construction time, hence their names are only chosen dynamically
/// @author Tiago Quintino
/// @author Andrea Lani
template <typename COMTYPE = LOCAL>
class DynamicDataSocketSet : public NamespaceMember {
public:

  /// storage type for list of sockets sinks
  typedef Common::GeneralStorage<BaseDataSocketSink> StorageTypeSink;
  /// storage type for list of sockets sources
  typedef Common::GeneralStorage<BaseDataSocketSource> StorageTypeSource;

public:

  /// Default constructor
  DynamicDataSocketSet();

  /// Default destructor
  ~DynamicDataSocketSet();

  /// Gets the socket sink by supplying its name
  /// @param name the name of the socket to try to get
  /// @throw NoSuchValueException if it does not exist
  template < typename TYPE>
  Common::SafePtr<DataSocketSink<TYPE, COMTYPE> > getSocketSink(const std::string& name) const;

  /// Add a socket sink to the set
  /// @param name the name of the socket to try to add
  /// @throw Common::StorageExistsException if it already exist
  template < typename TYPE >
  Common::SafePtr<DataSocketSink<TYPE, COMTYPE> > createSocketSink(const std::string& name);

  /// Add a socket sink to the set
  /// @param name the name of the socket to try to add
  /// @param isEssential a flag to tell if the sink socket needs to be connected to a source
  /// @throw Common::StorageExistsException if it already exist
  template < typename TYPE >
  Common::SafePtr<DataSocketSink<TYPE, COMTYPE> > createSocketSink(const std::string& name, const bool isEssential);

  /// Remove a socket sink from the set
  /// @param name the name of the socket to try to remove
  /// @throw NoSuchValueException if it does not exist
  void deleteSocketSink(const std::string& name);

  /// Gets the socket source by supplying its name
  /// @param name the name of the socket to try to get
  /// @throw NoSuchValueException if it does not exist
  template < typename TYPE >
  Common::SafePtr<DataSocketSource<TYPE, COMTYPE> > getSocketSource(const std::string& name);

  /// Add a socket source to the set
  /// @param name the name of the socket to try to add
  /// @throw Common::StorageExistsException if it already exist
  template < typename TYPE >
  Common::SafePtr<DataSocketSource<TYPE, COMTYPE> > createSocketSource(const std::string& name);

  /// Check if a certain source socket is present in the storage
  bool sourceSocketExists(const std::string& name);

  /// Check if a certain sink socket is present in the storage
  bool sinkSocketExists(const std::string& name);

  /// Remove a socket source from the set
  /// @param name the name of the socket to try to remove
  /// @throw NoSuchValueException if it does not exist
  void deleteSocketSource(const std::string& name);

  /// Delete all sockets
  void deleteAll();

  /// Gets all source sockets
  /// @return all the source sockets cast into BaseDataSocketSource
  std::vector<Common::SafePtr<BaseDataSocketSource> > getAllSourceSockets();

  /// Gets all sink sockets
  /// @return all the source sockets cast into BaseDataSocketSink
  std::vector<Common::SafePtr<BaseDataSocketSink> > getAllSinkSockets();

protected:

  /// Reset the namespace of all sockets inside
  virtual void setNamespace(const std::string& name)
  {
    if ( m_namespace != name )
    {
      m_namespace = name;
      StorageTypeSource::iterator src_itr = m_sources.begin();
      for (  ; src_itr != m_sources.end(); ++src_itr)
        (*src_itr).second->setParentNamespace(name);

      StorageTypeSink::iterator snk_itr = m_sinks.begin();
      for (  ; snk_itr != m_sinks.end(); ++snk_itr)
        (*snk_itr).second->setParentNamespace(name);
    }
  }

protected:

  /// storage of the DataSocketSink's
  StorageTypeSink m_sinks;
  /// storage of the DataSocketSource's
  StorageTypeSource m_sources;

}; // end of class DynamicDataSocketSet

//////////////////////////////////////////////////////////////////////////////

template < typename COMTYPE>
template < typename TYPE >
Common::SafePtr<DataSocketSink<TYPE, COMTYPE> >
DynamicDataSocketSet<COMTYPE>::getSocketSink(const std::string& name) const
{
  BaseDataSocketSink * ptr = m_sinks.getEntry(name);
  return dynamic_cast<DataSocketSink<TYPE,COMTYPE>*>(ptr);
}

//////////////////////////////////////////////////////////////////////////////

template < typename COMTYPE>
template < typename TYPE >
Common::SafePtr<DataSocketSink<TYPE,COMTYPE> >
DynamicDataSocketSet<COMTYPE>::createSocketSink(const std::string& name)
{
  DataSocketSink<TYPE,COMTYPE> * ptr = new DataSocketSink<TYPE,COMTYPE>(name);
  ptr->setParentNamespace(getNamespace());
  m_sinks.addEntry(name,ptr);
  return ptr;
}

//////////////////////////////////////////////////////////////////////////////

template < typename COMTYPE>
template < typename TYPE >
Common::SafePtr<DataSocketSink<TYPE,COMTYPE> >
DynamicDataSocketSet<COMTYPE>::createSocketSink(const std::string& name, const bool isEssential)
{
  DataSocketSink<TYPE,COMTYPE> * ptr = new DataSocketSink<TYPE,COMTYPE>(name, isEssential);
  ptr->setParentNamespace(getNamespace());
  m_sinks.addEntry(name,ptr);
  return ptr;
}

//////////////////////////////////////////////////////////////////////////////

template < typename COMTYPE>
template < typename TYPE >
Common::SafePtr<DataSocketSource<TYPE,COMTYPE> >
DynamicDataSocketSet<COMTYPE>::getSocketSource(const std::string& name)
{
  BaseDataSocketSource * ptr = m_sources.getEntry(name);
  return dynamic_cast<DataSocketSource<TYPE,COMTYPE>*>(ptr);
}

//////////////////////////////////////////////////////////////////////////////

template < typename COMTYPE>
template < typename TYPE >
Common::SafePtr<DataSocketSource<TYPE,COMTYPE> >
DynamicDataSocketSet<COMTYPE>::createSocketSource(const std::string& name)
{
  DataSocketSource<TYPE,COMTYPE> * ptr = new DataSocketSource<TYPE,COMTYPE>(name);
  ptr->setParentNamespace(getNamespace());
  m_sources.addEntry(name,ptr);
  return ptr;
}

//////////////////////////////////////////////////////////////////////////////

template < typename COMTYPE>
DynamicDataSocketSet<COMTYPE>::DynamicDataSocketSet() :
NamespaceMember()
{
}

//////////////////////////////////////////////////////////////////////////////
template < typename COMTYPE>
DynamicDataSocketSet<COMTYPE>::~DynamicDataSocketSet()
{
  deleteAll();
}

//////////////////////////////////////////////////////////////////////////////
template < typename COMTYPE>
void DynamicDataSocketSet<COMTYPE>::deleteAll()
{
  m_sinks.deleteAllEntries();
  m_sources.deleteAllEntries();
}

//////////////////////////////////////////////////////////////////////////////

template < typename COMTYPE>
void DynamicDataSocketSet<COMTYPE>::deleteSocketSink(const std::string& name)
{
  m_sinks.deleteEntry(name);
}

//////////////////////////////////////////////////////////////////////////////

template < typename COMTYPE>
void DynamicDataSocketSet<COMTYPE>::deleteSocketSource(const std::string& name)
{
  m_sources.deleteEntry(name);
}

//////////////////////////////////////////////////////////////////////////////

template < typename COMTYPE>
std::vector< Common::SafePtr<BaseDataSocketSource> >
DynamicDataSocketSet<COMTYPE>::getAllSourceSockets()
{
  std::vector< Common::SafePtr<BaseDataSocketSource> > result;

  result.reserve(m_sources.size());
  std::transform(m_sources.begin(),
                 m_sources.end(),
                 back_inserter(result),
                 Common::GeneralStorage<BaseDataSocketSource>::extract);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

template < typename COMTYPE>
std::vector< Common::SafePtr<BaseDataSocketSink> >
DynamicDataSocketSet<COMTYPE>::getAllSinkSockets()
{
  std::vector< Common::SafePtr<BaseDataSocketSink> > result;

  result.reserve(m_sinks.size());
  std::transform(m_sinks.begin(),
                 m_sinks.end(),
                 back_inserter(result),
                 Common::GeneralStorage<BaseDataSocketSink>::extract);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

template < typename COMTYPE>
bool DynamicDataSocketSet<COMTYPE>::sourceSocketExists(const std::string& name)
{
  return m_sources.checkEntry(name);
}

//////////////////////////////////////////////////////////////////////////////

template < typename COMTYPE>
bool DynamicDataSocketSet<COMTYPE>::sinkSocketExists(const std::string& name)
{
  return m_sinks.checkEntry(name);
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace COOLFluiD

} // namespace Framework

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_DynamicDataSocketSet_hh
