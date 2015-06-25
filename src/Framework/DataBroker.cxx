// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Common/OSystem.hh"
#include "Common/ProcessInfo.hh"

#include "Framework/DataBroker.hh"
#include "Framework/BaseDataSocketSink.hh"
#include "Framework/BaseDataSocketSource.hh"
#include "Framework/MeshData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
namespace Framework {

//////////////////////////////////////////////////////////////////////////////

template < typename SOCKET_T >
std::string dump_socket ( SOCKET_T* p )
{
  ostringstream res;
  res << p << "," << p->makeID().str();
  return res.str();
}

//////////////////////////////////////////////////////////////////////////////

std::string SocketID::str () const
{
  ostringstream res;
  res << first  << ","
      << second << ","
      << third  << ","
      << fourth;
  return res.str();
}

//////////////////////////////////////////////////////////////////////////////

std::string dump_contents ( DataBroker::map_source_t srcs, DataBroker::map_sink_t snks  )
{
  ostringstream res;
  res << "\n ------------------------------------------------------ \n";
  res << " ---- Total sources ---- \n";
  for ( DataBroker::map_source_t::iterator source = srcs.begin(); source != srcs.end(); ++source )
  {
      res << "\t\t source [" << dump_socket(source->second) << "]";
      source->second->isAllocated() ? res << " [allocated]\n" : res << " [not allocated]\n";
  }
  res << " ---- Total sinks ---- \n";
  for ( DataBroker::map_sink_t::iterator sink = snks.begin(); sink != snks.end(); ++sink )
  {
      res << "\t\t sink [" << dump_socket(sink->second) << "]";
      sink->second->isConnected() ? res << " [connected]\n" : res << " [not connected]\n";
  }
  res << " ------------------------------------------------------ \n\n";
  return res.str();
}

//////////////////////////////////////////////////////////////////////////////

DataBroker::DataBroker() : m_regsrcs (), m_regsnks ()
{
}

//////////////////////////////////////////////////////////////////////////////

DataBroker::~DataBroker()
{
}

//////////////////////////////////////////////////////////////////////////////

DataBroker& DataBroker::getInstance()
{
   static DataBroker singleton;
   return singleton;
}

//////////////////////////////////////////////////////////////////////////////

void DataBroker::connectSink ( BaseDataSocketSink* sink )
{
// std::cout << " +++ \n";

  const key_t lkey = sink->makeID();
  std::map < key_t, source_t >::iterator source = m_regsrcs.find ( lkey );
  if ( source == m_regsrcs.end() )
  {
    ostringstream msg;
    msg << "DataBroker: Trying to connect sink ["<< lkey.str() <<"] but matching source is not registered\n";
    if ( sink->isEssential() )
      throw Common::NoSuchStorageException ( FromHere(), msg.str() );
    else
    {
      CFLogDebugMin ( msg.str() );
      return;
    }
  }

  cf_assert ( lkey == source->second->makeID() );

  sink->connectTo( source->second );
}

//////////////////////////////////////////////////////////////////////////////

void DataBroker::registerSource ( BaseDataSocketSource* source )
{

  // std::cout << " +++ " << std::endl;
// std::cout << " +++ >>> trying to register source [" << dump_socket (source) << "]" << std::endl;
  
// std::cout << "\n" << OSystem::getInstance().getProcessInfo()->getBackTrace() << "\n" << std::endl;

  const key_t lkey = source->makeID();
  std::map < key_t, source_t >::const_iterator fitr = m_regsrcs.find ( lkey );
  bool is_namespace_not_null = !( source->getNamespace() == DataSocket::defaultNamespace() );
  if ( fitr != m_regsrcs.end() && is_namespace_not_null )
  {
    ostringstream msg;
    msg << "DataBroker: Trying to register source [" << dump_socket (source) << "] but it already registered\n";
    throw Common::StorageExistsException ( FromHere(), msg.str() );
  }

  m_regsrcs [ lkey ] = source;
// std::cout << " +++ >>> registering source [" << dump_socket (source) << "]" << std::endl;

  cf_assert (  m_regsrcs.find ( lkey ) != m_regsrcs.end() );

  if ( is_namespace_not_null )
  {
    // get all the matching sinks and connect to this source
    range_t range = m_regsnks.equal_range ( lkey );
    for ( map_sink_t::iterator sink = range.first; sink != range.second; ++sink )
    {
      if ( ! sink->second->isConnected() )
      {
        sink->second->connectTo( source );
  // std::cout << " +++ >>> connecting sink [" << dump_socket (sink->second) << "]\n";
      }
      else
      {
        ostringstream msg;
        msg << "DataBroker: found sink [" << dump_socket (sink->second)<< "] that was already connected to unknown source\n";
        Common::ShouldNotBeHereException ( FromHere() , msg.str() );
      }
    }
  }

  // loop on the MeshData's to allocate the socket in the MeshData with matching namespace
  typedef std::vector<Common::SafePtr<MeshData> > MDLst;
  MDLst mds = MeshDataStack::getInstance().getAllEntries();
  for( MDLst::iterator meshData = mds.begin(); meshData != mds.end(); ++meshData )
  {
      if ( (*meshData)->match( source->getNamespace() ) )
      {
        CFLog(VERBOSE, "DataBroker: allocating source [" << dump_socket (source) << "]\n" );
// std::cout << " +++ >>> allocating source [" << dump_socket (source) << "]\n";

        source->allocate( (*meshData)->getDataStorage(), source->getNamespace() );
      }
  }

// std::cout << dump_contents ( m_regsrcs, m_regsnks );
}

//////////////////////////////////////////////////////////////////////////////

void DataBroker::unregisterSource ( BaseDataSocketSource* source )
{
  // std::cout << " +++ " << std::endl;
// std::cout << " +++ >>> trying to unregister source [" << dump_socket (source) << "]" << std::endl;

  const key_t lkey = source->makeID();
  std::map < key_t, source_t >::iterator fitr = m_regsrcs.find ( lkey );
  if ( fitr == m_regsrcs.end() )
  {
    ostringstream msg;
    msg << "DataBroker: Trying to unregist source [" << dump_socket (source) << "] but it is not registered\n";
    throw Common::NoSuchStorageException ( FromHere(), msg.str() );
  }

  m_regsrcs.erase ( fitr );

// std::cout << " +++ >>> unregistering source ["<< dump_socket (source) <<"]\n";

  cf_assert (  m_regsrcs.find ( lkey ) == m_regsrcs.end() );

  // get all the matching sinks and unplug them
  range_t range = m_regsnks.equal_range ( lkey );
  for ( map_sink_t::iterator sink = range.first; sink != range.second; ++sink )
  {
    sink->second->unplug();
  }

//   typedef std::vector<Common::SafePtr<MeshData> > MDLst;
//   MDLst mds = MeshDataStack::getInstance().getAllEntries();
//   for( MDLst::iterator meshData = mds.begin(); meshData != mds.end(); ++meshData )
//   {
//     std::cerr << (*meshData)->getDataStorage()->dump() << std::endl;
//   }
//   std::cerr.flush();

  source->deallocate();

// std::cout << dump_contents ( m_regsrcs, m_regsnks );
}

//////////////////////////////////////////////////////////////////////////////

void DataBroker::registerSink ( BaseDataSocketSink* sink )
{
// std::cout << " +++ \n";
  const key_t lkey = sink->makeID();

  range_t range = m_regsnks.equal_range ( lkey );

  bool found = false;
  for ( map_sink_t::iterator itr = range.first; itr != range.second; ++itr )
  {
/* std::cout << "         >>> regist sink [" << sink << "] loop sink [" << itr->second << ":" << itr->second->makeID().str() <<"]\n"; */
    if ( itr->second == sink ) { found = true ; break; }
  }

  if ( !found ) // not registered so register
  {
    m_regsnks.insert ( std::pair< const key_t, sink_t >( lkey , sink ) );
// std::cout << " +++ >>> registering sink ["<< dump_socket (sink) <<"]\n";
  }
  else
  {
    ostringstream msg;
    msg << "DataBroker: Trying to regist sink ["<< dump_socket (sink) <<"] but it is already registered\n";
    throw Common::StorageExistsException ( FromHere(), msg.str() );
  }

  // find source and plug it
  std::map < key_t, source_t >::iterator source = m_regsrcs.find ( lkey );
  if ( source != m_regsrcs.end() )
  {
    bool is_namespace_not_null = !( source->second->getNamespace() == DataSocket::defaultNamespace() );
    if ( is_namespace_not_null )
      sink->connectTo(source->second);
// std::cout << " +++ >>> connecting sink ["<< dump_socket (sink) <<"] to source [" << dump_socket(source->second) << "]\n";
  }

// std::cout << dump_contents ( m_regsrcs, m_regsnks );
}

//////////////////////////////////////////////////////////////////////////////

void DataBroker::unregisterSink ( BaseDataSocketSink* sink )
{
// std::cout << " +++ \n";
  const key_t lkey = sink->makeID();

  range_t range = m_regsnks.equal_range ( lkey );

  bool found = false;
  for ( map_sink_t::iterator itr = range.first; itr != range.second; ++itr )
  {
/* std::cout << "         >>> loop sink [" << dump_socket(itr->second) <<"]\n"; */
    if ( itr->second == sink )
    {
// std::cout << " +++ >>> unregistering sink ["<< dump_socket(sink) <<"]\n";
        m_regsnks.erase ( itr );
        found = true;
        break;
    }
  }

  if ( !found ) // not registered so register
  {
    ostringstream msg;
    msg << "DataBroker: Trying to unregist sink ["<< dump_socket(sink) <<"] but it is not registered\n";
    throw Common::NoSuchStorageException ( FromHere(), msg.str() );
  }

  if ( sink->isConnected() ) sink->unplug();

// std::cout << dump_contents ( m_regsrcs, m_regsnks );
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

