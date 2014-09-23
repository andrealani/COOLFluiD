// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_DataBroker_hh
#define COOLFluiD_Framework_DataBroker_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/NonCopyable.hh"
#include "Common/Quartet.hh"

#include "Framework/Framework.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
namespace Framework {

  class BaseDataSocketSource;
  class BaseDataSocketSink;

//////////////////////////////////////////////////////////////////////////////

/// define a quartet of std::string
typedef Common::Quartet < std::string, std::string, std::string, std::string > StringQuartet;

/// Class to identify a socket
struct Framework_API SocketID : public StringQuartet
{
  SocketID ( const StringQuartet& in ) : StringQuartet(in) {}
  std::string str () const;
};

//////////////////////////////////////////////////////////////////////////////

/// This singleton class connects sockets sources to their sinks
/// @author Tiago Quintino
class Framework_API DataBroker : public Common::NonCopyable<DataBroker> {

public: // typedefs

    typedef SocketID key_t;

    typedef BaseDataSocketSource* source_t;
    typedef BaseDataSocketSink*   sink_t;

    typedef std::map < key_t , source_t >    map_source_t;
    typedef std::multimap < key_t , sink_t > map_sink_t;

    typedef std::pair< map_sink_t::iterator, map_sink_t::iterator > range_t;

public: // functions

  /// @return the instance of this singleton
  static DataBroker& getInstance();

  /// Connect the sink socket to the source socket
  void connectSink ( BaseDataSocketSink* sink );

  /// Register a new source
  /// @throw Common::StorageExistsException in case it already exists
  void registerSource ( BaseDataSocketSource* source );
  /// Unregister a source from the database
  /// @throw Common::NoSuchStorageException in case it does not exist
  void unregisterSource ( BaseDataSocketSource* source );

  /// Register a new sink
  /// @throw Common::StorageExistsException in case it already exists
  void registerSink ( BaseDataSocketSink* sink );
  /// Unregister a sink from the database
  /// @throw Common::NoSuchStorageException in case it does not exist
  void unregisterSink ( BaseDataSocketSink* sink );

private: // functions

  /// Priate constructor
  DataBroker();
  /// Private destructor
  ~DataBroker();

private: // data

  /// register for data sources
  map_source_t  m_regsrcs;
  /// register for data sinks
  map_sink_t  m_regsnks;

}; // end of class DataBroker

//////////////////////////////////////////////////////////////////////////////

} // namespace Framework
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_DataBroker_hh
