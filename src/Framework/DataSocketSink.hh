// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_DataSocketSink_hh
#define COOLFluiD_Framework_DataSocketSink_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/BaseDataSocketSink.hh"
#include "Framework/DataSocketSource.hh"
#include "Framework/SocketException.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This class is the base class for objects that need a DataSocket of a certain TYPE
/// @author Tiago Quintino
template < typename TYPE , typename STORAGETYPE = LOCAL >
class DataSocketSink : public BaseDataSocketSink {
public:

  /// Constructor
  explicit DataSocketSink(const std::string& name, bool isEssential = true);

  /// Constructor from a DataSocketSource
  DataSocketSink(DataSocketSource<TYPE,STORAGETYPE>& source);

  /// Copy constructor
  DataSocketSink<TYPE,STORAGETYPE>(const DataSocketSink<TYPE,STORAGETYPE>& sink) : BaseDataSocketSink (sink)
  {
    m_source = sink.m_source;
    m_isEssential = sink.m_isEssential;
    DataBroker::getInstance().registerSink ( this );
  }

  /// Operator =
  const DataSocketSink<TYPE,STORAGETYPE>& operator=(const DataSocketSink<TYPE,STORAGETYPE>& sink)
  {
    DataBroker::getInstance().unregisterSink ( this );
    BaseDataSocketSink::operator=(sink);
    m_isEssential = sink.m_isEssential;
    DataBroker::getInstance().registerSink ( this );
    return *this;
  }

  /// Default destructor.
  virtual ~DataSocketSink();

  /// Returns the DataHandle associated with the Socket
  /// @throw NoSuchStorageException if the socket is not connected, which may happen
  ///        if there is no matching source
  /// @return the DataHandle held by the source socket to which is connected
  DataHandle<TYPE,STORAGETYPE> getDataHandle()
  {
    cf_assert(checkConnection());
    return m_source->getDataHandle();
  }

  /// Returns the DataHandle associated with the Socket
  /// @throw NoSuchStorageException if the socket is not connected, which may happen
  ///        if there is no matching source
  /// @return the constant DataHandle held by the source socket to which is connected
  const DataHandle<TYPE,STORAGETYPE> getDataHandle() const
  {
    cf_assert(checkConnection());
    return m_source->getDataHandle();
  }
  
  /// Connects this Sink to the supplied Source
  /// @param source the DataSocket to connect to
  void connectTo (Common::SafePtr<BaseDataSocketSource> source);

  /// Unpluggs this socket
  void unplug () {  m_source.reset(CFNULL); }

  /// Indicates if this Sink is essential, meaning it must satisfied or else a
  /// consistency exception should be raised.
  bool isEssential () const {  return m_isEssential; }

  /// Indicates if this Sink is connected
  bool isConnected () const { return m_source.isNotNull(); }
  
  /// @return the global size of the underlying data array
  virtual CFuint getGlobalSize() const {return (m_isEssential) ? m_source->getDataHandle().getGlobalSize() : 0;}
  
  /// @return the local size of the underlying data array
  virtual CFuint getLocalSize() const {return (m_isEssential) ? m_source->getDataHandle().getLocalSize() : 0;}
  
private:
  
  /// Check if the socket is connected // this function is critical for performance and it is not inlinable
  /// alternatively implement a vaidate function that does some sanity check to all Data Sockets
  bool checkConnection() const 
  {
    if ( !isConnected() )
      throw SocketException ( FromHere(), "Data socket sink [" + this->makeID().str() + "] was is not connected to a matching source\n" );
    return true;
  }
  
private:

  /// the link to the connected DataSocketSource
  Common::SafePtr< DataSocketSource<TYPE,STORAGETYPE> > m_source;

  /// indicates if this Socket is essential and should be satisfied
  bool m_isEssential;

}; // end of class DataSocketSink

//////////////////////////////////////////////////////////////////////////////

template < typename TYPE , typename STORAGETYPE >
void deleteAllPtr(DataSocketSink<TYPE,STORAGETYPE>& socket)
{
}

//////////////////////////////////////////////////////////////////////////////

template < typename TYPE , typename STORAGETYPE >
void deleteAllPtr(DataSocketSink<TYPE*,STORAGETYPE>& socket)
{
  cf_assert ( socket.isConnected() );

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
DataSocketSink<TYPE,STORAGETYPE>::
DataSocketSink(const std::string& name, bool isEssential) :
  BaseDataSocketSink(name,DEMANGLED_TYPEID(STORAGETYPE),DEMANGLED_TYPEID(TYPE)),
  m_source(CFNULL),
  m_isEssential(isEssential)
{
  DataBroker::getInstance().registerSink ( this );
}

//////////////////////////////////////////////////////////////////////////////

template < typename TYPE , typename STORAGETYPE >
DataSocketSink<TYPE,STORAGETYPE>::
DataSocketSink(DataSocketSource<TYPE,STORAGETYPE>& source) :
  BaseDataSocketSink(source.getDataSocketName(),
                     DEMANGLED_TYPEID(STORAGETYPE),
                     DEMANGLED_TYPEID(TYPE)),
  m_source(&source),
  m_isEssential(true)
{
  m_namespace = source.getNamespace();
  DataBroker::getInstance().registerSink ( this );
}

//////////////////////////////////////////////////////////////////////////////

template < typename TYPE , typename STORAGETYPE >
DataSocketSink<TYPE,STORAGETYPE>::~DataSocketSink()
{
  DataBroker::getInstance().unregisterSink ( this );
}

//////////////////////////////////////////////////////////////////////////////

template < typename TYPE , typename STORAGETYPE >
void DataSocketSink<TYPE,STORAGETYPE>::
connectTo(Common::SafePtr<BaseDataSocketSource> source)
{
  unplug();

  CFLogDebugMin (  "CONNECTING socket: " << source->getDataSocketName()
                << " Namespace: "     << source->getNamespace()
                << " Storage: "       << source->getDataSocketStorage()
                << " Type: "          << source->getDataSocketType()
                << "\n");

  cf_assert (  this->makeID () == source->makeID() );

  m_source = source.template d_castTo<DataSocketSource<TYPE,STORAGETYPE> >();
}

//////////////////////////////////////////////////////////////////////////////

} // namespace Framework
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_DataSocketSink_hh
