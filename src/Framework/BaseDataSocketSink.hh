// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_BaseDataSocketSink_hh
#define COOLFluiD_Framework_BaseDataSocketSink_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/DataSocket.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

    class BaseDataSocketSource;

//////////////////////////////////////////////////////////////////////////////

/// This class is the base class for objects that provide a DataSocket
/// @author Tiago Quintino
class Framework_API BaseDataSocketSink : public DataSocket {
public:

  /// Default constructor without arguments.
  BaseDataSocketSink(const std::string& name, const std::string& storage, const std::string& type);

  /// Default destructor.
  virtual ~BaseDataSocketSink();

  /// Connects this Sink to the supplied Source
  /// @param source the DataSocket to connect to
  virtual void connectTo(Common::SafePtr<BaseDataSocketSource> source) = 0;

  /// Unplugs this sockets.
  virtual void unplug() = 0;

  /// Indicates if this Sink is essential, meaning it must satisfied or else a
  /// consistency exception should be raised.
  virtual bool isEssential() const = 0;

  /// Indicates if this Sink is connected
  virtual bool isConnected() const = 0;

  /// Reset the namespace and reregister in the DataBroker
  virtual void setNamespace(const std::string& name);

protected:

  /// copy construction is not allowed
  BaseDataSocketSink(const BaseDataSocketSink& sink);
  /// Operator =
  BaseDataSocketSink& operator= (const BaseDataSocketSink& sink);

}; // end of class BaseDataSocketSink

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_BaseDataSocketSink_hh
