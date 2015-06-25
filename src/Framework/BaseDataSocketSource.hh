// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_BaseDataSocketSource_hh
#define COOLFluiD_Framework_BaseDataSocketSource_hh

//////////////////////////////////////////////////////////////////////////////

#include "DataSocket.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
namespace Framework {

    class DataStorage;

//////////////////////////////////////////////////////////////////////////////

/// This class is the base class for objects that provide a DataSocket
/// @author Tiago Quintino
class Framework_API BaseDataSocketSource :
  public DataSocket,
  public Common::NonCopyable<BaseDataSocketSource>
{

public:

  /// Default constructor without arguments.
  BaseDataSocketSource(const std::string& name,
  	       const std::string& storage,
  	       const std::string& type);

  /// Default destructor.
  virtual ~BaseDataSocketSource();

  /// Allocation of the DataHandle in this DataSocket
  virtual void allocate(Common::SafePtr<DataStorage> storage, const std::string& nspaceName) = 0;
  
  /// Deallocation of the DataHandle in this DataSocket
  virtual void deallocate() = 0;
  
  /// Cheks if this socket has been allocated
  virtual bool isAllocated() const = 0;
  
  /// Reset the namespace and reregister in the DataBroker
  virtual void setNamespace(const std::string& name);

}; // end of class BaseDataSocketSource

//////////////////////////////////////////////////////////////////////////////

} // namespace Framework
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_BaseDataSocketSource_hh
