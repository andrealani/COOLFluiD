// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_DataSocket_hh
#define COOLFluiD_Framework_DataSocket_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/SafePtr.hh"

#include "Config/ConfigObject.hh"

#include "Framework/NamespaceMember.hh"
#include "Framework/Storage.hh"
#include "Framework/DataBroker.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This class is the base class for DataSocket's.
/// @author Tiago Quintino
class Framework_API DataSocket :
    public NamespaceMember,
    public Config::ConfigObject {

public: // functions

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Default constructor without arguments.
  DataSocket(const std::string& name, const std::string& storage, const std::string& type);

  /// Default destructor.
  virtual ~DataSocket();

  /// Configures this Socket
  /// @param args the arguments used for the configuration
  virtual void configure ( Config::ConfigArgs& args );

  /// Accessor to the DataSocket name
  std::string getDataSocketName() const {  return getName(); }

  /// Accessor to the DataSocket storage
  std::string getDataSocketStorage() const {  return m_dataSocketStorage; }

  /// Accessor to the DataSocket type
  std::string getDataSocketType() const {  return m_dataSocketType; }

  /// Set the DataSocket Namespace
  void setDataSocketNamespace(const std::string& name) {  setSelfNamespace(name); }

  /// Returns the name used to store the socket data in the DataStorage
  std::string getDataSocketFullStorageName() const;

  /// Checks if the supplyed DataSocket matchs this one
  DataBroker::key_t makeID () const;
  
   /// @return the global size of the underlying data array
  virtual CFuint getGlobalSize() const = 0;
  
  /// @return the local size of the underlying data array
  virtual CFuint getLocalSize() const = 0;
  
protected:

  /// Copy construction is not allowed
  DataSocket(const DataSocket& ss);
  /// Operator =
  DataSocket& operator= (const DataSocket& ss);

protected: // data

  /// the DataSocket storage
  std::string m_dataSocketStorage;
  /// the DataSocket type
  std::string m_dataSocketType;

}; // end of class DataSocket

//////////////////////////////////////////////////////////////////////////////

} // namespace Framework
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_DataSocket_hh
