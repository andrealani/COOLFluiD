// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Environment_SelfRegistry_hh
#define COOLFluiD_Environment_SelfRegistry_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/SafePtr.hh"
#include "Common/GeneralStorage.hh"
#include "Common/NonCopyable.hh"

#include "Environment/EnvironmentAPI.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Common { class ProviderBase; }

  namespace Environment {

    class ModuleRegisterBase;

//////////////////////////////////////////////////////////////////////////////

/// This class is a singleton object which serves as registry for all the
/// self-registration object providers.
/// The only instance of this object is held by the ModuleRegister
/// @see Common::ProviderBase
/// @see Common::Provider
/// @author Tiago Quintino
class Environment_API SelfRegistry :
  public Common::NonCopyable<SelfRegistry> {

  friend class Environment::ModuleRegisterBase;

public:

  /// Register a Object provider
  /// @param factory pointer to the provider
  void regist(Common::ProviderBase* provider);

  /// Remove a registered provider
  /// @param name of the provider
  void unregist(const std::string& name, const std::string& type);

  /// Remove a registered provider
  /// @param name of the provider
  void unregist(Common::ProviderBase* provider);

  /// Get a given provider by his name
  /// @param name of the provider
  /// @return a SafePtr to the provider
  Common::SafePtr<Common::ProviderBase> getProvider(const std::string& name, const std::string& type);

private: // methods

  /// Constructor is private to allow only the friend classes to build it
  SelfRegistry();

  /// Default destructor is private to allow only the friend classes to destroy it
  ~SelfRegistry();

private: // data

  std::map<std::string,Common::GeneralStorage<Common::ProviderBase> > m_store;

}; // end of class SelfRegistry

//////////////////////////////////////////////////////////////////////////////

  } // namespace Environment

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOFluiD_Environment_SelfRegistry_hh
