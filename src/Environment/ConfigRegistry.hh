// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Environment_ConfigRegistry_hh
#define COOLFluiD_Environment_ConfigRegistry_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/NonCopyable.hh"

#include "Common/COOLFluiD.hh"
#include "Common/SafePtr.hh"
#include "Common/NoSuchValueException.hh"

#include "Environment/EnvironmentAPI.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Environment {

    class ModuleRegisterBase;
    class ConfigRegisterBase;

//////////////////////////////////////////////////////////////////////////////

/// This class is a singleton object which serves as registry for all the
/// self-configuration object registers.
/// The only instance of this object is held by the CFEnv.
/// @see Common::ConfigObject
/// @author Tiago Quintino
class Environment_API ConfigRegistry :
  public Common::NonCopyable<ConfigRegistry> {

  friend class ModuleRegisterBase;

public:

  /// Dumps all the ConfigOption's registered
  void dumpOptions();

  /// Register a factory
  /// @param factory pointer to factory
  void regist(ConfigRegisterBase* configRegister);

  /// Remove a registered factory
  /// @throw Common::NoSuchValueException if does not exist
  void unregist(const std::string& configRegisterName);

  /// Get a given factory
  /// @throw Common::NoSuchValueException if does not exist
  Common::SafePtr<ConfigRegisterBase> getRegister(const std::string& configRegisterName);

protected:

  /// @return the map with the ConfigRegister's available to be used.
  static std::map<std::string, ConfigRegisterBase*>& getRegisterMap();

private: // methods

  /// Constructor is private to allow only the friend classes to build it
  ConfigRegistry();

  /// Default destructor is private to allow only the friend classes to destroy it
  ~ConfigRegistry();

}; // end of class ConfigRegistry

//////////////////////////////////////////////////////////////////////////////

  } // namespace Environment

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOFluiD_Environment_ConfigRegistry_hh
