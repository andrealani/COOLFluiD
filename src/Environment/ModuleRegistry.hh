// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Environment_ModuleRegistry_hh
#define COOLFluiD_Environment_ModuleRegistry_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/NonCopyable.hh"

#include "Common/SafePtr.hh"
#include "Common/GeneralStorage.hh"

#include "Environment/Environment.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Environment {

    class CFEnv;
    class ModuleRegisterBase;

//////////////////////////////////////////////////////////////////////////////

/// This class is a singleton object which serves as registry for all the
/// modules that are loaded
/// The only instance of this object is held by the CFEnv.
/// @author Tiago Quintino
class Environment_API ModuleRegistry :
  public Common::NonCopyable<ModuleRegistry> {

  friend class Environment::CFEnv;

public:

  /// Register a module
  /// @param module pointer to a ModuleRegisterBase to be added
  void regist(Environment::ModuleRegisterBase* module);

  /// Remove a registered module
  /// @param moduleName name of a ModuleRegisterBase to be removed
  void unregist(const std::string& moduleName);

  /// Checks that a module is registered
  /// @param moduleName name of a ModuleRegisterBase to be checked
  bool isRegistered(const std::string& moduleName);

  /// Get a given module
  /// @param moduleName name of a ModuleRegisterBase to be accessed
  /// @return a pointer to a ModuleRegisterBase if found or a null pointer if not found
  Common::SafePtr<Environment::ModuleRegisterBase> getModuleRegisterBase(const std::string& moduleName);

  /// Get all modules
  /// @return a vector with all the modules
  std::vector< Common::SafePtr<Environment::ModuleRegisterBase> > getAllModules();

private: // methods

  /// Constructor is private to allow only the friend classes to build it
  ModuleRegistry();

  /// Default destructor is private to allow only the friend classes to destroy it
  ~ModuleRegistry();

private: // data

  Common::GeneralStorage<Environment::ModuleRegisterBase> m_store;

}; // end of class ModuleRegistry

//////////////////////////////////////////////////////////////////////////////

  } // namespace Environment

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOFluiD_Environment_ModuleRegistry_hh
