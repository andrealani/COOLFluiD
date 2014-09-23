// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Environment_FactoryRegistry_hh
#define COOLFluiD_Environment_FactoryRegistry_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/SafePtr.hh"
#include "Common/GeneralStorage.hh"
#include "Common/NonCopyable.hh"

#include "Environment/Environment.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Environment {

    class FactoryBase;
    class CFEnv;

//////////////////////////////////////////////////////////////////////////////

/// This class is a singleton object which serves as registry for all the
/// Factory objects that are created
/// The only instance of this object is held by the CFEnv.
/// @author Tiago Quintino
class Environment_API FactoryRegistry :
  public Common::NonCopyable<FactoryRegistry> {

  friend class Environment::CFEnv;

public:

  /// Register a factory
  /// @param module pointer to a FactoryBase to be added
  void regist(Environment::FactoryBase* factory);

  /// Remove a registered factory
  /// @param name name of a FactoryBase to be removed
  void unregist(const std::string& name);

  /// Get a given factory
  /// @param name name of a FactoryBase to be accessed
  /// @return a pointer to a FactoryBase if found or a null pointer if not found
  Common::SafePtr<Environment::FactoryBase> getFactory(const std::string& name);

private: // methods

  /// Constructor is private to allow only the friend classes to build it
  FactoryRegistry();

  /// Default destructor is private to allow only the friend classes to destroy it
  ~FactoryRegistry();

private: // data

  Common::GeneralStorage<Environment::FactoryBase> m_store;

}; // end of class FactoryRegistry

//////////////////////////////////////////////////////////////////////////////

  } // namespace Environment

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOFluiD_Environment_FactoryRegistry_hh
