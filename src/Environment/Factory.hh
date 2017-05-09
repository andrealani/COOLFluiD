// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Common_Factory_hh
#define COOLFluiD_Common_Factory_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/FactoryBase.hh"
#include "Common/CFLog.hh"
#include "Common/SafePtr.hh"
#include "Common/NoSuchValueException.hh"
#include "Common/StorageExistsException.hh"
#include "Environment/FactoryRegistry.hh"
#include "Environment/Provider.hh"
#include "Environment/CFEnv.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Environment {

    template < class BASE > class Provider;

//////////////////////////////////////////////////////////////////////////////

/// Stores the Provider's that create self-registration objects polymorphically.
/// It is part of an abstract factory pattern implementation.
/// @author Andrea Lani
/// @author Tiago Quintino
template <class BASE>
class Factory : public FactoryBase {
public: // methods

  /// @return the instance of this singleton
  static Environment::Factory<BASE>& getInstance();

  /// Checks if a provider is registered
  /// @param name name of the provider
  bool exists(const std::string& name);

  /// Registers a provider
  /// @param provider pointer to the provider to be registered
  void regist(Environment::Provider<BASE>* provider);

  /// Remove a registered provider
  /// @throw Common::NoSuchValueException if the provider is not registered
  void unregist(const std::string& providerName);

  /// @return the name of the BASE of this factory
  std::string getTypeName() const { return BASE::getClassName(); }

  /// @return all the providers in this Factory
  std::vector<Common::ProviderBase*> getAllProviders();

  /// Get a given Provider
  /// @throw Common::NoSuchValueException if the Provider is not registered
  Common::SafePtr< typename BASE::PROVIDER > getProvider(const std::string& providerName);

protected: // helper function

  /// Constructor is protected because this is a Singleton.
  Factory();
  /// Destructor is protected because this is a Singleton.
  ~Factory();

  /// providers database
  std::map<std::string, Environment::Provider<BASE>*>& getProviderMap() { return m_map; }

private: // data

  /// providers database
  std::map<std::string, Provider<BASE>*> m_map;

}; // end of class Factory

//////////////////////////////////////////////////////////////////////////////

template <class BASE>
Factory<BASE>& Factory<BASE>::getInstance()
{
  static Environment::Factory<BASE> obj;
  return obj;
}

//////////////////////////////////////////////////////////////////////////////

template <class BASE>
Factory<BASE>::Factory()
{
  Environment::CFEnv::getInstance().getFactoryRegistry()->regist(this);
}

//////////////////////////////////////////////////////////////////////////////

template <class BASE>
Factory<BASE>::~Factory()
{
}

//////////////////////////////////////////////////////////////////////////////

template <class BASE>
void Factory<BASE>::regist(Provider<BASE>* provider)
{
  if (exists(provider->getName()))
  {
#ifndef CF_HAVE_ALLSTATIC
    throw Common::StorageExistsException (FromHere(),
      "In factory of [" + BASE::getClassName() +
      "] a provider with the name [" + provider->getName() +
      "] was found when trying to regist it" );
#endif
  }
  CFtrace << "Registering provider [" << provider->getName() << "] in factory of [" << BASE::getClassName() << "]\n";
  getProviderMap()[provider->getName()] = provider;
}

//////////////////////////////////////////////////////////////////////////////

template <class BASE>
bool Factory<BASE>::exists(const std::string& name)
{
  return (getProviderMap().count(name) > 0);
}

//////////////////////////////////////////////////////////////////////////////

template <class BASE>
void Factory<BASE>::unregist(const std::string& providerName)
{
  if (!exists(providerName))
  {
    throw Common::NoSuchValueException (FromHere(),
      "In factory of [" + BASE::getClassName() +
      "] a provider with the name [" + providerName +
      "] was not found while trying to unregist it" );
  }
  getProviderMap().erase(getProviderMap().find(providerName));
}

//////////////////////////////////////////////////////////////////////////////

template <class BASE>
Common::SafePtr< typename BASE::PROVIDER >
Factory<BASE>::getProvider(const std::string& providerName)
{
  if (!exists(providerName))
  {
    throw Common::NoSuchValueException (FromHere(),
      "In factory of [" + BASE::getClassName() +
      "] a provider with the name [" + providerName +
      "] was not found while trying to get the provider" );
  }

  return dynamic_cast<typename BASE::PROVIDER*>
    (getProviderMap().find(providerName)->second);
}

//////////////////////////////////////////////////////////////////////////////

template <class BASE>
std::vector<Common::ProviderBase*> Factory<BASE>::getAllProviders()
{
  std::vector<Common::ProviderBase*> result;
  typename std::map<std::string,Provider<BASE>*>::iterator itr;
  itr = getProviderMap().begin();
  for (; itr != getProviderMap().end(); ++itr) {
    result.push_back(itr->second);
  }
  return result;
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Environment
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Common_Factory_hh
