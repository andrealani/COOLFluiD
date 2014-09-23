// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Common/NoSuchValueException.hh"

#include "Framework/GeometricEntityFactory.hh"
#include "Framework/BaseGeometricEntityProvider.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

void GeometricEntityFactory::add(BaseGeometricEntityProvider* provider)
{
  getProviderMap()[provider->getName()] = provider;
}

//////////////////////////////////////////////////////////////////////////////

void GeometricEntityFactory::remove(const std::string& providerName)
{
  if (getProviderMap().count(providerName) == 0) {
    throw Common::NoSuchValueException (FromHere(),"GeometricEntityFactory::remove() "
             + providerName + " not registered");
  }
  getProviderMap().erase(getProviderMap().find(providerName));
}

//////////////////////////////////////////////////////////////////////////////

BaseGeometricEntityProvider* GeometricEntityFactory::getProvider ( const std::string& providerName )
{
  if (getProviderMap().count(providerName) == 0) {
    throw Common::NoSuchValueException (FromHere(),"GeometricEntityFactory::getProvider() "
            + providerName + " not registered");
  }
  return getProviderMap().find(providerName)->second;
}

//////////////////////////////////////////////////////////////////////////////

GeometricEntity* GeometricEntityFactory::create(const std::string& providerName)
{
  if (getProviderMap().count(providerName) == 0) {
    throw Common::NoSuchValueException (FromHere(),"GeometricEntityFactory::create() => shape function type "
				       + providerName + " not registered");
  }
  return (getProviderMap().find(providerName)->second)->create();
}

//////////////////////////////////////////////////////////////////////////////

std::map<std::string, BaseGeometricEntityProvider*>& GeometricEntityFactory::getProviderMap()
{
  /// Storages the provider's that can provide GeometricEntity's
  static std::map<std::string, BaseGeometricEntityProvider*> providers;
  return providers;
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD
