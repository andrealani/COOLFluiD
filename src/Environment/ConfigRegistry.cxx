// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Environment/ConfigRegistry.hh"
#include "Environment/ConfigRegisterBase.hh"
#include "Common/NotImplementedException.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Environment {

//////////////////////////////////////////////////////////////////////////////

ConfigRegistry::ConfigRegistry()
{
}

//////////////////////////////////////////////////////////////////////////////

ConfigRegistry::~ConfigRegistry()
{
}

//////////////////////////////////////////////////////////////////////////////

void ConfigRegistry::dumpOptions()
{
  throw Common::NotImplementedException (FromHere(),"ConfigRegistry::dumpOptions()");
}

//////////////////////////////////////////////////////////////////////////////

void ConfigRegistry::regist(ConfigRegisterBase* configRegister)
{
  getRegisterMap()[configRegister->getTypeName()] = configRegister;
}

//////////////////////////////////////////////////////////////////////////////

void ConfigRegistry::unregist(const std::string& configRegisterName)
{
  if (getRegisterMap().count(configRegisterName) == 0)
  {
    throw Common::NoSuchValueException (FromHere(),"ConfigRegistry::unregist() " + configRegisterName + " not registered");
  }
  getRegisterMap().erase(getRegisterMap().find(configRegisterName));
}

//////////////////////////////////////////////////////////////////////////////

Common::SafePtr<ConfigRegisterBase>
ConfigRegistry::getRegister(const std::string& configRegisterName)
{
  if (getRegisterMap().count(configRegisterName) == 0)
  {
    throw Common::NoSuchValueException (FromHere(),"ConfigRegistry::getFactory() " + configRegisterName + " not registered");
  }

  return getRegisterMap().find(configRegisterName)->second;
}

//////////////////////////////////////////////////////////////////////////////

std::map<std::string, ConfigRegisterBase*>& ConfigRegistry::getRegisterMap()
{
  static std::map<std::string, ConfigRegisterBase*> configRegistry;
  return configRegistry;
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Environment

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
