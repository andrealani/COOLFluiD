// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Environment/ModuleRegistry.hh"
#include "Environment/ModuleRegisterBase.hh"
#include "Common/CFLog.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Common;

namespace COOLFluiD {

  namespace Environment {

//////////////////////////////////////////////////////////////////////////////

ModuleRegistry::ModuleRegistry() {}

//////////////////////////////////////////////////////////////////////////////

ModuleRegistry::~ModuleRegistry() {}

//////////////////////////////////////////////////////////////////////////////

void ModuleRegistry::regist(Environment::ModuleRegisterBase* module)
{
  if ( !m_store.checkEntry( module->getName()) )
  {
    m_store.addEntry(module->getName(),module);
#ifndef CF_HAVE_CRAYSTATIC
    CFtrace << "Module " + module->getName() + " registered\n";
#endif
  }
  else
  {
#ifndef CF_HAVE_CRAYSTATIC
    CFtrace << "Module " + module->getName() + " already registered : skipping registration\n";
#endif 
 }
}

//////////////////////////////////////////////////////////////////////////////

void ModuleRegistry::unregist(const std::string& moduleName)
{
  if ( m_store.checkEntry( moduleName ) )
  {
    m_store.removeEntry(moduleName);
  }
  else
  {
    CFLogWarn("Module " + moduleName + " not registered : skipping removal\n");
  }
}

//////////////////////////////////////////////////////////////////////////////

bool ModuleRegistry::isRegistered(const std::string& moduleName)
{
  return m_store.checkEntry(moduleName);
}

//////////////////////////////////////////////////////////////////////////////

Common::SafePtr<Environment::ModuleRegisterBase>
ModuleRegistry::getModuleRegisterBase(const std::string& moduleName)
{
  if ( m_store.checkEntry( moduleName ) )
  {
    return m_store.getEntry(moduleName);
  }
  else
  {
    CFLogWarn("Module " + moduleName + " not registered : returning null pointer\n");
    return Common::SafePtr<Environment::ModuleRegisterBase>(CFNULL);
  }
}

//////////////////////////////////////////////////////////////////////////////

std::vector< Common::SafePtr<Environment::ModuleRegisterBase> >
ModuleRegistry::getAllModules()
{
  std::vector< SafePtr<ModuleRegisterBase> > all;
  all.reserve(m_store.size());
  std::transform(m_store.begin(),
                 m_store.end(),
                 back_inserter(all),
                 GeneralStorage<ModuleRegisterBase>::extract);
  return all;
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Environment

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
