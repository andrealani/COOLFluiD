// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include <iostream>

#include "Environment/ModuleRegisterBase.hh"
#include "Environment/CFEnv.hh"
#include "Environment/ModuleRegistry.hh"
#include "Common/CFLog.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace Environment {

//////////////////////////////////////////////////////////////////////////////

ModuleRegisterBase::ModuleRegisterBase(const std::string& name) :
 NamedObject(name),
 m_selfRegistry(),
 m_configRegistry(),
 m_init(false)
{
#ifdef CF_HAVE_LOG4CPP
  CFtrace << "ModuleRegisterBase::ModuleRegisterBase() => Registering module [" << name << "]\n";
#endif
  Environment::CFEnv::getInstance().getModuleRegistry()->regist(this);
}

//////////////////////////////////////////////////////////////////////////////

ModuleRegisterBase::~ModuleRegisterBase()
{
}

//////////////////////////////////////////////////////////////////////////////

Environment::SelfRegistry& ModuleRegisterBase::getSelfRegistry()
{
  return m_selfRegistry;
}

//////////////////////////////////////////////////////////////////////////////

Environment::ConfigRegistry& ModuleRegisterBase::getConfigRegistry()
{
  return m_configRegistry;
}

//////////////////////////////////////////////////////////////////////////////

void ModuleRegisterBase::initiate()
{
  if (!isInitialized())
  {
    // does nothing by default
    m_init = true;
  }
}

//////////////////////////////////////////////////////////////////////////////

void ModuleRegisterBase::terminate()
{
  if(isInitialized())
  {
    // does nothing by default
    m_init = false;
  }
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Environment

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
