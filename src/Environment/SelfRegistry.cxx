// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Common/ProviderBase.hh"

#include "Environment/SelfRegistry.hh"
#include "Common/CFLog.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Common;

namespace COOLFluiD {
  namespace Environment {

//////////////////////////////////////////////////////////////////////////////

SelfRegistry::SelfRegistry()
{
}

//////////////////////////////////////////////////////////////////////////////

SelfRegistry::~SelfRegistry()
{
}

//////////////////////////////////////////////////////////////////////////////

void SelfRegistry::regist(Common::ProviderBase* provider)
{
  const std::string& name = provider->getProviderName();
  const std::string& type = provider->getProviderType();
  if ( !m_store[type].checkEntry(name) )
  {
    m_store[type].addEntry(name,provider);
  }
  else
  {
#ifdef CF_HAVE_LOG4CPP
   CFtrace << "Provider [" + provider->getProviderName()
              + "] of type [" + provider->getProviderType()
              + "] already registered : skipping registration\n";
#endif  
  }
}

//////////////////////////////////////////////////////////////////////////////

void SelfRegistry::unregist(const std::string& name, const std::string& type)
{
  if ( m_store[type].checkEntry(name) )
  {
    m_store[type].removeEntry(name);
  }
  else
  {
#ifdef CF_HAVE_LOG4CPP
    CFtrace << "Provider [" + name
              + "] of type [" + type
              + "] not registered : skipping removal\n";
#endif 
 }
}

//////////////////////////////////////////////////////////////////////////////

void SelfRegistry::unregist(Common::ProviderBase* provider)
{

  const std::string& name = provider->getProviderName();
  const std::string& type = provider->getProviderType();
  if ( m_store[type].checkEntry(name) )
  {
    unregist(name,type);
  }
  else
  {
#ifdef CF_HAVE_LOG4CPP
    CFtrace << "Provider ["  + provider->getProviderName()
              + "] of type [" + provider->getProviderType()
              + "] not registered : skipping removal\n";
#endif
  }
}

//////////////////////////////////////////////////////////////////////////////

Common::SafePtr<Common::ProviderBase>
SelfRegistry::getProvider(const std::string& name, const std::string& type)
{
  if ( m_store[type].checkEntry(name) )
  {
    return m_store[type].getEntry(name);
  }
  else
  {
#ifdef CF_HAVE_LOG4CPP
    CFLog(VERBOSE, "Provider [" + name
              + "] of type [" + type
              + "] not registered : returning null pointer\n");
#endif
    return Common::SafePtr<Common::ProviderBase>(CFNULL);
 }
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Environment
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
