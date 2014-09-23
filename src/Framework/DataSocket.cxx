// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Common/CFLog.hh"
#include "Framework/DataSocket.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

void DataSocket::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< std::string >("Namespace","Namespace of this Socket.");
}

//////////////////////////////////////////////////////////////////////////////

DataSocket::DataSocket(const std::string& name, const std::string& storage, const std::string& type) :
NamespaceMember(),
ConfigObject(name),
m_dataSocketStorage(storage),
m_dataSocketType(type)
{
   addConfigOptionsTo(this);
   setParameter("Namespace", &m_namespace);
}

//////////////////////////////////////////////////////////////////////////////

DataSocket::~DataSocket()
{
}

//////////////////////////////////////////////////////////////////////////////

DataSocket::DataSocket (const DataSocket& ss) :
  NamespaceMember( ss ),
  ConfigObject( ss ),
  m_dataSocketStorage ( ss.m_dataSocketStorage ),
  m_dataSocketType ( ss.m_dataSocketType )
{
}

//////////////////////////////////////////////////////////////////////////////

DataSocket& DataSocket::operator=(const DataSocket& ss)
{
  NamespaceMember::operator=(ss);
  ConfigObject::operator=(ss);
  return *this;
}

//////////////////////////////////////////////////////////////////////////////

void DataSocket::configure ( Config::ConfigArgs& args )
{
  ConfigObject::configure(args);

  // configuring the socket namespace to which it belongs
  // this has precedence over the parent socket namespace
  if( m_namespace == NamespaceMember::defaultNamespace() )
  {
    setSelfNamespace (m_namespace);
  }
}

//////////////////////////////////////////////////////////////////////////////

DataBroker::key_t DataSocket::makeID () const
{
  return Common::make_Quartet ( getDataSocketName() , getDataSocketStorage() , getDataSocketType() , getNamespace() );
}

//////////////////////////////////////////////////////////////////////////////

std::string DataSocket::getDataSocketFullStorageName() const
{
  return getNamespace() + "_" + getDataSocketName();
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

