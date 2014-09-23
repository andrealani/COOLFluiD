// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/BaseDataSocketSink.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
namespace Framework {

//////////////////////////////////////////////////////////////////////////////

BaseDataSocketSink::BaseDataSocketSink(const std::string& name, const std::string& storage, const std::string& type) :
DataSocket(name,storage,type)
{
}

//////////////////////////////////////////////////////////////////////////////

BaseDataSocketSink::~BaseDataSocketSink()
{
}

//////////////////////////////////////////////////////////////////////////////

BaseDataSocketSink&
BaseDataSocketSink::operator=(const BaseDataSocketSink& sink)
{
  DataSocket::operator=(sink);
  return *this;
}

//////////////////////////////////////////////////////////////////////////////

BaseDataSocketSink::BaseDataSocketSink(const BaseDataSocketSink& sink) : DataSocket(sink)
{
}

//////////////////////////////////////////////////////////////////////////////

void BaseDataSocketSink::setNamespace (const std::string& name)
{
  if ( m_namespace != name )
  {
    DataBroker::getInstance().unregisterSink ( this );
    m_namespace = name;
    DataBroker::getInstance().registerSink ( this );
  }
}

//////////////////////////////////////////////////////////////////////////////

} // namespace Framework
} // namespace COOLFluiD

