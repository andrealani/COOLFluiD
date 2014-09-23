// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/BaseDataSocketSource.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

BaseDataSocketSource::BaseDataSocketSource(const std::string& name, const std::string& storage, const std::string& type) :
DataSocket(name,storage,type)
{
}

//////////////////////////////////////////////////////////////////////////////

BaseDataSocketSource::~BaseDataSocketSource()
{
}

//////////////////////////////////////////////////////////////////////////////

void BaseDataSocketSource::setNamespace (const std::string& name)
{
  if ( m_namespace != name )
  {
    DataBroker::getInstance().unregisterSource ( this );
    m_namespace = name;
    DataBroker::getInstance().registerSource ( this );
  }
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

