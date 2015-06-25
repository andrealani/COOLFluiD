// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/NamespaceMember.hh"
#include "Framework/NamespaceSwitcher.hh"
#include "Framework/SubSystemStatus.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

NamespaceMember::NamespaceMember(const std::string& name) :
m_namespace(name),
m_sameAsParent(true)
{
}

//////////////////////////////////////////////////////////////////////////////

NamespaceMember::~NamespaceMember() {}

//////////////////////////////////////////////////////////////////////////////

NamespaceMember& NamespaceMember::operator=(const NamespaceMember& other)
{
  m_namespace    =  other.m_namespace;
  m_sameAsParent = false;
  return *this;
}

//////////////////////////////////////////////////////////////////////////////

void NamespaceMember::setSelfNamespace(const std::string& name)
{
  m_sameAsParent = false;
  setNamespace( name );
}

//////////////////////////////////////////////////////////////////////////////

void NamespaceMember::setParentNamespace(const std::string& name)
{
  if ( m_sameAsParent )
    setNamespace( name );
}

//////////////////////////////////////////////////////////////////////////////

void NamespaceMember::setNamespace(const std::string& name)
{
    m_namespace = name;
}

//////////////////////////////////////////////////////////////////////////////

Common::SafePtr<Namespace> NamespaceMember::getNamespacePtr() const
{
  return NamespaceSwitcher::getInstance
    (SubSystemStatusStack::getCurrentName()).getNamespace(m_namespace);
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace COOLFluiD

} // namespace Framework

//////////////////////////////////////////////////////////////////////////////

