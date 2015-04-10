// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/QualifiedName.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

QualifiedName::QualifiedName(const std::string& namesp, const std::string& name) :
 NamedObject(name),
 m_namespace(namesp)
{
  m_qname = m_namespace + QualifiedName::separator() + getName();
}

//////////////////////////////////////////////////////////////////////////////

QualifiedName::~QualifiedName()
{
}

//////////////////////////////////////////////////////////////////////////////

bool QualifiedName::operator==(const QualifiedName& qnr) const
{
  return (m_namespace == qnr.m_namespace) && (m_qname == qnr.m_qname);
}

//////////////////////////////////////////////////////////////////////////////

CFchar QualifiedName::separator ()
{
  return ':';
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

