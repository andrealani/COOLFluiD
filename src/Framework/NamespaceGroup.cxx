// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.



#include "NamespaceGroup.hh"
#include "Config/DuplicateNameException.hh"
#include "Common/NoSuchValueException.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

NamespaceGroup::NamespaceGroup(const std::string& name) :
_namespaceList(0)
{
  _namespaceList.push_back(name);
}

//////////////////////////////////////////////////////////////////////////////

NamespaceGroup::~NamespaceGroup()
{
}

//////////////////////////////////////////////////////////////////////////////

void NamespaceGroup::addNamespace(const std::string& name)
{
  if(!match(name)) {
    _namespaceList.push_back(name);
  }
  else {
    std::string msg("Namespace ");
    msg += name;
    msg += " already exists";
    throw Config::DuplicateNameException (FromHere(),msg);
  }
}

//////////////////////////////////////////////////////////////////////////////

void NamespaceGroup::removeNamespace(const std::string& name)
{
  if(match(name)) {
    StorageType::iterator itr = find(_namespaceList.begin(),_namespaceList.end(),name);
    _namespaceList.erase(itr);
  }
  else {
    std::string msg("Namespace ");
    msg += name;
    msg += " does not exist";
    throw Common::NoSuchValueException (FromHere(),msg);
  }
}

//////////////////////////////////////////////////////////////////////////////

bool NamespaceGroup::match(const std::string& name)
{
  StorageType::iterator itr = find(_namespaceList.begin(),_namespaceList.end(),name);
  return (itr != _namespaceList.end());
}

//////////////////////////////////////////////////////////////////////////////

void NamespaceGroup::removeAll()
{
  _namespaceList.clear();
}

//////////////////////////////////////////////////////////////////////////////

std::string NamespaceGroup::getPrimaryNamespace()
{
  cf_assert(!_namespaceList.empty());
  return _namespaceList[0];
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

