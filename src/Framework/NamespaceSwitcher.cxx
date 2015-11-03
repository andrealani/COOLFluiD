// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/NamespaceSwitcher.hh"
#include "Framework/MeshData.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/SubSystemStatus.hh"
#include "Common/PE.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

NamespaceSwitcher& NamespaceSwitcher::getInstance(const std::string& subSystemName)
{
  // somebody will have to delete those entries
  static Common::GeneralStorage<NamespaceSwitcher> nsmap;
  if (nsmap.checkEntry(subSystemName)) {return *nsmap.getEntry(subSystemName);}
  
  // add a new namespace if none corresponds to the given subsystem name
  nsmap.addEntry(subSystemName, new NamespaceSwitcher());
  return *nsmap.getEntry(subSystemName);
}
    
//////////////////////////////////////////////////////////////////////////////
    
NamespaceSwitcher::NamespaceSwitcher() :
m_isEnabled(false)
{
}

//////////////////////////////////////////////////////////////////////////////

NamespaceSwitcher::~NamespaceSwitcher()
{
  deleteAllNamespaces();
}

//////////////////////////////////////////////////////////////////////////////

Common::SafePtr<Namespace>
NamespaceSwitcher::createUniqueNamespace(const std::string& name)
{
  Namespace * ptr;
  if(m_Storage.checkEntry(name)) {
    ptr = m_Storage.getEntry(name);
  }
  else {
    ptr = new Namespace(name);
    m_Storage.addEntry(name,ptr);
  }
  return ptr;
}

//////////////////////////////////////////////////////////////////////////////

void NamespaceSwitcher::deleteAllNamespaces()
{
  // clear the stack
  while (!m_Stack.empty()) {
    m_Stack.pop();
  }

  m_Storage.deleteAllEntries();
}

//////////////////////////////////////////////////////////////////////////////

Common::SafePtr<Namespace> NamespaceSwitcher::getNamespace(const std::string& name)
{
  return m_Storage.getEntry(name);
}

//////////////////////////////////////////////////////////////////////////////

void NamespaceSwitcher::pushNamespace(const std::string& name)
{
  cf_assert(isEnabled());
  Namespace * ptr = m_Storage.getEntry(name);

  m_Stack.push(ptr);

  MeshDataStack::getInstance().push(ptr);
  PhysicalModelStack::getInstance().push(ptr);
  SubSystemStatusStack::getInstance().push(ptr);
}

//////////////////////////////////////////////////////////////////////////////

Common::SafePtr<Namespace> NamespaceSwitcher::getCurrentNamespace()
{
  cf_assert(isEnabled());
  cf_assert(!m_Stack.empty());
  return m_Stack.top();
}

//////////////////////////////////////////////////////////////////////////////

Common::SafePtr<Namespace> NamespaceSwitcher::popNamespace()
{
  cf_assert(isEnabled());
  cf_assert(!m_Stack.empty());
  Namespace * ptr = m_Stack.top();

  Common::SafePtr<MeshData> mdp = MeshDataStack::getInstance().pop();
  cf_assert(mdp->getName() == ptr->getMeshDataName());

  Common::SafePtr<PhysicalModel> pmp = PhysicalModelStack::getInstance().pop();
  cf_assert(pmp->getName() == ptr->getPhysicalModelName());

  Common::SafePtr<SubSystemStatus> ssp =
    SubSystemStatusStack::getInstance().pop();
  cf_assert(ssp->getName() == ptr->getSubSystemStatusName());

  m_Stack.pop();
  return ptr;
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<Namespace> >
NamespaceSwitcher::getAllNamespaces()
{
  std::vector<Common::SafePtr<Namespace> > ret;
  ret.reserve(m_Storage.size());
  std::transform(m_Storage.begin(),
                 m_Storage.end(),
                 back_inserter(ret),
                 Common::GeneralStorage<Namespace>::extract);
  return ret;
}

//////////////////////////////////////////////////////////////////////////////

bool NamespaceSwitcher::isEnabled()
{
  return m_isEnabled;
}

//////////////////////////////////////////////////////////////////////////////

void NamespaceSwitcher::setEnabled(bool isEnabled)
{
  m_isEnabled = isEnabled;
}

//////////////////////////////////////////////////////////////////////////////

std::string NamespaceSwitcher::getName
(std::const_mem_fun_t<std::string, Namespace> fun, const bool filterCoupling) 
{
  typedef std::vector<Common::SafePtr<Namespace> > NspVec;
  NspVec lst = getAllNamespaces();
  
  cf_assert(!lst.empty());
  const int rank = Common::PE::GetPE().GetRank("Default");
  for(NspVec::iterator nsp = lst.begin(); nsp != lst.end(); ++nsp) {
    const std::string nspaceName = (*nsp)->getName();
    if (Common::PE::GetPE().isRankInGroup(rank, nspaceName)) {
      if (!filterCoupling) {
	return fun(&*(*nsp)); // Namespace* from SafePtr<Namespace>
      }
      else {
	cf_assert(filterCoupling);
	if (!(*nsp)->isForCoupling()) {
	  return fun(&*(*nsp)); // Namespace* from SafePtr<Namespace>
	}
      }
    }
  }
  cf_assert(false);
  return "";
}
    
//////////////////////////////////////////////////////////////////////////////
 
CFuint NamespaceSwitcher::getID(const bool filterCoupling)
{
  typedef std::vector<Common::SafePtr<Namespace> > NspVec;
  NspVec lst = getAllNamespaces();
  
  cf_assert(!lst.empty());
  const int rank = Common::PE::GetPE().GetRank("Default");
  
  CFuint idx = 0;
  for(NspVec::iterator nsp = lst.begin(); nsp != lst.end(); ++nsp, ++idx) {
    const std::string nspaceName = (*nsp)->getName();
    if (Common::PE::GetPE().isRankInGroup(rank, nspaceName)) {
      if (!filterCoupling) {
	return idx;
      }
      else {
	cf_assert(filterCoupling);
	if (!(*nsp)->isForCoupling()) {
	  return idx; // Namespace* from SafePtr<Namespace>
	}
      }
    }
  }
  
  return 0;
}

//////////////////////////////////////////////////////////////////////////////
 
  } // namespace COOLFluiD

} // namespace Framework

//////////////////////////////////////////////////////////////////////////////

