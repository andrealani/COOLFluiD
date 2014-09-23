// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/MethodRegistry.hh"
#include "Framework/Method.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

MethodRegistry::~MethodRegistry()
{
  m_storebytype.deleteAllEntries();
}

//////////////////////////////////////////////////////////////////////////////

MethodRegistry& MethodRegistry::getInstance()
{
  static MethodRegistry aInstance;
  return aInstance;
}

//////////////////////////////////////////////////////////////////////////////

Common::SafePtr<Method>
MethodRegistry::getMethod(const std::string& name)
{
  CFLog(VERBOSE, "Warning: you are using a deprecated function. File " << __FILE__ << " Line " << __LINE__ << " Function " << __FUNCTION__ << "\n");
  for ( MtdStorage::iterator itr = m_storebyqname.begin();
        itr != m_storebyqname.end(); ++itr)
        {
          Method * ptr = MtdStorage::extract(*itr);
          if (ptr->getName() == name)
            return ptr;
        }
  throw Common::NoSuchStorageException (FromHere(),"No method with  name " + name + " exists.");
}

//////////////////////////////////////////////////////////////////////////////

Common::SafePtr<Method>
MethodRegistry::getMethod(const QualifiedName& qname)
{
  return m_storebyqname.getEntry(qname.str());
}

//////////////////////////////////////////////////////////////////////////////

Common::SafePtr<Method>
MethodRegistry::getMethod(const std::string& type, const QualifiedName& qname)
{
  insure_mstorage(type);

  QNameStorage * qnst = m_storebytype.getEntry(type);
  QNameStorage::iterator fitr = qnst->find(qname);
  if (fitr != qnst->end())
  {
    cf_assert( *fitr == qname);
    return m_storebyqname.getEntry(qname.str());
  }
  else
    throw Common::NoSuchStorageException (FromHere(),"No method with qualified name " + qname.str() + " and type " + type + " exists.");
}

//////////////////////////////////////////////////////////////////////////////

void MethodRegistry::insure_mstorage(const std::string& key)
{
  if (!m_storebytype.checkEntry(key))
    m_storebytype.addEntry( key, new QNameStorage() );
}

//////////////////////////////////////////////////////////////////////////////

std::vector< Common::SafePtr<Method> >
MethodRegistry::getAllMethods()
{
  std::vector< Common::SafePtr<Method> > all;

  MtdStorage::iterator itr = m_storebyqname.begin();
  std::transform(m_storebyqname.begin(),
                 m_storebyqname.end(),
                 back_inserter(all),
                 MtdStorage::extract);
  return all;
}

//////////////////////////////////////////////////////////////////////////////

std::vector< Common::SafePtr<Method> >
MethodRegistry::getAllMethods(const std::string& type, const std::string& namesp)
{
  insure_mstorage(type);

  QNameStorage * qnst = m_storebytype.getEntry(type);
  std::vector< Common::SafePtr<Method> > all;
  QNameStorage::iterator itr = qnst->begin();
  for (; itr != qnst->end(); ++itr)
  {
    if (itr->getNamespace() == namesp)
    {
      Method * ptr = m_storebyqname.getEntry(itr->str());
      all.push_back(ptr);
      CFLog(DEBUG_MED, "Returning method " << ptr->getName() << "\n");
    }
  }
  return all;
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
