// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/MethodData.hh"
#include "Framework/BadFormatException.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

CollaboratorAccess::CollaboratorAccess(Common::SafePtr<Method> owner) :
  m_owner(owner),
  m_CollaboratorNames(std::vector<std::string>())
{
}

//////////////////////////////////////////////////////////////////////////////

CollaboratorAccess::~CollaboratorAccess()
{
}

//////////////////////////////////////////////////////////////////////////////

void CollaboratorAccess::setup()
{
  // verify that all the collaborator names have methods associated
  std::vector<QualifiedName>::iterator itr = m_QNames.begin();
  for (; itr != m_QNames.end(); ++itr)
  {
    try
    {
      Common::SafePtr<Method> mtd = MethodRegistry::getInstance().getMethod(*itr);
    }
    catch (Common::Exception& e)
    {
      throw CollaboratorException (FromHere(),"Check to the collaborator names in method " + m_owner->getName() + " failed. Original problem was :" + e.str());
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void CollaboratorAccess::unsetup()
{
}

//////////////////////////////////////////////////////////////////////////////

void CollaboratorAccess::process()
{
CFAUTOTRACE;

  // build the qualifed name list
  CFLogDebugMin("Method " << m_owner->QName().str());
  CFLogDebugMin(" processing " << m_CollaboratorNames.size() << " collabs.\n");
  CFLogDebugMin("  Adding ... ");
  std::vector<std::string>::iterator itr = m_CollaboratorNames.begin();
  for (; itr != m_CollaboratorNames.end(); ++itr)
  {
    std::vector<std::string> words = Common::StringOps::getWords(*itr,QualifiedName::separator());
    if (words.size() == 1)
    {
      QualifiedName qname = QualifiedName(m_owner->getNamespace(),words[0]);
      m_QNames.push_back(qname);
      CFLogDebugMin(" " << qname.str());
    }
    else
    {
      if (words.size() == 2)
      {
        QualifiedName qname = QualifiedName(words[0],words[1]);
        m_QNames.push_back(qname);
        CFLogDebugMin(" " << qname.str());
      }
      else
      {
        throw BadFormatException (FromHere(),"Collaborator name is not in the correct format [namespace:method]. Found " + (*itr) );
      }
    }
  }
  CFLogDebugMin(" done\n");
}

//////////////////////////////////////////////////////////////////////////////

bool CollaboratorAccess::isCollaborator(const QualifiedName& qname)
{
  CFAUTOTRACE;

  const std::vector<QualifiedName>& collaborators = getCollaboratorQNames();

  std::vector<QualifiedName>::const_iterator collab = collaborators.begin();
  for(; collab != collaborators.end(); ++collab) {
    if(*collab == qname) return true;
  }

  return false;
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
