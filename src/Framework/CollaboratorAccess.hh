// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_CollaboratorAccess_hh
#define COOLFluiD_Framework_CollaboratorAccess_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/SafePtr.hh"
#include "Framework/CollaboratorException.hh"
#include "Framework/MethodRegistry.hh"
#include "Framework/Method.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This class takes manages the access to the collaborators of a Method.
/// These collaborators are themselves  Method's.
/// @author Tiago Quintino
class Framework_API CollaboratorAccess {

public: // functions

  /// Constructor
  CollaboratorAccess(Common::SafePtr<Method> owner);

  /// Default destructor
  ~CollaboratorAccess();

  /// Setup function
  void setup();

  /// Unsetup function
  void unsetup();

  /// Answers if the method described by its qualified name, namespace plus name is a collaborator.
  /// @param qualified_name is the namespace prepended to the name od the method "namespace::name".
  /// @return if the described Method is a collaborator of this Method
  bool isCollaborator(const QualifiedName& qname);

  /// Gets names of the collaborator objects needed for the mapping during set up
  const std::vector<QualifiedName>& getCollaboratorQNames() const
  {
    return m_QNames;
  }

  /// This function will return the collaborator of the method owning the object,
  /// if this collaborator is unique in TYPE.
  /// The namespace assumed is the current active one, which should be the same
  /// as the owner Method.
  /// @throw COllaboratorException if there no unique collaborator of this TYPE
  template < typename TYPE >
  Common::SafePtr <TYPE> getCollaborator();

  /// This function will return the all collaborators of the method owning the object
  /// with the TYPE specified.
  /// The namespace assumed is the current active one, which should be the same
  /// as the owner Method.
  template < typename TYPE >
  std::vector < Common::SafePtr <TYPE> > getAllCollaborators();

  /// This function will return the collaborator of the method owning the object,
  /// if this collaborator is unique in TYPE within the supplied namespace.
  /// @param namesp is the Namespace in which to search for the collaborator.
  /// @throw CollaboratorException if there no unique collaborator of this TYPE
  template < typename TYPE >
  Common::SafePtr <TYPE> getCollaborator(const std::string namesp);

  /// This function will return all the collaborators of the method owning the object.
  /// @param namesp is the Namespace in which to search for the collaborator.
  template < typename TYPE >
  std::vector < Common::SafePtr <TYPE> > getAllCollaborators(const std::string namesp);

  /// Acessor to the own method for MethodCommand's
  Common::SafePtr<Method> getOwnMethod() { return m_owner; }

protected: // helper method

  /// Processes the CollaboratorNames and creates the QualifiedName list
  void process ();

protected: // data

  /// list of QualifiedNames of collaborators
  /// is is postprocessed from the m_CollaboratorNames
  std::vector<QualifiedName> m_QNames;

  /// Acquaintance of the owning method
  Common::SafePtr<Method> m_owner;

  /// Names of the collaborators objects needed for the mapping during set up
  std::vector<std::string> m_CollaboratorNames;

}; // end of class CollaboratorAccess

//////////////////////////////////////////////////////////////////////////////

template < typename TYPE >
Common::SafePtr <TYPE> CollaboratorAccess::getCollaborator()
{
  typedef std::vector< Common::SafePtr< TYPE > > VecMet;
  VecMet all = MethodRegistry::getInstance().template getAllMethods< TYPE >();

  if (all.size() == 0)
    throw CollaboratorException (FromHere(),"From method " + m_owner->getName() + ", no collaborating methods of type " + TYPE::getClassName() + " were found in any namespace.");

  if ((all.size() == 1) && (all[0]->getNamespace() == m_owner->getNamespace()))
  {
    CFLog(DEBUG_MED, "Found only one Collaborator of type " + TYPE::getClassName() + " in same namespace as owner. Accessing collaborator " << all[0]->getName() << "\n");
    return all[0];
  }
  else
  {
    // although their are multiple methods of this type
    // check that there is only one collaborator defined by the user

    VecMet matches;
    typename VecMet::iterator itr = all.begin();
    for (; itr != all.end(); ++itr)
    {
      if ( std::count_if(m_QNames.begin(),
                         m_QNames.end(),
                         std::bind2nd( std::equal_to<QualifiedName>(), (*itr)->QName()) ) )
        matches.push_back(*itr);
    }
    if (matches.size() == 0)
      throw CollaboratorException (FromHere(),"From method " + m_owner->getName() +  ", no collaborator of type " + TYPE::getClassName() + " found in any namespace. Must define appropriate collaborator name.");

    if (matches.size() > 1)
      throw CollaboratorException (FromHere(),"From method " + m_owner->getName() +  ", non unique collaborator of type " + TYPE::getClassName() + " found across namespaces. Probably this method should only have one collaborator of this type.");

    CFLog(DEBUG_MED, "Found only one Collaborator of type " + TYPE::getClassName() + ".Accessing collaborator " << matches[0]->getName() << "\n");
    return matches[0];
  }
}

//////////////////////////////////////////////////////////////////////////////

template < typename TYPE >
Common::SafePtr <TYPE> CollaboratorAccess::getCollaborator(const std::string namesp)
{
  typedef std::vector< Common::SafePtr< TYPE > > VecMet;
  VecMet all = MethodRegistry::getInstance().template getAllMethods<TYPE>(namesp);

  if (all.size() == 0)
    throw CollaboratorException (FromHere(),"From method " + m_owner->getName() + ", no collaborating methods of type " + TYPE::getClassName() + " were found in namespace " + namesp);

  if (all.size() == 1)
  {
    CFLog(DEBUG_MED, "Found only one Collaborator of type " + TYPE::getClassName() + ". Accessing collaborator " << all[0]->getName() << "\n");
    return all[0];
  }
  else
  {
    // although their are multiple methods of this type
    // check that there is only one collaborator defined by the user

    VecMet matches;
    typename VecMet::iterator itr = all.begin();
    for (; itr != all.end(); ++itr)
    {
      if ( std::count_if(m_QNames.begin(),
                         m_QNames.end(),
                         std::bind2nd( std::equal_to<QualifiedName>(), (*itr)->QName()) ) )
        matches.push_back(*itr);
    }
    if (matches.size() == 0)
      throw CollaboratorException (FromHere(),"From method " + m_owner->getName() +  ", no collaborator of type " + TYPE::getClassName() + " found in namespace " + namesp + ". Must define appropriate collaborator name.");

    if (matches.size() > 1)
      throw CollaboratorException (FromHere(),"From method " + m_owner->getName() +  ", non unique collaborator of type " + TYPE::getClassName() + " found in namespace " + namesp + ". Probably this method should only have one collaborator of this type.");

    CFLog(DEBUG_MED, "Found only one Collaborator of type " + TYPE::getClassName() + ".Accessing collaborator " << matches[0]->getName() << "\n");
    return matches[0];
  }
}

//////////////////////////////////////////////////////////////////////////////

template < typename TYPE >
std::vector < Common::SafePtr <TYPE> >
CollaboratorAccess::getAllCollaborators()
{
  return getAllCollaborators<TYPE>(m_owner->getNamespace());
}

//////////////////////////////////////////////////////////////////////////////

template < typename TYPE >
std::vector < Common::SafePtr <TYPE> >
CollaboratorAccess::getAllCollaborators(const std::string namesp)
{
  typedef std::vector< Common::SafePtr< TYPE > > VecMet;

  VecMet all = MethodRegistry::getInstance().template getAllMethods< TYPE >(namesp);

  VecMet matches;
  typename VecMet::iterator itr = all.begin();
  for (; itr != all.end(); ++itr)
  {
    CFuint nbmatches = std::count_if(m_QNames.begin(),
                                     m_QNames.end(),
                                     std::bind2nd(std::equal_to<QualifiedName>(),(*itr)->QName()));
    if ( nbmatches > 1)
      throw CollaboratorException (FromHere(),"Duplicate name of method " + (*itr)->QName().getName() + " found in Collaborators names.");

    if ( nbmatches == 1)
    {
      CFLog(DEBUG_MIN, "Accessing collaborator " << (*itr)->getName() << "\n");
      matches.push_back(*itr);
    }
  }

  if (matches.size() == 0)
    throw CollaboratorException (FromHere(),"No collaborator of type " + TYPE::getClassName() + " found. Must define appropriate collaborators names.");

  return matches;
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOFluiD_Framework_CollaboratorAccess_hh
