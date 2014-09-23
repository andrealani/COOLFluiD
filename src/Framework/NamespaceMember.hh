// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_NamespaceMember_hh
#define COOLFluiD_Framework_NamespaceMember_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/StringOps.hh"
#include "Common/SafePtr.hh"

#include "Framework/Framework.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

  class Namespace;

//////////////////////////////////////////////////////////////////////////////

/// This class represents a member of a Namespace
/// @author Tiago Quintino
class Framework_API NamespaceMember {
public:

  /// Default constructor without arguments
  NamespaceMember ( const std::string& name = "NULL" );

  /// Default destructor
  virtual ~NamespaceMember();

  /// Copy constructor based on operator=
  NamespaceMember ( const NamespaceMember& other) { *this = other; }

  /// Assignment operator=
  NamespaceMember& operator= ( const NamespaceMember& other);

  /// Accessor to the Namespace string
  /// @return CFstring with namespace
  std::string getNamespace() const { return m_namespace; }

  /// Retrive the Namespace ptr
  /// @return SafePtr pointing to the Namespace
  /// @throw Common::NoSuchStorageException if a namespace with this name does not exist
  Common::SafePtr<Namespace> getNamespacePtr() const;

  /// Mutator to the Namespace string to be used by the owner (parent) object
  void setParentNamespace(const std::string& name);

  /// Namespaces are set to this by default
  /// Later it should be modified by the owning Methods
  static std::string defaultNamespace() { return "NULL"; }

protected:

  /// Mutator to the Namespace
  /// Virtual to allow derived classes to perform specialized operations when namespace is changed
  virtual void setNamespace(const std::string& name);

  /// Mutator to the Namespace string to be used by the derived classes
  void setSelfNamespace(const std::string& name);

protected:

  /// the namespace
  std::string m_namespace;

private:

  /// indicates that the socket name should be the same as the parent
  bool m_sameAsParent;

}; // end of class NamespaceMember

//////////////////////////////////////////////////////////////////////////////

  } // namespace COOLFluiD

} // namespace Framework

//////////////////////////////////////////////////////////////////////////////

// #include "NamespaceMember.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_NamespaceMember_hh
