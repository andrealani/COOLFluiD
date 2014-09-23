// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_NamespaceGroup_hh
#define COOLFluiD_Framework_NamespaceGroup_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/StringOps.hh"
#include "Framework/Framework.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a list of Namespaces
/// @author Tiago Quintino
class Framework_API NamespaceGroup {
public:

  /// storage type inside this list
  typedef std::vector<std::string> StorageType;

  /// constante iterastor to the storage inside this list
  typedef std::vector<std::string>::const_iterator NamespaceGroupIterator;

public:

  /// Default constructor
  /// @param name the primary namespace name
  NamespaceGroup(const std::string& name);

  /// Default destructor
  ~NamespaceGroup();

  /// Add a namespace to the list
  /// @param name the name of the namespace to try to add
  /// @throw DuplicateNameException if it does not exist
  void addNamespace(const std::string& name);

  /// Remove a namespace from the list
  /// @param name the name of the namespaces to try to remove
  /// @throw NoSuchValueException if it does not exist
  void removeNamespace(const std::string& name);

  /// Check if the supplied namespace exists in the list
  /// @param name the name of the namespace to try to match
  /// @return false if it does not exist
  bool match(const std::string& name);

  /// Remove all namespaces from the list
  void removeAll();

  /// Get Primary Namespace
  std::string getPrimaryNamespace();

protected:

  /// the storage of the list of namespaces
  StorageType _namespaceList;

}; // end of class NamespaceGroup

//////////////////////////////////////////////////////////////////////////////

  } // namespace COOLFluiD

} // namespace Framework

//////////////////////////////////////////////////////////////////////////////

// #include "NamespaceGroup.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_NamespaceGroup_hh
