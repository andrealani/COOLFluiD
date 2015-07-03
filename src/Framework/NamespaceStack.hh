// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_NamespaceStack_hh
#define COOLFluiD_Framework_NamespaceStack_hh

//////////////////////////////////////////////////////////////////////////////

#include <stack>

#include "Common/SafePtr.hh"
#include "Common/NonCopyable.hh"
#include "Common/GeneralStorage.hh"

#include "Framework/Namespace.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This class represents an object that holds a stack of objects of
/// the same TYPE which are switched according to the active Namespace.
/// @see NamespaceSwitcher
/// @author Tiago Quintino
template < typename TYPE >
class NamespaceStack : public Common::NonCopyable< NamespaceStack<TYPE> > {

public: // methods

  /// Default constructor without arguments
  NamespaceStack();

  /// Destructor
  virtual ~NamespaceStack();

  /// Set the flag to say that the objects of
  /// the stack can be accessed
  void setEnabled(bool isEnabled);

  /// Returns a pointer to the entry associated to the supplied name
  /// @param name the name of the object to get
  /// @throw Common::NoSuchStorageException if the entry does not exist in the storage
  /// @return SafePtr to the object of the stack
  Common::SafePtr<TYPE> getEntry(const std::string& name);

  /// Returns a pointer to the entry associated to the supplied name
  /// @param name the name of the namespace to which the object is associated
  /// @throw Common::NoSuchStorageException if the namespace does not exist in the storage
  /// @return SafePtr to the object of the stack
  Common::SafePtr<TYPE> getEntryByNamespace(const Common::SafePtr<Namespace>& nsp);

  /// Creates an object with a unique name and puts it on the internal
  ///  storage
  /// @param name the name of the object to create
  /// @return SafePtr to the object created
  Common::SafePtr<TYPE> createUnique(const std::string& name);

  /// Deletes all entries in the Stack and Storage
  void deleteAllEntries();
 
  /// Deletes all entries in the Stack and Storage
  void deleteEntry(const std::string name);

  /// Gets all entries in the Stack
  std::vector<Common::SafePtr<TYPE> > getAllEntries();

  /// Pushs an object defined by the Namespace
  /// into the stack
  /// @param name the name of the Namespace to push into the stack
  /// @throw Common::NoSuchStorageException if there is no object as defined by the
  ///        namespace
  void push(const Common::SafePtr<Namespace>& nsp);

  /// Pops the top object from the stack
  /// @pre m_stack is not empty
  /// @return SafePtr to the object that has been poped from the stack
  Common::SafePtr<TYPE> pop();

  /// Returns the top object on the stack
  /// @pre m_stack is not empty
  /// @return SafePtr to the object that on the top of the stack
  Common::SafePtr<TYPE> top();

protected: // helper functions

  /// Gets the name of the object form the Namespace
  /// @param nsp the Namespace from where to get te object name
  virtual std::string
  getObjectName(const Common::SafePtr<Namespace>& nsp) const = 0;

  /// Creates an Object of the type
  /// @param name the object name
  virtual TYPE *
  createObject(const std::string& name) = 0;

  /// Check if the objects of the stack can be accessed
  bool isEnabled();

private: // member data

  /// Storage of the Objects
  Common::GeneralStorage<TYPE> m_storage;

  /// stack of the Object pointers
  std::stack<TYPE*> m_stack;

  /// flag to determine if objects of the stack can be accessed
  bool m_isEnabled;

}; // end of class NamespaceStack

//////////////////////////////////////////////////////////////////////////////

  } // namespace COOLFluiD

} // namespace Framework

//////////////////////////////////////////////////////////////////////////////

#include "NamespaceStack.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_NamespaceStack_hh
