// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_NamespaceSwitcher_hh
#define COOLFluiD_Framework_NamespaceSwitcher_hh

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

/// This class is responsible to switch the Active Namespace
/// upon request by the Method's.
/// This is a Singleton object.
/// @author Tiago Quintino
class Framework_API NamespaceSwitcher : public Common::NonCopyable<NamespaceSwitcher> {
public:
  
  /// Returns the instance of the NamespaceSwither Singleton
  static NamespaceSwitcher& getInstance(const std::string& subSystemName);
  
  /// Creates a Namespace with a unique name and puts it on the internal
  ///  storage
  /// @param name the name of the Namespace to create
  Common::SafePtr<Namespace> createUniqueNamespace(const std::string& name);

  /// Removes and deletes all the Namespaces
  void deleteAllNamespaces();

  /// Pushs a Namespace identified by the name into the stack
  /// @param name the name of the Namespace to push into the stack
  void pushNamespace(const std::string& name);

  /// Gets a Namespace identified by the name
  /// @param name the name of the Namespace to get
  /// @return SafePtr to the Namespace
  Common::SafePtr<Namespace> getNamespace(const std::string& name);

  /// Pops the top Namespace from the stack
  /// @pre stack is not empty
  /// @return SafePtr to the Namespace that has been poped from the stack
  Common::SafePtr<Namespace> popNamespace();

  /// Gets the top Namespace from the stack
  /// @pre stack is not empty
  /// @return SafePtr to the Namespace that is on the top of the stack
  Common::SafePtr<Namespace> getCurrentNamespace();

  /// Gets the list of Namespace's in this SubSystem
  /// @return vector with pointers to all Namespace's
  std::vector<Common::SafePtr<Namespace> > getAllNamespaces();
  
  /// Set the flag to say that the objects of
  /// the stack can be accessed
  void setEnabled(bool isEnabled);
  
  /// Get the name of the given object belonging to the namespace
  /// @param filterCoupling  flag telling whether to discard coupling namespaces
  std::string getName(std::const_mem_fun_t<std::string, Namespace> fun,
		      bool filterCoupling);
  
  /// Get the ID corresponding to the namespace 
  /// @param filterCoupling  flag telling whether to discard coupling namespaces
  CFuint getID(const bool filterCoupling);
  
private:

  /// Default constructor without arguments
  NamespaceSwitcher();

  /// Destructor
  ~NamespaceSwitcher();

protected: // helper functions

  /// Check if the objects of the stack can be accessed
  bool isEnabled();

private: // member data

  /// Storage of the Namespaces
  Common::GeneralStorage<Namespace> m_Storage;

  /// stack of the Namespace pointers
  std::stack<Namespace*> m_Stack;

  /// flag to determine if objects of the stack can be accessed
  bool m_isEnabled;

}; // end of class NamespaceSwitcher

//////////////////////////////////////////////////////////////////////////////

  } // namespace COOLFluiD

} // namespace Framework

//////////////////////////////////////////////////////////////////////////////

// #include "NamespaceSwitcher.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_NamespaceSwitcher_hh
