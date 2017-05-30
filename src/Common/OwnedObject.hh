// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Common_OwnedObject_hh
#define COOLFluiD_Common_OwnedObject_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/COOLFluiD.hh"
#include "Common/Common.hh"
#include "Common/FactoryRegistry.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Common {

//////////////////////////////////////////////////////////////////////////////

/// Base class for objects that will be held by SharedPtr or SelfRegistPtr
/// @see SharedPtr
/// @see SelfRegistPtr
/// @author Andrea Lani
/// @author Tiago Quintino

class Common_API OwnedObject {

 public: // functions
  
  /// Default constructor
  OwnedObject() : m_fr(CFNULL), m_owners(0) {}
  
  /// Virtual destructor
  virtual ~OwnedObject() {}
  
  /// Set the factory registry
  virtual void setFactoryRegistry(Common::SafePtr<Common::FactoryRegistry> fr) 
  {m_fr = fr;}
  
  /// Get the factory registry
  Common::SafePtr<Common::FactoryRegistry> getFactoryRegistry() 
  {
#ifdef CF_HAVE_SINGLE_EXEC
    cf_assert(m_fr != CFNULL);
#endif
    return m_fr;
  }  
  
  /// Add an owner.
  void addOwner() { m_owners++; }

  /// Remove an owner.
  void removeOwner() { m_owners--; }

  /// Check if the object is not owned anymore.
  bool hasNoOwner() { return (m_owners == 0) ? true : false; }

  /// Get the number of owners
  CFuint getNbOwners() const { return m_owners; }

private: // data
  
  /// factory registry to allow polymorphic creation of objects
  SafePtr<FactoryRegistry> m_fr;
  
  /// Counter for the number of owners of this object
  CFuint m_owners;
  
}; // end class OwnedObject

//////////////////////////////////////////////////////////////////////////////

  } // namespace Common

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Common_OwnedObject_hh
