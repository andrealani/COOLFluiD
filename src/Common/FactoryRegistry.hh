// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Common_FactoryRegistry_hh
#define COOLFluiD_Common_FactoryRegistry_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/SafePtr.hh"
#include "Common/GeneralStorage.hh"
#include "Common/NonCopyable.hh"
#include "Common/Common.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Common {
    class FactoryBase;
    
//////////////////////////////////////////////////////////////////////////////

/// This class is a singleton object which serves as registry for all the
/// Factory objects that are created
///
/// @author Tiago Quintino
class Common_API FactoryRegistry : public NonCopyable<FactoryRegistry> {
 
 public:
  /// Default constructor
  FactoryRegistry();
  
  /// Default destructor 
  ~FactoryRegistry();
  
  /// Register a factory
  /// @param module pointer to a FactoryBase to be added
  void regist(FactoryBase* factory);
  
  /// Remove a registered factory
  /// @param name name of a FactoryBase to be removed
  void unregist(const std::string& name);
  
  /// Get a given factory
  /// @param name name of a FactoryBase to be accessed
  /// @return a pointer to a FactoryBase if found or a null pointer if not found
  SafePtr<FactoryBase> getFactory(const std::string& name);
  
 private: // data
  
  GeneralStorage<FactoryBase> m_store;

}; // end of class FactoryRegistry

//////////////////////////////////////////////////////////////////////////////

  } // namespace Common

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
  
#endif // COOFluiD_Common_FactoryRegistry_hh
