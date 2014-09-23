// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Environment_Registry_hh
#define COOLFluiD_Environment_Registry_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/SafePtr.hh"
#include "Common/GeneralStorage.hh"
#include "Common/NonCopyable.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Environment {

//////////////////////////////////////////////////////////////////////////////

/// This class defines a General Registry class
/// @author Tiago Quintino
template < typename TYPE >
class Registry : public Common::NonCopyable<Registry> {

public:

  /// Constructor is private to allow only the friend classes to build it
  Registry() {}

  /// Default destructor is private to allow only the friend classes to destroy it
  ~Registry() {}

  /// Register an entry
  /// @param entry pointer to a TYPE to be added
  void regist(Common::TYPE* entry);

  /// Remove a registered entry
  /// @param entryName name of a TYPE to be removed
  void unregist(const std::string& entryName);

  /// Checks that a entry is registered
  /// @param entryName name of a TYPE to be checked
  bool isRegistered(const std::string& entryName);

  /// Get a given entry
  /// @param entryName name of a TYPE to be accessed
  /// @return a pointer to a TYPE if found or a null pointer if not found
  Common::SafePtr<Common::TYPE> getEntry(const std::string& entryName);

  /// Get all entrys
  /// @return a vector with all the entrys
  std::vector< Common::SafePtr<Common::TYPE> > getAllModules();

private: // data

  Common::GeneralStorage<Common::TYPE> m_store;

}; // end of class Registry

//////////////////////////////////////////////////////////////////////////////

  } // namespace Environment

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOFluiD_Environment_Registry_hh
