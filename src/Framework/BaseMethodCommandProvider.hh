// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_BaseMethodCommandProvider_hh
#define COOLFluiD_Framework_BaseMethodCommandProvider_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/StringOps.hh"
#include "Common/SharedPtr.hh"
#include "Common/SelfRegistPtr.hh"
#include "Environment/Provider.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This class is used  as the base class for Providers who create a Commands
/// belonging to a Method.
/// It is an factory method pattern implementation.
/// This class in an abstract class.
/// @see Provider
/// @see CommandProvider
/// @author Tiago Quintino
/// @author Andrea Lani
template <class DATA, class BASECOMMAND>
class BaseMethodCommandProvider : public Environment::Provider<BASECOMMAND> {

public:

  /// Constructor
  /// @param name String defining the element type to be created
  explicit BaseMethodCommandProvider(const std::string& name)
    : Environment::Provider<BASECOMMAND>(name)
  {
  }

  /// Default destructor
  virtual ~BaseMethodCommandProvider()
  {
  }

  /// Creates a new MethodCommand with the supplied Data object.
  /// @param data the Data object.
  /// @return Pointer to the new MethodCommand
  Common::SelfRegistPtr<BASECOMMAND> create(const Common::SharedPtr<DATA>& data)
  {
    return create(BaseMethodCommandProvider
    <DATA, BASECOMMAND>::getName(), data);
  }

  /// Creates a new MethodCommand with the supplied Data object.
  /// @param name the name of the object to be used in configuration.
  /// @param data the Data object.
  /// @return Pointer to the new MethdCommand
  virtual Common::SelfRegistPtr<BASECOMMAND> create(const std::string& name,
      const Common::SharedPtr<DATA>& data) = 0;

    /// Free an instance created by this factory.
    /// (warning: If a provider is capable of instantiating multiple kinds of
    /// objects, the freeInstance method should make sure to call the right
    /// delete operator!)
  virtual void freeInstance (void* ptr) = 0;

}; // class BaseMethodCommandProvider

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_BaseMethodCommandProvider_hh
