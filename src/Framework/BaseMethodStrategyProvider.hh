// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_BaseMethodStrategyProvider_hh
#define COOLFluiD_Framework_BaseMethodStrategyProvider_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/StringOps.hh"
#include "Common/SharedPtr.hh"
#include "Common/SelfRegistPtr.hh"
#include "Environment/Provider.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This class is used  as the base class for Providers who create a Strategy's
/// belonging to a Method.
/// It is an factory method pattern implementation.
/// This class in an abstract class.
/// @see Provider
/// @see StrategyProvider
/// @author Tiago Quintino
template <class DATA, class BASESTRATEGY>
class BaseMethodStrategyProvider : public Environment::Provider<BASESTRATEGY> {
public:

  /// Constructor
  /// @param name String defining the element type to be created
  explicit BaseMethodStrategyProvider(const std::string& name)
    : Environment::Provider<BASESTRATEGY>(name)
  {
  }

  /// Default destructor
  virtual ~BaseMethodStrategyProvider()
  {
  }

  /// Creates a new MethodStrategy with the supplied Data object.
  /// @param data the Data object.
  /// @return Pointer to the new MethodStrategy
  Common::SelfRegistPtr<BASESTRATEGY> create(const Common::SharedPtr<DATA>& data)
  {
    return create(BaseMethodStrategyProvider <DATA, BASESTRATEGY>::getName(), data);
  }

  /// Creates a new MethodStrategy with the supplied Data object.
  /// @param name the name of the object to be used in configuration.
  /// @param data the Data object.
  /// @return Pointer to the new MethdStrategy
  virtual Common::SelfRegistPtr<BASESTRATEGY> create(const std::string& name,
                                                 const Common::SharedPtr<DATA>& data) = 0;

  /// Free an instance created by this factory
  /// @param ptr pointer to be freed
  virtual void freeInstance (void* ptr) = 0;

}; // class BaseMethodStrategyProvider

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_BaseMethodStrategyProvider_hh
