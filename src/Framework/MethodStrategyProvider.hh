// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_MethodStrategyProvider_hh
#define COOLFluiD_Framework_MethodStrategyProvider_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/CFLog.hh"
#include "Environment/ModuleRegister.hh"
#include "Framework/BaseMethodStrategyProvider.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This class is used to create a Strategy belonging to a Method.
/// It is an factory method pattern implementation.
/// @see BaseMethodStrategyProvider
/// @see MethodStrategy
/// @author Tiago Quintino
template <class STRATEGY, class DATA, class BASESTRATEGY, class MODULE>
class MethodStrategyProvider :
  public BaseMethodStrategyProvider< DATA, BASESTRATEGY > {

public:

  /// Constructor.
  /// @param name String defining the element type to be created
  explicit MethodStrategyProvider(const std::string& name)
    : BaseMethodStrategyProvider<DATA, BASESTRATEGY >(name)
  {
#ifndef CF_HAVE_SINGLE_EXEC 
    Environment::ModuleRegister<MODULE>::getInstance().getSelfRegistry().regist(this);
#endif
  }

  /// Default destructor.
  ~MethodStrategyProvider()
  {
//     Environment::ModuleRegister<MODULE>::getInstance().getSelfRegistry().unregist(this);
  }

  /// Creates a new StrategyProvider with the supplied Data object.
  /// @param name name of the strategy
  /// @param data the method data object.
  /// @return Pointer to the new StrategyProvider
  Common::SelfRegistPtr< BASESTRATEGY > create (const std::string& name,
                                                const Common::SharedPtr<DATA>& data)
  {
    Common::SelfRegistPtr< BASESTRATEGY > ptr(new STRATEGY(name), this);
    ptr->setMethodData(data);
    return ptr;
  }

  /// Creates a new StrategyProvider with the supplied Data object.
  /// @param name name of the strategy
  /// @param arg1 extra parameter
  /// @param data the method data object.
  /// @return Pointer to the new StrategyProvider
  template < typename ARG1 >
  Common::SelfRegistPtr< BASESTRATEGY > create (const std::string& name,
                                                const ARG1& arg1,
                                                const Common::SharedPtr<DATA>& data)
  {
    Common::SelfRegistPtr< BASESTRATEGY > ptr(new STRATEGY(name,arg1), this);
    ptr->setMethodData(data);
    return ptr;
  }

  /// Creates a new StrategyProvider with the supplied Data object.
  /// @param name name of the strategy
  /// @param arg1 extra parameter
  /// @param data the method data object.
  /// @return Pointer to the new StrategyProvider
  template < typename ARG1, typename ARG2 >
  Common::SelfRegistPtr< BASESTRATEGY > create (const std::string& name,
                                                const ARG1& arg1,
                                                const ARG2& arg2,
                                                const Common::SharedPtr<DATA>& data)
  {
    Common::SelfRegistPtr< BASESTRATEGY > ptr(new STRATEGY(name,arg1,arg2), this);
    ptr->setMethodData(data);
    return ptr;
  }

  /// Free an instance created by this factory.
  /// (warning: If a provider is capable of instantiating multiple kinds of
  /// objects, the freeInstance method should make sure to call the right
  /// delete operator!)
  void freeInstance (void* ptr)
  {
    cf_assert(ptr != CFNULL);
    STRATEGY* com = reinterpret_cast<STRATEGY*>(ptr);
    cf_assert(com != CFNULL);
    delete com;
  }

}; // class MethodStrategyProvider

//////////////////////////////////////////////////////////////////////////////

  } //  namespace Framework

} //  namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_MethodStrategyProvider_hh
