// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_MethodCommandProvider_hh
#define COOLFluiD_Framework_MethodCommandProvider_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/CFLog.hh"
#include "Environment/ModuleRegister.hh"
#include "Framework/BaseMethodCommandProvider.hh"
#include "Framework/MethodCommand.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This class is used to create a command belonging to a Method.
/// It is an factory method pattern implementation.
/// @see BaseMethodCommandProvider
/// @see MethodCommand
/// @author Tiago Quintino
template <class COMMAND, class DATA, class MODULE>
class MethodCommandProvider : public BaseMethodCommandProvider<DATA, MethodCommand<DATA> > {
public:

  /// Constructor.
  /// @param name String defining the element type to be created
  explicit MethodCommandProvider(const std::string& name)
    : BaseMethodCommandProvider<DATA, MethodCommand<DATA> >(name)
  {
#ifndef CF_HAVE_SINGLE_EXEC   
    Environment::ModuleRegister<MODULE>::getInstance().getSelfRegistry().regist(this);
#endif 
  }

  /// Default destructor.
  ~MethodCommandProvider()
  {
//     Environment::ModuleRegister<MODULE>::getInstance().getSelfRegistry().unregist(this);
  }

  /// Creates a new CommandProvider with the supplied Data object.
  /// @param data the Data object.
  /// @return Pointer to the new CommandProvider
  Common::SelfRegistPtr< MethodCommand<DATA> > create(const std::string& name,
    				   const Common::SharedPtr<DATA>& data)
  {
    Common::SelfRegistPtr< MethodCommand<DATA> > ptr(new COMMAND(name), this);
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
    COMMAND* com = reinterpret_cast<COMMAND*>(ptr);
    cf_assert(com != CFNULL);
    delete com;
  }

}; // class MethodCommandProvider

//////////////////////////////////////////////////////////////////////////////

  } //  namespace Framework

} //  namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_MethodCommandProvider_hh
