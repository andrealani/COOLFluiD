// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Environment_ObjectProvider_hh
#define COOLFluiD_Environment_ObjectProvider_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/COOLFluiD.hh"
#include "Common/CFLog.hh"

#include "Common/SelfRegistPtr.hh"
#include "Environment/ModuleRegister.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace Environment {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a concrete templated provider to create
/// whatever kind of polymorphic objects. By default it takes two
/// template parameters and it builds a non configurable object
/// whose create() method accepts no arguments.
/// @author Andrea Lani
/// @author Tiago Quintino
template <class CONCRETE, class BASE, class MODULE, int NBARGS = 0>
class ObjectProvider : public BASE::PROVIDER {
public:

  /// Constructor
  explicit ObjectProvider(const std::string& name) : BASE::PROVIDER(name)
  {
#ifdef CF_HAVE_LOG4CPP
    CFtrace << "Creating provider [" << name << "] of type [" << BASE::getClassName() << "]\n";
#endif
#ifndef CF_HAVE_SINGLE_EXEC 
   Environment::ModuleRegister<MODULE>::getInstance().getSelfRegistry().regist(this);
#endif
  }

  /// Destructor
  /// @todo unregistration is currently not working
  ~ObjectProvider()
  {
    // Environment::ModuleRegister<MODULE>::getInstance().getSelfRegistry().unregist(this);
  }

  /// Polymorphic function to create objects of dynamical type BASE
  /// @return SelfRegistPtr olding the created object
  Common::SelfRegistPtr<BASE> create()
  {
    return Common::SelfRegistPtr<BASE>(new CONCRETE(), this);
  }

  /// Polymorphic function to create objects of dynamical type BASE
  /// @return SelfRegistPtr olding the created object
  Common::SelfRegistPtr<BASE>* createPtr()
  {
   return new Common::SelfRegistPtr<BASE>(new CONCRETE(), this);
  }

  /// Free an instance created by this factory
  /// @param ptr pointer to be freed
  void freeInstance ( void* ptr )
  {
    cf_assert(ptr != CFNULL);
    CONCRETE* obj = reinterpret_cast<CONCRETE*>(ptr);

    cf_assert(obj != CFNULL);
    deletePtr<CONCRETE>( obj );
  }

}; // end of class ObjectProvider

//////////////////////////////////////////////////////////////////////////////

/// This class represents a concrete templated provider to create
/// whatever kind of polymorphic objects. It takes four
/// template parameters and it builds a non configurable object
/// whose create() method accepts one arguments.
/// @author Andrea Lani
/// @author Tiago Quintino
template <class CONCRETE, class BASE, class MODULE>
class ObjectProvider<CONCRETE,BASE,MODULE, 1> : public BASE::PROVIDER {
public:

  /// Constructor
  explicit ObjectProvider(const std::string& name) : BASE::PROVIDER(name)
  {
#ifdef CF_HAVE_LOG4CPP
    CFtrace << "Creating provider \'" << name << "\' of type \'" << BASE::getClassName() << "\'\n";
#endif
#ifndef CF_HAVE_SINGLE_EXEC  
    MODULE::getInstance().getSelfRegistry().regist(this);
#endif   
  }

  /// Destructor
  /// @todo unregistration is currently not working
  ~ObjectProvider()
  {
    // Environment::ModuleRegister<MODULE>::getInstance().getSelfRegistry().unregist(this);
  }

  /// Polymorphic function to create objects of dynamical type BASE
  /// @param arg1 first parameter
  /// @return SelfRegistPtr olding the created object
  Common::SelfRegistPtr<BASE> create(typename BASE::ARG1 arg)
  {
    CFtrace << "ObjectProvider::create() => 1" << BASE::getClassName() << "\n";
    CONCRETE* obj = new CONCRETE(arg); 
    CFtrace << "ObjectProvider::create() => 2" << BASE::getClassName() << "\n";
    Common::SelfRegistPtr<BASE> ptr (obj, dynamic_cast<Common::ProviderBase*>(this));
    CFtrace << "ObjectProvider::create() => ptr created\n";
    return ptr;
    // return Common::SelfRegistPtr<BASE>(new CONCRETE(arg), this);
  }

  /// Polymorphic function to create objects of dynamical type BASE
  /// @param arg1 first parameter
  /// @return SelfRegistPtr olding the created object
  Common::SelfRegistPtr<BASE>* createPtr(typename BASE::ARG1 arg)
  {
    CONCRETE* obj = new CONCRETE(arg);
    return new Common::SelfRegistPtr<BASE>(obj, dynamic_cast<Common::ProviderBase*>(this));
  }

  /// Free an instance created by this factory
  /// @param ptr pointer to be freed
  void freeInstance ( void* ptr )
  {
    cf_assert(ptr != CFNULL);
    CONCRETE* obj = reinterpret_cast<CONCRETE*>(ptr);

    cf_assert(obj != CFNULL);
    deletePtr<CONCRETE>( obj );
  }

}; // end of class ObjectProvider

//////////////////////////////////////////////////////////////////////////////

/// This class represents a concrete templated provider to create
/// whatever kind of polymorphic objects. It takes four
/// template parameters and it builds a non configurable object
/// whose create() method accepts one arguments.
/// @author Andrea Lani
/// @author Tiago Quintino
template <class CONCRETE, class BASE, class MODULE>
class ObjectProvider<CONCRETE,BASE,MODULE, 2> : public BASE::PROVIDER {
public:

  /// Constructor
  explicit ObjectProvider(const std::string& name) :
    BASE::PROVIDER(name)
  {
#ifdef CF_HAVE_LOG4CPP
    CFtrace << "Creating provider \'" << name << "\' of type \'" << BASE::getClassName() << "\'\n";
#endif
#ifndef CF_HAVE_SINGLE_EXEC  
   MODULE::getInstance().getSelfRegistry().regist(this);
#endif	 
  }

  /// Destructor
  /// @todo unregistration is currently not working
  ~ObjectProvider()
  {
    // Environment::ModuleRegister<MODULE>::getInstance().getSelfRegistry().unregist(this);
  }

  /// Polymorphic function to create objects of dynamical type BASE
  /// @param arg1 first parameter
  /// @param arg2 first parameter
  /// @return SelfRegistPtr olding the created object
  Common::SelfRegistPtr<BASE> create(typename BASE::ARG1 arg1,
                                    typename BASE::ARG2 arg2)
  {
    return Common::SelfRegistPtr<BASE>(new CONCRETE(arg1, arg2), this);
  }

  /// Polymorphic function to create objects of dynamical type BASE
  /// @param arg1 first parameter
  /// @param arg2 first parameter
  /// @return SelfRegistPtr olding the created object
  Common::SelfRegistPtr<BASE>* createPtr(typename BASE::ARG1 arg1,
                                         typename BASE::ARG2 arg2)
  {
    return new Common::SelfRegistPtr<BASE>(new CONCRETE(arg1, arg2), this);
  }

  /// Free an instance created by this factory
  /// @param ptr pointer to be freed
  void freeInstance ( void* ptr )
  {
    cf_assert(ptr != CFNULL);
    CONCRETE* obj = reinterpret_cast<CONCRETE*>(ptr);

    cf_assert(obj != CFNULL);
    deletePtr<CONCRETE>( obj );
  }

}; // end of class ObjectProvider

//////////////////////////////////////////////////////////////////////////////

  } // namespace Environment
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Environment_ObjectProvider_hh
