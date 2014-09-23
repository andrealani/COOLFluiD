// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_DynamicBalancerMethod_hh
#define COOLFluiD_Framework_DynamicBalancerMethod_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/ConcreteProvider.hh"

#include "Framework/Storage.hh"
#include "Framework/MultiMethodHandle.hh"
#include "Framework/Method.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

    class DynamicBalancerMethodData;

//////////////////////////////////////////////////////////////////////////////

/// This class represents a DynamicBalancerMethod.
/// @author
class Framework_API DynamicBalancerMethod : public Method,
                          public Common::DynamicFunctionCaller<DynamicBalancerMethod> {

public: // typedefs

  /// Type for the provider of this abstract class
  typedef Environment::ConcreteProvider<DynamicBalancerMethod,1> PROVIDER;
  typedef const std::string& ARG1;

public: // static methods

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Gets the Class name
  static std::string getClassName()
  {
    return "DynamicBalancerMethod";
  }

public: // std methods

  /// Default constructor without arguments
  DynamicBalancerMethod(const std::string& name);

  /// Default destructor
  virtual ~DynamicBalancerMethod();

  /// Configures this Method.
  /// @param args the arguments used for the configuration
  virtual void configure ( Config::ConfigArgs& args );

public: // methods

  /// Does dynamic load balancing
  void doDynamicBalance();

  /// Run the function defined by the function name
  /// @param func name of the function to run. It should be void function with nor parameters.
  virtual void run_function(const std::string & func)
  {
    Common::DynamicFunctionCaller<DynamicBalancerMethod>::run_dynamic_function(func);
  }

protected: // interface implementation functions

  /// Sets the data of the method.
  /// @see Method::setMethod()
  virtual void setMethodImpl();

  /// Sets the data of the method.
  /// @see Method::setMethod()
  virtual void unsetMethodImpl();

  /// Do Repartitioning
  virtual void doDynamicBalanceImpl() = 0;

protected: // functions

  /// Adds the ActionListener's of this EventListener to the EventHandler
  virtual void registActionListeners();

  /// Declares which functions can be called dynamically
  virtual void build_dynamic_functions();

private: // member data

}; // end DynamicBalancerMethod

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

Framework_Factory(DynamicBalancerMethod) // define the factoty instance

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_DynamicBalancerMethod_hh
