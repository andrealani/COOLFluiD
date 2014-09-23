// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_ErrorEstimatorMethod_hh
#define COOLFluiD_Framework_ErrorEstimatorMethod_hh

//////////////////////////////////////////////////////////////////////////////

#include "Method.hh"
#include "Environment/ConcreteProvider.hh"
#include "MultiMethodHandle.hh"
#include "SpaceMethod.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a ErrorEstimatorMethod.
/// The ErrorEstimatorMethod estimates the error in the provided solution.
/// @author Tiago Quintino
/// @author Jurek Majevski
class Framework_API ErrorEstimatorMethod : public Method,
                    public Common::DynamicFunctionCaller<ErrorEstimatorMethod> {
public: // typedefs

  /// Type for the provider of this abstract class
  typedef Environment::ConcreteProvider<ErrorEstimatorMethod,1> PROVIDER;
  typedef const std::string& ARG1;

public: // methods

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Default constructor without arguments
  ErrorEstimatorMethod(const std::string& name);

  /// Default destructor
  virtual ~ErrorEstimatorMethod();

  /// Configures this Method.
  /// @param args the arguments used for the configuration
  virtual void configure ( Config::ConfigArgs& args );

  /// Estimate the error
  void estimate();

  /// Sets the data of the method.
  /// @see Method::setMethod()
  virtual void setMethodImpl();

  /// Sets the data of the method.
  /// @see Method::setMethod()
  virtual void unsetMethodImpl();

  /// Run the function defined by the function name
  /// @param func name of the function to run. It should be void function with nor parameters.
  virtual void run_function(const std::string & func)
  {
    Common::DynamicFunctionCaller<ErrorEstimatorMethod>::run_dynamic_function(func);
  }

  /// Gets the Class name
  static std::string getClassName() { return "ErrorEstimatorMethod"; }

protected: // abstract interface implementations

  /// Estimate the error
  /// This is the abstract function that the concrete methods must implement.
  virtual void estimateImpl() = 0;

protected: // functions

  /// Adds the ActionListener's of this EventListener to the EventHandler
  virtual void registActionListeners();

  /// Declares which functions can be called dynamically
  virtual void build_dynamic_functions();

protected: // member data

  /// Estimation rate for the error estimator
  CFuint  _estimateRate;

}; // end ErrorEstimatorMethod

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

Framework_Factory(ErrorEstimatorMethod) // define the factoty instance

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_ErrorEstimatorMethod_hh
