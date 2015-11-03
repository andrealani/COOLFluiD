// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_CouplerMethod_hh
#define COOLFluiD_Framework_CouplerMethod_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/Method.hh"
#include "Framework/SpaceMethod.hh"
#include "Framework/MultiMethodHandle.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a CouplerMethod.
/// @author Thomas Wuilbaut
class Framework_API CouplerMethod : public Method,
                      public Common::DynamicFunctionCaller<CouplerMethod> {

public: // typedefs

  /// Type for the provider of this abstract class
  typedef Environment::ConcreteProvider<CouplerMethod,1> PROVIDER;
  typedef const std::string& ARG1;

public: // methods

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Default constructor without arguments
  CouplerMethod(const std::string& name);

  /// Default destructor
  virtual ~CouplerMethod();

  /// Configures this Method.
  /// @param args the arguments used for the configuration
  virtual void configure ( Config::ConfigArgs& args );

  /// Write info for the other subsystems
  virtual void setInfoToOtherSubSystem();

  /// Read info from the other subsystems
  virtual void getInfoFromOtherSubSystem();

  /// Finish configuring this Method.
  virtual void postConfigure( Config::ConfigArgs& args );

  /// Writes the data necessary for the preProcessing
  /// @post pushs and pops the Namespace to which this Method belongs
  void preProcessWrite();

  /// Reads the data necessary for the preProcessing
  /// @post pushs and pops the Namespace to which this Method belongs
  void preProcessRead();

  /// Writes the data necessary for the mesh matching
  /// @post pushs and pops the Namespace to which this Method belongs
  void meshMatchingWrite();

  /// Reads the data necessary for the mesh matching
  /// @post pushs and pops the Namespace to which this Method belongs
  void meshMatchingRead();

  /// Reads the data from the Coupled SubSystems
  /// @post pushs and pops the Namespace to which this Method belongs
  void dataTransferRead();

  /// Writes the data for the Coupled SubSystems
  /// @post pushs and pops the Namespace to which this Method belongs
  void dataTransferWrite();

  /// Finalize the coupling procedure
  /// @post pushs and pops the Namespace to which this Method belongs
  void finalize();

  /// Run the function defined by the function name
  /// @param func name of the function to run. It should be void function with nor parameters.
  virtual void run_function(const std::string & func)
  {
    Common::DynamicFunctionCaller<CouplerMethod>::run_dynamic_function(func);
  }

  /// Gets the Class name
  static std::string getClassName()
  {
    return "CouplerMethod";
  }

protected: // abstract interface implementations

  /// Sets the data of the method.
  /// @see Method::setMethod()
  virtual void setMethodImpl();

  /// Sets the data of the method.
  /// @see Method::setMethod()
  virtual void unsetMethodImpl();

  /// Writes the data necessary for the preProcessing
  /// This is the abstract function that the concrete methods must implement.
  virtual void preProcessWriteImpl() = 0;

  /// Reads the data necessary for the preProcessing
  /// This is the abstract function that the concrete methods must implement.
  virtual void preProcessReadImpl() = 0;

  /// Writes the data necessary for the mesh matching
  /// This is the abstract function that the concrete methods must implement.
  virtual void meshMatchingWriteImpl() = 0;

  /// Reads the data necessary for the mesh matching
  /// This is the abstract function that the concrete methods must implement.
  virtual void meshMatchingReadImpl() = 0;

  /// Writes the data for the Data Transfer
  /// This is the abstract function that the concrete methods must implement.
  virtual void dataTransferWriteImpl() = 0;

  /// Reads the data from the Coupled SubSystems
  /// This is the abstract function that the concrete methods must implement.
  virtual void dataTransferReadImpl() = 0;
  
  /// Finalize the coupling
  /// This is the abstract function that the concrete methods must implement.
  virtual void finalizeImpl() = 0;

protected: // functions

  /// Adds the ActionListener's of this EventListener to the EventHandler
  virtual void registActionListeners();

  /// Declares which functions can be called dynamically
  virtual void build_dynamic_functions();

}; // class CouplerMethod

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

Framework_Factory(CouplerMethod) // define the factoty instance

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_CouplerMethod_hh
