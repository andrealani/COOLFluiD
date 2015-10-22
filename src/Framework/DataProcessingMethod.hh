// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_DataProcessingMethod_hh
#define COOLFluiD_Framework_DataProcessingMethod_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/Method.hh"
#include "Environment/ConcreteProvider.hh"
#include "Framework/MultiMethodHandle.hh"
#include "Common/SafePtr.hh"
#include "Framework/SpaceMethod.hh"
#include "Framework/ConvergenceMethod.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

    class LinearSystemSolver;

//////////////////////////////////////////////////////////////////////////////

/// This class represents a DataProcessingMethod.
/// @author Thomas Wuilbaut
class Framework_API DataProcessingMethod : public Framework::Method,
                    public Common::DynamicFunctionCaller<DataProcessingMethod> {

public: // typedefs

  /// Type for the provider of this abstract class
  typedef Environment::ConcreteProvider<DataProcessingMethod,1> PROVIDER;
  typedef const std::string& ARG1;

public: // methods

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Default constructor without arguments
  DataProcessingMethod(const std::string& name);

  /// Default destructor
  virtual  ~DataProcessingMethod();

  /// Configures this Method.
  /// @param args the arguments used for the configuration
  virtual void configure( Config::ConfigArgs& args );
  
  /// Execute the data processing on the trs
  bool needsInitialization() {return m_needsInitialization;}
  
  /// Execute the data processing on the trs
  void processData();

  /// Clear the commands
  void clearComs();
    
  /// Sets the LinearSystemSolver
  void setCollaborator(MultiMethodHandle<LinearSystemSolver> lssMtd)
  {
    m_lssMtd = lssMtd;
  }
  
  /// Sets the ConvergenceMethod
  void setCollaborator(MultiMethodHandle<ConvergenceMethod> convergenceMtd)
  {
    m_convergenceMtd = convergenceMtd;
  }

  /// Checks if the ConvergenceMethod for this DataProcessing has been set.
  bool isConvergenceMethodSet() const
  {
    return (m_convergenceMtd.isNotNull());
  }

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
    Common::DynamicFunctionCaller<DataProcessingMethod>::run_dynamic_function(func);
  }

  /// Gets the Class name
  static std::string getClassName()
  {
    return "DataProcessing";
  }

protected: // functions

  /// Adds the ActionListener's of this EventListener to the EventHandler
  virtual void registActionListeners();

  /// Declares which functions can be called dynamically
  virtual void build_dynamic_functions();

protected: // abstract interface implementations

  /// Execute the data processing on the trs
  /// This is the abstract function that the concrete methods must implement.
  virtual void processDataImpl() = 0;

protected: // helper functions

  /// Gets the ConvergenceMethod which this DataProcessing uses
  /// @return pointer to the ConvergenceMethod
  MultiMethodHandle<ConvergenceMethod> getConvergenceMethod() const
  {
    return m_convergenceMtd;
  }

  /// Gets the LinearSystemSolver which this DataProcessing uses
  /// @return pointer to the LinearSystemSolver
  MultiMethodHandle<LinearSystemSolver> getLinearSystemSolver() const
  {
    return m_lssMtd;
  }

protected: // member data

  /// Rate to process the data.
  CFuint  m_processRate;

  /// Iteration at which processing starts
  CFuint  m_startIter;

  /// Iteration at which processing stops
  CFuint  m_stopIter;
  
  /// Flag telling if this processing needs initialization
  bool  m_needsInitialization;
  
private: // member data

  /// Convergence Method
  MultiMethodHandle<ConvergenceMethod> m_convergenceMtd;

  /// Linear system solver
  MultiMethodHandle<LinearSystemSolver> m_lssMtd;

}; // class DataProcessingMethod

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

Framework_Factory(DataProcessingMethod) // define the factory instance

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_DataProcessingMethod_hh
