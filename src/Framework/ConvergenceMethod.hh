// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_ConvergenceMethod_hh
#define COOLFluiD_Framework_ConvergenceMethod_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/Stopwatch.hh"

#include "Environment/ConcreteProvider.hh"

#include "Framework/Storage.hh"
#include "Framework/MultiMethodHandle.hh"
#include "Framework/CFL.hh"
#include "Framework/Method.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

    class SpaceMethod;
    class LinearSystemSolver;
    class ConvergenceMethodData;
    class StopConditionController;

//////////////////////////////////////////////////////////////////////////////

/// This class represents a ConvergenceMethod.
/// The ConvergenceMethod makes use of the SpaceMethod for discretizing the
/// spacial terms of the equation.
/// @todo Implement automatic monitoring of CFL for convergenceMethods.
/// @author Tiago Quintino
/// @author Andrea Lani
class Framework_API ConvergenceMethod :
    public Method,
    public Common::DynamicFunctionCaller<ConvergenceMethod> {

public: // typedefs

  /// Type for the provider of this abstract class
  typedef Environment::ConcreteProvider<ConvergenceMethod,1> PROVIDER;
  typedef const std::string& ARG1;

public: // methods

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Default constructor without arguments
  ConvergenceMethod(const std::string& name);

  /// Default destructor
  virtual ~ConvergenceMethod();

  /// Gets the Data aggregator of this Convergence method
  /// @return SafePtr to the ConvergenceMethodData
  virtual Common::SafePtr<Framework::ConvergenceMethodData> getConvergenceMethodData() = 0;

  /// Configures this Method.
  /// @param args the arguments used for the configuration
  virtual void configure ( Config::ConfigArgs& args );

  /// Take one timestep
  /// @post pushs and pops the Namespace to which this Method belongs
  void takeStep();

  /// Write the convergence information on screen
  void writeOnScreen();

  /// Sets the LinearSystemSolver
  void setCollaborator(MultiMethodHandle<LinearSystemSolver> lss)
  {
    m_lss = lss;
  }

  /// Sets the SpaceMethod
  virtual void setCollaborator(MultiMethodHandle<SpaceMethod> sm) {}
  
  /// Checks if the LinearSystemSolver for this ConvergenceMethod has been set.
  bool isLinearSystemSolverSet() const
  {
    return (m_lss.isNotNull());
  }

  /// Run the function defined by the function name
  /// @param func name of the function to run. It should be void function with nor parameters.
  virtual void run_function(const std::string & func)
  {
    Common::DynamicFunctionCaller<ConvergenceMethod>::run_dynamic_function(func);
  }

  /// Get the CFL
  Common::SafePtr<Framework::CFL> getCFL();

  /// Gets the Class name
  static std::string getClassName() { return "ConvergenceMethod";  }

  bool outputSpaceResidual() const {return m_outputSpaceResidual;}

protected: // interface implementation functions

  /// Set the convergence method
  virtual void setMethodImpl();

  /// Set the convergence method
  virtual void unsetMethodImpl();

  /// Take one timestep
  /// This is the abstract function that the concrete methods must implement.
  virtual void takeStepImpl() = 0;

protected: // helper functions

  /// function which indicates if we have to update the convergence file
  /// for this processor
  bool hasToUpdateConv();

  /// Gets the LinearSystemSolver which this ConvergenceMethod uses
  /// @return pointer to the ConvergenceMethod
  MultiMethodHandle<LinearSystemSolver>& getLinearSystemSolver()
  {
    return m_lss;
  }

  /// Syncronize the all and compute the residual
  void syncAllAndComputeResidual(const bool computeResidual);

  /// Syncronize the states and compute the residual
  void syncGlobalDataComputeResidual(const bool computeResidual);

  /// Prepare the convergence file
  void prepareConvergenceFile();

  /// Update the convergence file
  void updateConvergenceFile();
  
protected: // functions

  /// Adds the ActionListener's of this EventListener to the EventHandler
  virtual void registActionListeners();

  /// Declares which functions can be called dynamically
  virtual void build_dynamic_functions();
  
  /// @return the convergence file name
  std::string getFileName(const std::string name) const;
  
 protected: // member data

  /// configuration string for the StopCondition
  std::string m_stopConditionStr;

  /// controller for stoping the iterative process
  std::auto_ptr<Framework::StopConditionController> m_stopCondControler;

private: // member data

  /// LinearSystemSolver used to solve the linear system, if one is present
  MultiMethodHandle<LinearSystemSolver> m_lss;
  
  /// stopwatch to keep the track of the time spent converging
  Common::Stopwatch<Common::WallTime>  m_stopwatch;

  /// Name of Convergence File where to write the convergence history.
  std::string m_nameConvergenceFile;
  
  /// Name of Space Residual File where to write the space residual history.
  std::string m_nameSpaceResidualFile;

  /// Rate to save convegence info to convergence file.
  CFuint  m_convRate;

  /// Rate to show convergence message to stdout.
  CFuint  m_showRate;

  /// Precision for outputting to screen
  CFint m_precision;

  /// flag to indicate if only the processor 0 should update the convergence file
  bool m_onlyP0;
  
  /// flag to indicate if spatial residual must be computed separately (for implicit methods)
  bool m_outputSpaceResidual;
  
}; // end ConvergenceMethod

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

Framework_Factory(ConvergenceMethod) // define the factoty instance

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_ConvergenceMethod_hh
