// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_SpaceMethod_hh
#define COOLFluiD_Framework_SpaceMethod_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/SafePtr.hh"

#include "Environment/ConcreteProvider.hh"

#include "Framework/Method.hh"
#include "Framework/MultiMethodHandle.hh"

#include "MathTools/RealVector.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

    class LinearSystemSolver;
    class ConvergenceMethod;
    class SpaceMethodData;
    class VolumeIntegrator;
    class MeshDataBuilder;
    class GlobalJacobianSparsity;

//////////////////////////////////////////////////////////////////////////////

/// This class represents a SpaceMethod.
/// @author Tiago Quintino
/// @author Andrea Lani
class Framework_API SpaceMethod :
    public Method,
    public Common::DynamicFunctionCaller<SpaceMethod> {

public: // typedefs

  /// Type for the provider of this abstract class
  typedef Environment::ConcreteProvider<SpaceMethod,1> PROVIDER;
  typedef const std::string& ARG1;

public: // methods

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Default constructor without arguments
  SpaceMethod(const std::string& name);

  /// Default destructor
  virtual  ~SpaceMethod();

  /// Configures this Method.
  /// @param args the arguments used for the configuration
  virtual void configure( Config::ConfigArgs& args );

  /// Gets the Data aggregator of this space method
  /// @return SafePtr to the SpaceMethodData
  virtual Common::SafePtr<Framework::SpaceMethodData> getSpaceMethodData() = 0;

  /// Sets the LinearSystemSolver for this SpaceMethod to use
  /// @pre the pointer to LinearSystemSolver is not constant to
  ///      allow dynamic_casting
  virtual void setCollaborator(MultiMethodHandle<LinearSystemSolver> lss) = 0;

  /// Sets the ConvergenceMethod for this SpaceMethod to use
  /// @pre the pointer to LinearSystemSolver is not constant to
  ///      allow dynamic_casting
  virtual void setCollaborator(MultiMethodHandle<ConvergenceMethod> convMtd) = 0;

  /// Creates a MeshDataBuilder that suits this SpaceMethod
  /// This is only a build method. Destruction of the allocated object is
  /// the responsability of the client code.
  Common::SelfRegistPtr<MeshDataBuilder> createMeshBuilder();

  /// Get the Glocal Jacobian Sparsity
  /// Used by LinearSystemSolver's to allocate correctly the global matrix,
  /// possibly in parallel.
  /// @post Destruction of the allocated object is the responsability of the client code.
  Common::SelfRegistPtr<GlobalJacobianSparsity> createJacobianSparsity();

  /// Initialize the solution before starting the computation.
  /// @post pushs and pops the Namespace to which this Method belongs
  virtual void initializeSolution();

  /// Prepare to compute.
  /// Typically reset matrix and solutions to zero.
  /// @post pushs and pops the Namespace to which this Method belongs
  void prepareComputation();

  /// Compute the residual and the jacobian contributions
  /// coming from the space discretization
  /// @param factor is used to multiply the residual and jacobians
  ///               according to the different time stepping schemes
  /// @pre Each iteration must be called after prepareComputation()
  /// @post pushs and pops the Namespace to which this Method belongs
  void computeSpaceResidual(CFreal factor);

  /// Compute the contribution to the residual or jacobian
  /// coming from the time discretization
  /// @param factor is used to multiply the residual and jacobians
  ///               according to the different time stepping schemes
  /// @pre Each iteration must be called after prepareComputation()
  /// @post pushs and pops the Namespace to which this Method belongs
  void computeTimeResidual(CFreal factor);

  /// Apply the boundary conditions
  /// @post pushs and pops the Namespace to which this Method belongs
  void applyBC();

  /// PostProcess the solution.
  /// For instance for the application of a limiter or a filter.
  /// @post pushs and pops the Namespace to which this Method belongs
  void postProcessSolution();

  /// PreProcess the solution.
  /// For instance for the application of a limiter or a filter.
  /// @post pushs and pops the Namespace to which this Method belongs
  void preProcessSolution();

  /// Compute the diagonal block jacobian contributions
  /// coming from the space discretization
  /// @param factor is used to multiply the residual and jacobians
  ///               according to the different time stepping schemes
  /// @pre Each iteration must be called after prepareComputation()
  /// @post pushs and pops the Namespace to which this Method belongs
//   void computeSpaceDiagBlockJacobContrib(CFreal factor);

  /// Compute the diagonal block jacobian contributions
  /// coming from the time discretization
  /// @param factor is used to multiply the residual and jacobians
  ///               according to the different time stepping schemes
  /// @pre Each iteration must be called after prepareComputation()
  /// @post pushs and pops the Namespace to which this Method belongs
//   void computeTimeDiagBlockJacobContrib(CFreal factor);

  /// Compute the rhs for the given set of states
  /// coming from the space discretization
  /// @param factor is used to multiply the residual and jacobians
  ///               according to the different time stepping schemes
  /// @pre Each iteration must be called after prepareComputation()
  /// @post pushs and pops the Namespace to which this Method belongs
  void computeSpaceRhsForStatesSet(CFreal factor);

  /// Compute the rhs for the given set of states
  /// coming from the time discretization
  /// @param factor is used to multiply the residual and jacobians
  ///               according to the different time stepping schemes
  /// @pre Each iteration must be called after prepareComputation()
  /// @post pushs and pops the Namespace to which this Method belongs
  void computeTimeRhsForStatesSet(CFreal factor);

  /// Extrapolates the states to the node positions
  /// @post pushs and pops the Namespace to which this Method belongs
  void extrapolateStatesToNodes();

  /// Set the flag telling to compute the jacobians
  void setComputeJacobianFlag(bool flag);

  /// Set the flag telling to preprocess the solution only once
  void setOnlyPreprocessSolution(bool flag);
  
  /// Get the volume integrator of the space method.
  virtual Common::SafePtr<Framework::VolumeIntegrator> getVolumeIntegrator() {return CFNULL;}
  
  /// Action which is executed by the ActionLinstener for the "CF_ON_MESHADAPTER_BEFOREMESHUPDATE" Event
  /// @param eBefore the event which provoked this action
  /// @return an Event with a reply message in its body
  /// @post pushs and pops the Namespace to which this Method belongs
  Common::Signal::return_t beforeMeshUpdateAction( Common::Signal::arg_t eBefore);

  /// Action which is executed by the ActionLinstener for the "CF_ON_MESHADAPTER_AFTERMESHUPDATE" Event
  /// @param eAfter the event which provoked this action
  /// @return an Event with a reply message in its body
  /// @post pushs and pops the Namespace to which this Method belongs
  Common::Signal::return_t afterMeshUpdateAction( Common::Signal::arg_t eAfter);

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
    Common::DynamicFunctionCaller<SpaceMethod>::run_dynamic_function(func);
  }

  /// Gets the Class name
  static std::string getClassName() { return "SpaceMethod"; }

protected: // abstract interface implementations

  /// Action which is executed by the ActionLinstener for the "CF_ON_MAESTRO_MODIFYRESTART" Event
  /// @param eModifyRestart the event which provoked this action
  /// @return an Event with a message in its body
  Common::Signal::return_t modifyRestartAction( Common::Signal::arg_t eModifyRestart);

  /// Initialize the solution before starting the computation.
  /// This is the abstract function that the concrete methods must implement.
  virtual void initializeSolutionImpl(bool isRestart) = 0;

  /// Prepare to compute.
  /// Typically reset matrix and solutions to zero.
  /// This is the abstract function that the concrete methods must implement.
  virtual void prepareComputationImpl() = 0;

  /// Compute the residual and the jacobian contributions
  /// coming from the space discretization
  /// This is the abstract function that the concrete methods must implement.
  /// @param factor is used to multiply the residual and jacobians
  ///               according to the different time stepping schemes
  /// @pre Each iteration must be called after prepareComputation()
  virtual void computeSpaceResidualImpl(CFreal factor) = 0;

  /// Compute the contribution to the residual or jacobian
  /// coming from the time discretization
  /// This is the abstract function that the concrete methods must implement.
  /// @param factor is used to multiply the residual and jacobians
  ///               according to the different time stepping schemes
  /// @pre Each iteration must be called after prepareComputation()
  virtual void computeTimeResidualImpl(CFreal factor) = 0;

  /// Apply the boundary conditions
  /// This is the abstract function that the concrete methods must implement.
  virtual void applyBCImpl() = 0;

  /// Postprocess the solution.
  /// For instance for the application of a limiter or a filter..
  /// This is the abstract function that the concrete methods must implement.
  virtual void postProcessSolutionImpl() {}

  /// Preprocess the solution.
  /// For instance for the application of a limiter or a filter..
  /// This is the abstract function that the concrete methods must implement.
  virtual void preProcessSolutionImpl() {}
  
  /// Compute the diagonal block jacobian contributions
  /// coming from the space discretization
  /// This function should be overwritten by the concrete method.
  /// @param factor is used to multiply the residual and jacobians
  ///               according to the different time stepping schemes
  /// @pre Each iteration must be called after prepareComputation()
  /// @post pushs and pops the Namespace to which this Method belongs
//   virtual void computeSpaceDiagBlockJacobContribImpl(CFreal factor);

  /// Compute the diagonal block jacobian contributions
  /// coming from the time discretization
  /// This function should be overwritten by the concrete method.
  /// @param factor is used to multiply the residual and jacobians
  ///               according to the different time stepping schemes
  /// @pre Each iteration must be called after prepareComputation()
  /// @post pushs and pops the Namespace to which this Method belongs
//   virtual void computeTimeDiagBlockJacobContribImpl(CFreal factor);

  /// Compute the rhs for the given set of states
  /// coming from the space discretization
  /// This function should be overwritten by the concrete method.
  /// @param factor is used to multiply the residual and jacobians.
  ///               according to the different time stepping schemes
  /// @pre Each iteration must be called after prepareComputation()
  /// @post pushs and pops the Namespace to which this Method belongs
  virtual void computeSpaceRhsForStatesSetImpl(CFreal factor);

  /// Compute the rhs for the given set of states
  /// coming from the time discretization
  /// This function should be overwritten by the concrete method.
  /// @param factor is used to multiply the residual and jacobians.
  ///               according to the different time stepping schemes
  /// @pre Each iteration must be called after prepareComputation()
  /// @post pushs and pops the Namespace to which this Method belongs
  virtual void computeTimeRhsForStatesSetImpl(CFreal factor);

  /// Extrapolates the states to the node positions
  /// This is the abstract function that the concrete methods must implement.
  virtual void extrapolateStatesToNodesImpl() = 0;

  /// Action which is executed by the ActionLinstener for the "CF_ON_MESHADAPTER_BEFOREMESHUPDATE" Event
  /// This is the abstract function that the concrete methods must implement.
  /// @param eBefore the event which provoked this action
  /// @return an Event with a reply message in its body
  virtual Common::Signal::return_t beforeMeshUpdateActionImpl ( Common::Signal::arg_t eBefore) = 0;

  /// Action which is executed by the ActionLinstener for the "CF_ON_MESHADAPTER_AFTERMESHUPDATE" Event
  /// This is the abstract function that the concrete methods must implement.
  /// @param eAfter the event which provoked this action
  /// @return an Event with a reply message in its body
  virtual Common::Signal::return_t afterMeshUpdateActionImpl ( Common::Signal::arg_t eAfter) = 0;

protected: // functions

  /// Adds the ActionListener's of this EventListener to the EventHandler
  virtual void registActionListeners();

  /// Declares which functions can be called dynamically
  virtual void build_dynamic_functions();

protected: // member data

  /// name of the MeshDataBuilder used by the concrete SpaceMethod
  std::string m_builder;
  /// name of the GlobalJacobianSparsity computer used by the concrete SpaceMethod
  std::string m_sparsity;

  /// stored configuration arguments
  /// @todo this should be avoided and removed.
  ///       It is currently only a quick fix for the MeshBuilder
  Config::ConfigArgs m_stored_args;

private: // data

  /// do you want to restart from a previous solution
  bool m_restart;

}; // class SpaceMethod

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

Framework_Factory(SpaceMethod)

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_SpaceMethod_hh
