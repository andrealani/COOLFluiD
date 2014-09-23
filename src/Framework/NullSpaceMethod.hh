// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_NullSpaceMethod_hh
#define COOLFluiD_Framework_NullSpaceMethod_hh

//////////////////////////////////////////////////////////////////////////////

#include "SpaceMethod.hh"
#include "SpaceMethodData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

    class LinearSystemSolver;

//////////////////////////////////////////////////////////////////////////////

/// This class represents a NullSpaceMethod.
/// @author Tiago Quintino
class Framework_API NullSpaceMethod : public SpaceMethod {
public:

  /// Constructor.
  explicit NullSpaceMethod(const std::string& name);

  /// Destructor
  ~NullSpaceMethod();

  /// Sets the LinearSystemSolver for this SpaceMethod to use
  /// @pre the pointer to LinearSystemSolver is not constant to
  ///      allow dynamic_casting
  void setCollaborator(MultiMethodHandle<LinearSystemSolver> lss);

  /// Sets the ConvergenceMethod for this SpaceMethod to use
  /// @pre the pointer to ConvergenceMethod is not constant to
  ///      allow dynamic_casting
  void setCollaborator(Framework::MultiMethodHandle<Framework::ConvergenceMethod> convMtd);

  /// Checks if this object is a Null object.
  /// @return true since this is NullSpaceMethod
  bool isNull() const {  return true; }

  /// Get the volume integrator of the space method.
  virtual Common::SafePtr<Framework::VolumeIntegrator> getVolumeIntegrator() {   return CFNULL;  }

  /// Gets the Data aggregator of this method
  /// @return SafePtr to the MethodData
  virtual Common::SafePtr< Framework::MethodData > getMethodData () const;

  /// Gets the Data aggregator of this space method
  /// @return SafePtr to the SpaceMethodData
  virtual Common::SafePtr<Framework::SpaceMethodData> getSpaceMethodData();
  
  /// Initialize the solution before starting the computation.
  /// @post pushs and pops the Namespace to which this Method belongs
  virtual void initializeSolution();
  
protected: // interface implementation functions

  /// Sets up the data, commands and strategies of this Method
  /// @see Method::setMethod()
  virtual void setMethodImpl();

  /// Unsets the data, commands and strategies of this Method
  /// @see Method::unsetMethod()
  virtual void unsetMethodImpl();

  /// Initialize the solution before starting the computation.
  /// @see SpaceMethod::initializeSolution()
  void initializeSolutionImpl(bool isRestart);

  /// Extrapolates the states to the node positions
  void extrapolateStatesToNodesImpl();

  /// Compute the residual and the jacobian contributions
  /// coming from the space discretization
  void computeSpaceResidualImpl(CFreal factor);

  /// Compute the residual coming from the time discretization
  void computeTimeResidualImpl(CFreal factor);

  /// Apply the boundary conditions
  void applyBCImpl();

  /// Prepare to compute.
  void prepareComputationImpl();

  /// Postprocess the solution.
  void postProcessSolutionImpl() {};

  /// Action which is executed by the ActionLinstener for the "CF_ON_MESHADAPTER_BEFOREMESHUPDATE" Event
  /// @param eBefore the event which provoked this action
  /// @return an Event with a message in its body
  Common::Signal::return_t beforeMeshUpdateActionImpl(Common::Signal::arg_t eBefore);

  /// Action which is executed by the ActionLinstener for the "CF_ON_MESHADAPTER_AFTERMESHUPDATE" Event
  /// @param eAfter the event which provoked this action
  /// @return an Event with a message in its body
  Common::Signal::return_t afterMeshUpdateActionImpl(Common::Signal::arg_t eAfter);

}; // class NullSpaceMethod

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_NullSpaceMethod_hh
