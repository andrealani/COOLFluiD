// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_FluxReconstructionMethod_FluxReconstructionSolver_hh
#define COOLFluiD_FluxReconstructionMethod_FluxReconstructionSolver_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/SpaceMethod.hh"
#include "FluxReconstructionMethod/FluxReconstructionSolverData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework { class NumericalCommand; }

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/// This class implements a FluxReconstruction solver
/// @author Alexander Papen
/// @author Ray Vandenhoeck
class FluxReconstructionSolver : public Framework::SpaceMethod {

public: // functions

  /// Constructor
  explicit FluxReconstructionSolver(const std::string& name);

  /// Destructor
  ~FluxReconstructionSolver();

  /// Defines the Config Option's of this class
  /// @param options an OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Configures the method, by allocating its dynamic members
  virtual void configure ( Config::ConfigArgs& args );

  /// Gets the Data aggregator of this space method
  /// @return SafePtr to the SpaceMethodData
  virtual Common::SafePtr< Framework::SpaceMethodData > getSpaceMethodData()
  {
    return m_data.getPtr();
  }

  /// Gets the Data aggregator of this space method
  /// @return SafePtr to the MethodData
  virtual Common::SafePtr< Framework::MethodData > getMethodData () const
  {
    return m_data.getPtr();
  }

  /// Sets the LinearSystemSolver for this SpaceMethod to use
  /// @pre pointer to LSS is not constant to allow dynamic casting
  void setCollaborator( Framework::MultiMethodHandle< Framework::LinearSystemSolver > lss );

  /// Sets the ConvergenceMethod for this SpaceMethod to use
  /// @pre the pointer to ConvergenceMethod is not constant to
  ///      allow dynamic_casting
  void setCollaborator(Framework::MultiMethodHandle<Framework::ConvergenceMethod> convMtd);


  /// Gets the volume integrator of the space method.
  virtual Common::SafePtr< Framework::VolumeIntegrator > getVolumeIntegrator()
  {
    return CFNULL;
  }

protected: // interface implementation functions

  /// Sets up data, commands and strategies of this Method
  virtual void setMethodImpl();

  /// Unsets the data, commands and strategies of this Method
  virtual void unsetMethodImpl();

  /// Extrapolates the states to the node positions
  virtual void extrapolateStatesToNodesImpl();

  /// Initialize the solution before starting the computation
  virtual void initializeSolutionImpl(bool isRestart);

  /// Set matrix, right hand side and solve system
  virtual void computeSpaceResidualImpl(CFreal factor);

  /// Compute the time contribution to residual
  virtual void computeTimeResidualImpl(CFreal factor);

  /// Apply boundary conditions
  virtual void applyBCImpl();

  /// Prepare to compute
  virtual void prepareComputationImpl();

  /// Postprocess the solution.
  virtual void postProcessSolutionImpl() {};

  /// Executed on "CF_ON_MESHADAPTER_BEFOREMESHUPDATE" event
  Common::Signal::return_t beforeMeshUpdateActionImpl(Common::Signal::arg_t eBefore);

  /// Executed on "CF_ON_MESHADAPTER_AFTERMESHUPDATE" event
  Common::Signal::return_t afterMeshUpdateActionImpl(Common::Signal::arg_t eAfter);

private: // data

  ///The Setup command to use
  Common::SelfRegistPtr< FluxReconstructionSolverCom > m_setup;

  ///The UnSetup command to use
  Common::SelfRegistPtr< FluxReconstructionSolverCom > m_unsetup;

  ///The solve command
  Common::SelfRegistPtr< FluxReconstructionSolverCom > m_solve;

  ///The Setup string for configuration
  std::string m_setupStr;

  ///The UnSetup string for configuration
  std::string m_unsetupStr;

  /// The string for configuration of the _prepare command
  std::string m_solveStr;

  ///The data to share between FluxReconstructionSolverCom commands
  Common::SharedPtr< FluxReconstructionSolverData > m_data;

//////////////////////////////////////////////////////////////////////////////

}; // class FluxReconstructionSolver

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionMethod_FluxReconstructionSolver_hh

