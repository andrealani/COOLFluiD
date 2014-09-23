#ifndef COOLFluiD_DiscontGalerkin_DiscontGalerkinSolver_hh
#define COOLFluiD_DiscontGalerkin_DiscontGalerkinSolver_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/SpaceMethod.hh"
#include "DiscontGalerkin/DiscontGalerkinSolverData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework { class NumericalCommand; }

  namespace DiscontGalerkin {

//////////////////////////////////////////////////////////////////////////////

/// This class implements an discontinuousGalerkinSolver
/// @author Martin Holik
/// @author Vit Dolejsi
class DiscontGalerkinSolver : public Framework::SpaceMethod {

public: // functions

  /// Constructor
  explicit DiscontGalerkinSolver(const std::string& name);

  /// Destructor
  ~DiscontGalerkinSolver();

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
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
  virtual Common::SafePtr< Framework::VolumeIntegrator > getVolumeIntegrator()  {  return CFNULL;  }

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

  /// Executed on "CF_ON_MESHADAPTER_BEFOREMESHUPDATE" event
  Common::Signal::return_t beforeMeshUpdateActionImpl(Common::Signal::arg_t eBefore);

  /// Executed on "CF_ON_MESHADAPTER_AFTERMESHUPDATE" event
  Common::Signal::return_t afterMeshUpdateActionImpl(Common::Signal::arg_t eAfter);

  /// Executed after solving linear system
  virtual void postProcessSolutionImpl();

  /// Clears (deletes) the solution initializing commands
  void clearInitComs();

  /// Clears (deletes) the stored dynamic data for BC commands
  void clearBCComs();

private: // data

  ///The Setup command to use
  Common::SelfRegistPtr< DiscontGalerkinSolverCom > m_setup;

  ///The UnSetup command to use
  Common::SelfRegistPtr< DiscontGalerkinSolverCom > m_unsetup;

  ///The solveCells command
  Common::SelfRegistPtr< DiscontGalerkinSolverCom > m_solveCells;

  ///The solveFaces command
  Common::SelfRegistPtr< DiscontGalerkinSolverCom > m_solveFaces;

  ///The set residual command
  Common::SelfRegistPtr< DiscontGalerkinSolverCom > m_setResidual;

  ///The stabilization command
  Common::SelfRegistPtr< DiscontGalerkinSolverCom > m_stabilization;

  ///The time Dependencies command
  // Common::SelfRegistPtr< DiscontGalerkinSolverCom > m_timeDependencies;

  ///The Setup string for configuration
  std::string m_setupStr;

  ///The UnSetup string for configuration
  std::string m_unsetupStr;

  /// The string for configuration of the solve cells command
  std::string m_solveCellsStr;

  /// The string for configuration of the set residual command
  std::string m_setResidualStr;

  /// The string for configuration of the solvefaces command
  std::string m_solveFacesStr;

  /// The string for configuration of the stabilization command
  std::string m_stabilizationStr;

  /// The string for configuration of the timeDependencies command
  //std::string m_timeDependentStr;

  ///The data to share between DiscontGalerkinSolverCom commands
  Common::SharedPtr< DiscontGalerkinSolverData > m_data;

  /// The command to use for initializing the solution.
  std::vector<Common::SelfRegistPtr<DiscontGalerkinSolverCom> > m_inits;

  /// The command to use for computing the voundary conditions.
  std::vector<Common::SelfRegistPtr<DiscontGalerkinSolverCom> > m_bcs;

  /// The solution initializing command types
  std::vector<std::string> m_initTypeStr;

  /// The solution initializing command names for configuration
  std::vector<std::string> m_initNameStr;

  /// The bc command types
  std::vector<std::string> m_bcTypeStr;

  /// The bc command names for configuration
  std::vector<std::string> m_bcNameStr;

//////////////////////////////////////////////////////////////////////////////

}; // class DiscontGalerkinSolver

//////////////////////////////////////////////////////////////////////////////

  }  // namespace DiscontGalerkin
}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_DiscontGalerkin_DiscontGalerkinSolver_hh

