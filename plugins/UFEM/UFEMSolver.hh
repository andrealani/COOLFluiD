#ifndef COOLFluiD_UFEMSolver_UFEMSolver_hh
#define COOLFluiD_UFEMSolver_UFEMSolver_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/SpaceMethod.hh"
#include "UFEM/UFEMSolverData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework { class NumericalCommand; }

  namespace UFEM {

//////////////////////////////////////////////////////////////////////////////

/// This class implements an Unsteady FEM solver
/// @author Tiago Quintino
class UFEM_API UFEMSolver : public Framework::SpaceMethod {

public: // functions

  /// Constructor
  explicit UFEMSolver(const std::string& name);

  /// Destructor
  ~UFEMSolver();

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Configures the method, by allocating its dynamic members
  virtual void configure ( Config::ConfigArgs& args );

  /// Gets the Data aggregator of this space method
  /// @return SafePtr to the SpaceMethodData
  virtual Common::SafePtr< Framework::SpaceMethodData > getSpaceMethodData() {  return m_data.getPtr(); }

  /// Gets the Data aggregator of this space method
  /// @return SafePtr to the MethodData
  virtual Common::SafePtr< Framework::MethodData > getMethodData () const {  return m_data.getPtr(); }

  /// Sets the LinearSystemSolver for this SpaceMethod to use
  /// @pre pointer to LSS is not constant to allow dynamic casting
  void setCollaborator( Framework::MultiMethodHandle< Framework::LinearSystemSolver > lss );

  /// Sets the ConvergenceMethod for this SpaceMethod to use
  /// @pre the pointer to ConvergenceMethod is not constant to
  ///      allow dynamic_casting
  void setCollaborator(Framework::MultiMethodHandle<Framework::ConvergenceMethod> convMtd);


  /// Gets the volume integrator of the space method.
  virtual Common::SafePtr< Framework::VolumeIntegrator > getVolumeIntegrator() { return CFNULL; }

protected: // interface implementation functions

  /// Sets up data, commands and strategies of this Method
  virtual void setMethodImpl();

  /// Unsets the data, commands and strategies of this Method
  virtual void unsetMethodImpl();

  /// Extrapolates the states to the node positions
  void extrapolateStatesToNodesImpl();

  /// Initialize the solution before starting the computation
  void initializeSolutionImpl(bool isRestart);

  /// Set matrix, right hand side and solve system
  void computeSpaceResidualImpl(CFreal factor);

  /// Compute the time contribution to residual
  void computeTimeResidualImpl(CFreal factor);

  /// Apply boundary conditions
  void applyBCImpl();

  /// Prepare to compute
  void prepareComputationImpl();

  /// Executed on "CF_ON_MESHADAPTER_BEFOREMESHUPDATE" event
  Common::Signal::return_t beforeMeshUpdateActionImpl(Common::Signal::arg_t eBefore);

  /// Executed on "CF_ON_MESHADAPTER_AFTERMESHUPDATE" event
  Common::Signal::return_t afterMeshUpdateActionImpl(Common::Signal::arg_t eAfter);

  /// Clears the solution initializing commands
  void clearInitComs();
  /// Clears the boundary condition commands
  void clearBCComs();

private: // data

  ///The Setup command(s) to use
  std::vector< Common::SelfRegistPtr< UFEMSolverCom > > m_setup;

  ///The UnSetup command to use
  Common::SelfRegistPtr< UFEMSolverCom > m_unsetup;

  ///The solve command
  Common::SelfRegistPtr< UFEMSolverCom > m_solve;

  ///The clean command
  Common::SelfRegistPtr< UFEMSolverCom > m_clean;

  ///The command to use for initializing the solution.
  std::vector< Common::SelfRegistPtr<UFEMSolverCom> > m_inits;

  ///The command to use for computing the boundary conditions.
  std::vector<Common::SelfRegistPtr<UFEMSolverCom> > m_bcs;

  ///The Setup string for configuration
  std::vector< std::string > m_setupStr;

  ///The UnSetup string for configuration
  std::string m_unsetupStr;

  /// The string for configuration of the _prepare command
  std::string m_solveStr;

  /// The string for configuration of the clean command
  std::string m_cleanStr;

  ///The solution initializing command types
  std::vector<std::string> m_initTypeStr;

  ///The solution initializing command names for configuration
  std::vector<std::string> m_initNameStr;

  ///The bc command types
  std::vector<std::string> m_bcTypeStr;

  ///The bc command names for configuration
  std::vector<std::string> m_bcNameStr;

  ///The data to share between UFEMSolverCom commands
  Common::SharedPtr< UFEMSolverData > m_data;

//////////////////////////////////////////////////////////////////////////////

}; // class UFEMSolver

//////////////////////////////////////////////////////////////////////////////

  }  // namespace UFEM
}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_UFEMSolver_UFEMSolver_hh

