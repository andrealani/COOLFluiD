#ifndef COOLFluiD_Muffin_Muffin_hh
#define COOLFluiD_Muffin_Muffin_hh

#include "Framework/SpaceMethod.hh"
#include "Muffin/MuffinData.hh"

namespace COOLFluiD {

  namespace Framework {  class NumericalCommand;  }

  namespace Muffin {


/// This is implementation of a steady-state RDS Navier-Stokes solver
/// @author Pedro Maciel
class Muffin : public Framework::SpaceMethod {

 public:

  /// Constructor
  explicit Muffin(const std::string& name);

  /// Destructor
  ~Muffin()
  {}

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Configures the method, by allocating its dynamic members
  void configure(Config::ConfigArgs& args);

  /// Gets the Data aggregator of this space method
  /// @return SafePtr to the SpaceMethodData
  inline Common::SafePtr< Framework::SpaceMethodData > getSpaceMethodData()
  {
    return m_data.getPtr();
  }

  /// Gets the Data aggregator of this space method
  /// @return SafePtr to the MethodData
  inline Common::SafePtr< Framework::MethodData > getMethodData() const
  {
    return m_data.getPtr();
  }

  /// Sets the LinearSystemSolver for this SpaceMethod to use
  /// @pre pointer to LSS is not constant to allow dynamic casting
  void setCollaborator(Framework::MultiMethodHandle< Framework::LinearSystemSolver > lss)
  {
    m_data->setLinearSystemSolver(lss);
  }

  /// Sets the ConvergenceMethod for this SpaceMethod to use
  void setCollaborator(Framework::MultiMethodHandle< Framework::ConvergenceMethod > convMtd)
  {
    // convergence method collaborator is made available to the commands
    // through the method data
    m_data->setConvergenceMethod(convMtd);
  }

  // Get the volume integrator of the space method.
  inline Common::SafePtr< Framework::VolumeIntegrator > getVolumeIntegrator()
  {
    return CFNULL;
  }


 protected:  // interface implementation functions

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

  // Compute the time contribution to residual
  void computeTimeResidualImpl(CFreal factor);

  /// Apply boundary conditions
  void applyBCImpl();

  /// Prepare to compute
  void prepareComputationImpl();

  /// Executed on "CF_ON_MESHADAPTER_BEFOREMESHUPDATE" event
  Common::Signal::return_t beforeMeshUpdateActionImpl(Common::Signal::arg_t eBefore);

  /// Executed on "CF_ON_MESHADAPTER_AFTERMESHUPDATE" event
  Common::Signal::return_t afterMeshUpdateActionImpl(Common::Signal::arg_t eAfter);


 private:  // data

  /// The data to share between commands
  Common::SharedPtr< MuffinData > m_data;

  /// Standard command (name) to link Systems to coupling and boundary conditions (configurable)
  std::string m_linkStr;

  /// Standard command (name) to build method required data (configurable)
  std::string m_setupStr;

  /// Standard command (name) to solve testcase (configurable)
  std::string m_solveStr;

  /// Standard command (name) to dispose of built data (configurable)
  std::string m_unsetupStr;

  /// Standard command to link Systems to coupling and boundary conditions
  Common::SelfRegistPtr< MuffinCom > m_link;

  /// Standard command to build method required data
  Common::SelfRegistPtr< MuffinCom > m_setup;

  /// Standard command to solve testcase
  Common::SelfRegistPtr< MuffinCom > m_solve;

  /// Standard command to dispose of built data
  Common::SelfRegistPtr< MuffinCom > m_unsetup;

  /// Special commands that can be called within loops (names/types, configurable)
  std::vector< std::string > m_specialcomds;

  /// Special commands that can be called within loops
  std::vector< Common::SelfRegistPtr< MuffinCom > > m_vcomm_specials;

  /// Loop commands (types, configurable)
  std::vector< std::string > m_loopcomds;

  /// System commands (types, configurable)
  std::vector< std::string > m_syscomds;

  /// Coupling condition commands (types, configurable)
  std::vector< std::string > m_cccomds;

  /// Boundary condition commands (types, configurable)
  std::vector< std::string > m_bccomds;

  /// Loop commands (names, configurable)
  std::vector< std::string > m_loopnames;

  /// System commands (names, configurable)
  std::vector< std::string > m_sysnames;

  /// Coupling condition commands (names, configurable)
  std::vector< std::string > m_ccnames;

  /// Boundary condition commands (names, configurable)
  std::vector< std::string > m_bcnames;

  /// Loop commands
  std::vector< Common::SelfRegistPtr< MuffinCom > > m_vcomm_loops;

  /// System commands
  std::vector< Common::SelfRegistPtr< MuffinCom > > m_vcomm_sys;

  /// Coupling condition commands
  std::vector< Common::SelfRegistPtr< MuffinCom > > m_vcomm_cc;

  /// Boundary condition commands
  std::vector< Common::SelfRegistPtr< MuffinCom > > m_vcomm_bc;

}; // class Muffin


  }  // namespace Muffin
}  // namespace COOLFluiD

#endif // COOLFluiD_Muffin_Muffin_hh

