#ifndef COOLFluiD_Numerics_FiniteElement_FiniteElementMethod_hh
#define COOLFluiD_Numerics_FiniteElement_FiniteElementMethod_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/SpaceMethod.hh"
#include "FiniteElementMethodData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework { class NumericalCommand; }

  namespace Numerics {

    namespace FiniteElement {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a Finite Element spacial discretization method.
/// @author Tiago Quintino
/// @author Thomas Wuilbaut
class FiniteElementMethod : public Framework::SpaceMethod {

public: // functions

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Constructor.
  /// @param name missing documentation
  explicit FiniteElementMethod(const std::string& name);


  /// Default destructor.
  ~FiniteElementMethod();

  /// Configures the method, by allocating its dynamic members.
  /// @param args missing documentation
  virtual void configure ( Config::ConfigArgs& args );

  /// Sets the LinearSystemSolver for this SpaceMethod to use
  /// @pre the pointer to LinearSystemSolver is not constant to
  ///      allow dynamic_casting
  void setCollaborator(Framework::MultiMethodHandle<Framework::LinearSystemSolver> lss);

  /// Sets the ConvergenceMethod for this SpaceMethod to use
  /// @pre the pointer to ConvergenceMethod is not constant to
  ///      allow dynamic_casting
  void setCollaborator(Framework::MultiMethodHandle<Framework::ConvergenceMethod> convMtd);

  /// Gets a vector with all the NumericalStrategy's this method will use.
  /// @return vector with the strategy pointers.
  virtual std::vector<Common::SafePtr<Framework::NumericalStrategy> > getStrategyList() const;

  /// Get the volume integrator of the space method.
  virtual Common::SafePtr<Framework::VolumeIntegrator> getVolumeIntegrator() { return m_data->getVolumeIntegrator(); }

  /// Gets the Data aggregator of this method
  /// @return SafePtr to the MethodData
  virtual Common::SafePtr< Framework::MethodData > getMethodData () const;

  /// Gets the Data aggregator of this space method
  /// @return SafePtr to the SpaceMethodData
  virtual Common::SafePtr<Framework::SpaceMethodData> getSpaceMethodData();

protected: // interface implementation functions

  /// Sets up the data, commands and strategies of this Method
  /// @see Method::setMethod()
  virtual void setMethodImpl();

  /// Unsets the data, commands and strategies of this Method
  /// @see Method::unsetMethod()
  virtual void unsetMethodImpl();

  /// Extrapolates the states to the node positions
  void extrapolateStatesToNodesImpl();

  /// Initialize the solution before starting the computation.
  void initializeSolutionImpl(bool isRestart);

  /// Compute the space contribution to the residual,
  /// sets the sysstem, right hand side and solves it.
  void computeSpaceResidualImpl(CFreal factor);

  /// Compute the time contribution to the residual.
  void computeTimeResidualImpl(CFreal factor);

  /// Apply the boundary conditions
  void applyBCImpl();

  /// Prepare to compute.
  void prepareComputationImpl();

  /// Action which is executed by the ActionLinstener for the "CF_ON_MESHADAPTER_BEFOREMESHUPDATE" Event
  /// @param eBefore the event which provoked this action
  /// @return an Event with a message in its body
  Common::Signal::return_t beforeMeshUpdateActionImpl(Common::Signal::arg_t eBefore);

  /// Action which is executed by the ActionLinstener for the "CF_ON_MESHADAPTER_AFTERMESHUPDATE" Event
  /// @param eAfter the event which provoked this action
  /// @return an Event with a message in its body
  Common::Signal::return_t afterMeshUpdateActionImpl(Common::Signal::arg_t eAfter);

protected: // helper functions

  /// Clears (deletes) the solution initializing commands
  void clearInitComs();

  /// Clears (deletes) the stored dynamic data for BC commands
  void clearBCComs();

  /// Checks if the system matrix shouldbe frozen
  void checkMatrixFrozen() const;

private: // data

  ///The Setup command to use
  Common::SelfRegistPtr<FiniteElementMethodCom> m_setup;

  ///The UnSetup command to use
  Common::SelfRegistPtr<FiniteElementMethodCom> m_unSetup;

  ///The Prepare command, used to prepare the computation
  Common::SelfRegistPtr<FiniteElementMethodCom> m_prepare;

  ///This command sets the nodal states
  Common::SelfRegistPtr<FiniteElementMethodCom> m_extrapolateStates;

  ///The command to use for computing the spacial contribution to the residual
  Common::SelfRegistPtr<FiniteElementMethodCom> m_computeSpaceResidual;

  ///The command to use for computing the time contribution to the residual
  Common::SelfRegistPtr<FiniteElementMethodCom> m_computeTimeResidual;

  ///The command to use for initializing the solution.
  std::vector<Common::SelfRegistPtr<FiniteElementMethodCom> > m_inits;

  ///The command to use for computing the boundary conditions.
  std::vector<Common::SelfRegistPtr<FiniteElementMethodCom> > m_bcs;

  ///The Setup string for configuration
  std::string m_setupStr;

  ///The UnSetup string for configuration
  std::string m_unSetupStr;

  /// The string for configuration of the m_prepare command
  std::string m_prepareStr;

  /// The string for configuration of the m_extrapolateStates command
  std::string m_extrapolateStatesStr;

  ///The ComputeSpaceResidual string for configuration
  std::string m_computeSpaceResidualStr;

  ///The ComputeTimeResidual string for configuration
  std::string m_computeTimeResidualStr;

  ///The solution initializing command types
  std::vector<std::string> m_initTypeStr;

  ///The solution initializing command names for configuration
  std::vector<std::string> m_initNameStr;

  ///The bc command types
  std::vector<std::string> m_bcTypeStr;

  ///The bc command names for configuration
  std::vector<std::string> m_bcNameStr;

  ///The data to share between FiniteElementMethodCom commands
  Common::SharedPtr<FiniteElementMethodData> m_data;

}; // class FiniteElementMethod

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteElement

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteElement_FiniteElementMethod_hh
