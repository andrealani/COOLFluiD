#ifndef COOLFluiD_Numerics_FiniteVolume_CellCenterFVM_hh
#define COOLFluiD_Numerics_FiniteVolume_CellCenterFVM_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/SpaceMethod.hh"
#include "CellCenterFVMData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
    class NumericalCommand;
  }

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a Finite Volume spatial discretization method.
/// @author Andrea Lani
/// @author Tiago Quintino
class CellCenterFVM : public Framework::SpaceMethod {
public:

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Constructor.
  /// @param name missing documentation
  explicit CellCenterFVM(const std::string& name);


  /// Default destructor.
  ~CellCenterFVM();

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
  virtual std::vector<Common::SafePtr<Framework::NumericalStrategy> > getStrategyList () const;
  
  /// Data for the CellCenterFVM method
  Common::SafePtr<CellCenterFVMData> getData() const
  {
    return _data.getPtr();
  }

  /// Gets the Data aggregator of this space method
  /// @return SafePtr to the SpaceMethodData
  virtual Common::SafePtr<Framework::SpaceMethodData> getSpaceMethodData();

  /// Gets the Data aggregator of this method
  /// @return SafePtr to the MethodData
  virtual Common::SafePtr< Framework::MethodData > getMethodData () const;

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

  /// Compute the residual and the jacobian contributions
  /// coming from the space discretization
  void computeSpaceResidualImpl(CFreal factor);

  /// Compute the residual coming from the time discretization
  void computeTimeResidualImpl(CFreal factor);

  /// Compute the rhs for the given set of states
  /// coming from the space discretization
  /// @param factor is used to multiply the residual and jacobians
  ///               according to the different time stepping schemes
  /// @pre Each iteration must be called after prepareComputation()
  /// @post pushs and pops the Namespace to which this Method belongs
  void computeSpaceRhsForStatesSetImpl(CFreal factor);

  /// Compute the rhs for the given set of states
  /// coming from the time discretization
  /// @param factor is used to multiply the residual and jacobians
  ///               according to the different time stepping schemes
  /// @pre Each iteration must be called after prepareComputation()
  /// @post pushs and pops the Namespace to which this Method belongs
  void computeTimeRhsForStatesSetImpl(CFreal factor);

  /// Apply the boundary conditions
  void applyBCImpl();

  /// Prepare to compute.
  void prepareComputationImpl();

  /// Postprocess the solution.
  void postProcessSolutionImpl();

  /// Preprocess the solution.
  /// For instance for the application of a limiter or a filter..
  /// This is the abstract function that the concrete methods must implement.
  void preProcessSolutionImpl();
  
  /// Action which is executed by the ActionLinstener for the "CF_ON_MESHADAPTER_BEFOREMESHUPDATE" Event
  /// @param eBefore the event which provoked this action
  /// @return an Event with a message in its body
  Common::Signal::return_t beforeMeshUpdateActionImpl(Common::Signal::arg_t eBefore);

  /// Action which is executed by the ActionLinstener for the "CF_ON_MESHADAPTER_AFTERMESHUPDATE" Event
  /// @param eAfter the event which provoked this action
  /// @return an Event with a message in its body
  Common::Signal::return_t afterMeshUpdateActionImpl(Common::Signal::arg_t eAfter);

protected: //  helper functions

  /// Clears (deletes) the solution initializing commands
  void clearInitComs();

  /// Clears (deletes) the stored dynamic data for BC commands
  void clearBCComs();

  /// Clears (deletes) the stored dynamic data for Setup commands
  void clearSetupComs();

  /// Clears (deletes) the stored dynamic data for UnSetup commands
  void clearUnSetupComs();

  /// Clears (deletes) the stored dynamic data for Pre-process commands
  void clearPreProcessComs();
  
  /// Checks if the system matrix shouldbe frozen
  void checkMatrixFrozen() const;

private: // member data

  /// The Setup command to use
  std::vector<Common::SelfRegistPtr<CellCenterFVMCom> > _setups;

  /// The UnSetup command to use
  std::vector<Common::SelfRegistPtr<CellCenterFVMCom> > _unSetups;

  /// The pre-process command to use
  std::vector<Common::SelfRegistPtr<CellCenterFVMCom> > _preProcess;

  /// This command sets the nodal states
  Common::SelfRegistPtr<CellCenterFVMCom> _extrapolateStates;

  /// The command to use for computing the space rhs.
  Common::SelfRegistPtr<CellCenterFVMCom> _computeSpaceRHS;

  /// The command to use for computing the time rhs.
  Common::SelfRegistPtr<CellCenterFVMCom> _computeTimeRHS;

  /// The command to use for initializing the solution.
  std::vector<Common::SelfRegistPtr<CellCenterFVMCom> > _inits;

  /// The command to use for computing the voundary conditions.
  std::vector<Common::SelfRegistPtr<CellCenterFVMCom> > _bcs;

  /// The command to use for the action before the mesh is updated
  Common::SelfRegistPtr<CellCenterFVMCom> _beforeMeshUpdate;

  /// The command to use for the action after the mesh has been updated
  Common::SelfRegistPtr<CellCenterFVMCom> _afterMeshUpdate;

  /// The command that computes the contribution of the spatial discretization
  /// to the rhs of a given set of states (== a cell), for use with LU-SGS.
  Common::SelfRegistPtr< CellCenterFVMCom > _spaceRHSForGivenCell;

  /// The command that computes the contribution of the temporal discretization
  /// to the rhs of a given set of states (== a cell), for use with LU-SGS.
  Common::SelfRegistPtr< CellCenterFVMCom > _timeRHSForGivenCell;
  
  /// The data to share between CellCenterFVMCom commands
  Common::SharedPtr<CellCenterFVMData> _data;
  
  /// Flag telling if BC's have been applied already
  bool _isBcApplied;
  
  /// Flag telling if to use only the initialization commands for initialization
  bool _useOnlyInitComs;

  /// The Setup string for configuration
  std::vector<std::string> _setupStr;

  /// The Setup Name string for configuration
  std::vector<std::string> _setupNameStr;

  /// The UnSetup string for configuration
  std::vector<std::string> _unSetupStr;

  /// The UnSetup Name string for configuration
  std::vector<std::string> _unSetupNameStr;
  
  /// The Pre-process string for configuration
  std::vector<std::string> _preProcessStr;
  
  /// The Pre-process Name string for configuration
  std::vector<std::string> _preProcessNameStr;
  
  /// The string for configuration of the _extrapolateStates command
  std::string _extrapolateStatesStr;
  
  /// The CompRHS string for configuration
  std::string _computeSpaceRHSStr;

  /// The CompRHS string for configuration
  std::string _computeTimeRHSStr;

  /// The solution initializing command types
  std::vector<std::string> _initTypeStr;

  /// The solution initializing command names for configuration
  std::vector<std::string> _initNameStr;

  /// The bc command types
  std::vector<std::string> _bcTypeStr;

  /// The bc command names for configuration
  std::vector<std::string> _bcNameStr;

  /// The string for the configuration of the action before mesh update
  std::string _beforeMeshUpdateStr;

  /// The string for the configuration of the action after mesh update
  std::string _afterMeshUpdateStr;

  /// The string for configuration of the m_spaceRHSForGivenCell command
  std::string _spaceRHSForGivenCellStr;

  /// The string for configuration of the m_timeRHSForGivenCell command
  std::string _timeRHSForGivenCellStr;

}; // class CellCenterFVM

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_CellCenterFVM_hh
