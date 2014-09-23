#ifndef COOLFluiD_Numerics_FluctSplit_FluctuationSplit_hh
#define COOLFluiD_Numerics_FluctSplit_FluctuationSplit_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/SpaceMethod.hh"
#include "FluctSplit/FluctuationSplitData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework  { class NumericalCommand; }

    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a Fluctuation Split spacial discretization method.
/// Fluctuation Splitting can also be called Residual Distribution.
/// @author Tiago Quintino
/// @author Andrea Lani
class FluctSplit_API FluctuationSplit : public Framework::SpaceMethod {
public:

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Constructor.
  /// @param name missing documentation
  explicit FluctuationSplit(const std::string& name);

  /// Default destructor.
  ~FluctuationSplit();

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

  /// Get the volume integrator of the space method.
  virtual Common::SafePtr<Framework::VolumeIntegrator> getVolumeIntegrator()  {  return _data->getVolumeIntegrator(); }

  /// Gets the Data aggregator of this method
  /// @return SafePtr to the MethodData
  virtual Common::SafePtr< Framework::MethodData > getMethodData () const;

  /// Gets the Data aggregator of this space method
  /// @return SafePtr to the SpaceMethodData
  virtual Common::SafePtr<Framework::SpaceMethodData> getSpaceMethodData();

  /// Data for the FluctuationSplit method
  Common::SafePtr<FluctuationSplitData> getData() const { return _data.getPtr(); }

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

  /// Clears (deletes) the stored dynamic data for Setup commands
  void clearSetupComs();

  /// Clears (deletes) the stored dynamic data for UnSetup commands
  void clearUnSetupComs();

  /// Configures the setup commands
  void configureSetupComs(Config::ConfigArgs& args);

  /// Configures the unsetup commands
  void configureUnsetupComs(Config::ConfigArgs& args);

private: // member data

  /// The Setup command to use
  std::vector<Common::SelfRegistPtr<FluctuationSplitCom> > _setups;

  /// The UnSetup command to use
  std::vector<Common::SelfRegistPtr<FluctuationSplitCom> > _unSetups;

  /// This command  extrapolates the states to the nodes
  Common::SelfRegistPtr<FluctuationSplitCom> _extrapolateStates;

  /// The command to use for computing the space rhs.
  Common::SelfRegistPtr<FluctuationSplitCom> _computeSpaceRHS;

  /// The command to use for computing the time rhs.
  Common::SelfRegistPtr<FluctuationSplitCom> _computeTimeRHS;

  /// The command to use for initializing the solution.
  std::vector<Common::SelfRegistPtr<FluctuationSplitCom> > _inits;

  /// The command to use for computing the voundary conditions.
  std::vector<Common::SelfRegistPtr<FluctuationSplitCom> > _bcs;

  /// The command to use for the action before the mesh is updated
  Common::SelfRegistPtr<FluctuationSplitCom> _beforeMeshUpdate;

  /// The command to use for the action after the mesh has been updated
  Common::SelfRegistPtr<FluctuationSplitCom> _afterMeshUpdate;

  /// The data to share between FluctuationSplitCom commands
  Common::SharedPtr<FluctuationSplitData> _data;

  /// The Setup string for configuration
  std::vector<std::string> _setupStr;

  /// The Setup Name string for configuration
  std::vector<std::string> _setupNameStr;

  /// The UnSetup string for configuration
  std::vector<std::string> _unSetupStr;

  /// The UnSetup Name string for configuration
  std::vector<std::string> _unSetupNameStr;

  /// The string for configuration of the _extrapolateStates command
  std::string _extrapolateStatesStr;

  /// The compute space RHS string for configuration
  std::string _computeSpaceRHSStr;

  /// The compute time RHS string for configuration
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

}; // class FluctuationSplit

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_FluctuationSplit_hh
