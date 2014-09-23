#ifndef COOLFluiD_SpectralFD_SpectralFDMethod_hh
#define COOLFluiD_SpectralFD_SpectralFDMethod_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/SpaceMethod.hh"
#include "SpectralFD/SpectralFDMethodData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
    class NumericalCommand;
  }

  namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

    // forward declarations
    class BndFaceTermRHSSpectralFD;
    class ConvBndFaceTermRHSSpectralFD;
    class DiffBndFaceTermRHSSpectralFD;

//////////////////////////////////////////////////////////////////////////////

/// This class implements a SpectralFD solver
/// @author Kris Van Den Abeele
class SpectralFDMethod : public Framework::SpaceMethod {

public: // functions

  /// Constructor
  explicit SpectralFDMethod(const std::string& name);

  /// Destructor
  ~SpectralFDMethod() {}

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

  /// Gets the Data aggregator of this space method
  /// @return SafePtr to the MethodData
  virtual Common::SafePtr< SpectralFDMethodData > getData()
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

  /// Get the volume integrator of the space method.
  virtual Common::SafePtr< Framework::VolumeIntegrator > getVolumeIntegrator()  { return CFNULL; }

  /// Defined the strategy list of this Method
  /// @todo Remove this function and configure strategies using the configureStrategy template function
  ///       provided by MethodData
  std::vector< Common::SafePtr< Framework::NumericalStrategy > > getStrategyList() const;

protected: // interface implementation functions

  /// Sets up data, commands and strategies of this Method
  virtual void setMethodImpl()
  {
    SpaceMethod::setMethodImpl();

    setupCommandsAndStrategies();
    cf_assert(m_setup.isNotNull());
    m_setup->execute();
  }

  /// Unsets the data, commands and strategies of this Method
  virtual void unsetMethodImpl()
  {
    cf_assert(m_unsetup.isNotNull());
    m_unsetup->execute();
    unsetupCommandsAndStrategies();

    SpaceMethod::unsetMethodImpl();
  }

  /// Extrapolates the states to the node positions
  void extrapolateStatesToNodesImpl();

  /// Initialize the solution before starting the computation.
  /// @see SpaceMethod::initializeSolution()
  void initializeSolutionImpl(bool isRestart);

  /// Set matrix, right hand side and solve system
  void computeSpaceResidualImpl(CFreal factor);

  /// Compute the time contribution to residual
  void computeTimeResidualImpl(CFreal factor);

  /// Apply boundary conditions
  void applyBCImpl();

  /// Apply boundary conditions for diffusive terms
  void applyBCDiffImpl();

  /// add source terms
  void addSourceTermsImpl();

  /// Prepare to compute.
  /// This should clear the rhs and the update coefficients
  void prepareComputationImpl();

  /// Postprocess the solution.
  void postProcessSolutionImpl();

  /// Compute the rhs for the given set of states.
  void computeSpaceRhsForStatesSetImpl(CFreal factor);

  /// Compute the rhs for the given set of states.
  void computeTimeRhsForStatesSetImpl(CFreal factor);

  /// Executed on "CF_ON_MESHADAPTER_BEFOREMESHUPDATE" event
  Common::Signal::return_t beforeMeshUpdateActionImpl(Common::Signal::arg_t eBefore)
  {
    CFAUTOTRACE;
    return Common::Signal::return_t ();
  }

  /// Executed on "CF_ON_MESHADAPTER_AFTERMESHUPDATE" event
  Common::Signal::return_t afterMeshUpdateActionImpl(Common::Signal::arg_t eAfter)
  {
    CFAUTOTRACE;
    return Common::Signal::return_t ();
  }

private: // functions

  void configureSourceTermCommands( Config::ConfigArgs& args );

  void configureInitCommands( Config::ConfigArgs& args );

  void configureBcCommands( Config::ConfigArgs& args );

private: // data

  ///The Setup command to use
  Common::SelfRegistPtr< SpectralFDMethodCom > m_setup;

  ///The UnSetup command to use
  Common::SelfRegistPtr< SpectralFDMethodCom > m_unsetup;

  /// The extrapolate command
  Common::SelfRegistPtr< SpectralFDMethodCom > m_extrapolate;

  /// Command used to prepare the computation
  Common::SelfRegistPtr< SpectralFDMethodCom > m_prepare;

  /// Command used to limit a solution
  Common::SelfRegistPtr< SpectralFDMethodCom > m_limiter;

  /// The command that computes the volume terms of the discretization of the convective terms
  Common::SelfRegistPtr< SpectralFDMethodCom > m_convVolTerm;

  /// The command that computes the face terms of the discretization of the convective terms
  Common::SelfRegistPtr< SpectralFDMethodCom > m_convFaceTerm;

  /// The command that computes the volume terms of the discretization of the diffusive terms
  Common::SelfRegistPtr< SpectralFDMethodCom > m_diffVolTerm;

  /// The command that computes the face terms of the discretization of the diffusive terms
  Common::SelfRegistPtr< SpectralFDMethodCom > m_diffFaceTerm;

  /// The command that divides rhs and updateCoeffs by the cell volumes
  Common::SelfRegistPtr< SpectralFDMethodCom > m_divideRHSByCellVol;

  /// The command that computes the contribution of the time discretization to the rhs
  Common::SelfRegistPtr< SpectralFDMethodCom > m_timeRHSJacob;

  /// The command that computes the contribution of the spatial discretization
  /// to the rhs of a given set of states (== a cell), for use with LU-SGS.
  Common::SelfRegistPtr< SpectralFDMethodCom > m_spaceRHSForGivenCell;

  /// The command that computes the contribution of the temporal discretization
  /// to the rhs of a given set of states (== a cell), for use with LU-SGS.
  Common::SelfRegistPtr< SpectralFDMethodCom > m_timeRHSForGivenCell;

  /// The Setup string for configuration
  std::string m_setupStr;

  /// The UnSetup string for configuration
  std::string m_unsetupStr;

  /// The string for configuration of the extrapolate command
  std::string m_extrapolateStr;

  /// The string for configuration of the m_prepare command
  std::string m_prepareStr;

  /// The string for configuration of the m_limiter command
  std::string m_limiterStr;

  /// The string for configuration of the commands that compute
  /// RHS (and optionally the Jacobian)
  std::string m_spaceRHSJacobStr;

  /// The string for configuration of the command that divides the RHS and updateCoeffs by the volumes
  std::string m_divideRHSByCellVolStr;

  /// The string for configuration of the m_timeRHSJacob command
  std::string m_timeRHSJacobStr;

  /// The string for configuration of the m_spaceRHSForGivenCell command
  std::string m_spaceRHSForGivenCellStr;

  /// The string for configuration of the m_timeRHSForGivenCell command
  std::string m_timeRHSForGivenCellStr;

  /// The commands for the source terms
  std::vector< Common::SelfRegistPtr< SpectralFDMethodCom > > m_srcTerms;

  /// The source term command types
  std::vector<std::string> m_srcTermTypeStr;

  /// The source term command names for configuration
  std::vector<std::string> m_srcTermNameStr;

  /// The commands to use for initializing the solution.
  std::vector< Common::SelfRegistPtr< SpectralFDMethodCom > > m_inits;

  /// The solution initializing command types
  std::vector<std::string> m_initTypeStr;

  /// The solution initializing command names for configuration
  std::vector<std::string> m_initNameStr;

  /// The commands to use for applying the boundary conditions for the convective terms,
  /// with ConvBndFaceTermRHSSpectralFD as type
  std::vector< Common::SafePtr< ConvBndFaceTermRHSSpectralFD > > m_bcs;

  /// The commands to use for applying the boundary conditions for the convective terms
  std::vector< Common::SelfRegistPtr< SpectralFDMethodCom > > m_bcsComs;

  /// The commands to use for applying the boundary conditions for the diffusive terms,
  /// with DiffBndFaceTermRHSSpectralFD as type
  std::vector< Common::SafePtr< DiffBndFaceTermRHSSpectralFD > > m_bcsDiff;

  /// The commands to use for applying the boundary conditions for the diffusive terms
  std::vector< Common::SelfRegistPtr< SpectralFDMethodCom > > m_bcsDiffComs;

  /// The commands to use for applying the boundary conditions for both convective and diffusive terms,
  /// with BndFaceTermRHSSpectralFD as type
  std::vector< Common::SafePtr< BndFaceTermRHSSpectralFD > > m_bcsConvDiff;

  /// The boundary condition command names for configuration
  std::vector<std::string> m_bcNameStr;

  ///The data to share between SpectralFDMethodCom commands
  Common::SharedPtr< SpectralFDMethodData > m_data;

}; // class SpectralFDMethod

//////////////////////////////////////////////////////////////////////////////

  }  // namespace SpectralFD
}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_SpectralFD_SpectralFDMethod_hh

