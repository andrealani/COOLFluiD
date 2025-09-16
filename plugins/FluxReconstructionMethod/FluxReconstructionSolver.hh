// Copyright (C) 2016 KU Leuven, Belgium
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
    
    // Forward declarations
    class ConvBndCorrectionsRHSFluxReconstruction;
    class DiffBndCorrectionsRHSFluxReconstruction;
    
//////////////////////////////////////////////////////////////////////////////  
    

/**
 * This class implements a FluxReconstruction solver
 * 
 * @author Alexander Papen
 * @author Ray Vandenhoeck
 */
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
  
  /// Gets the Data aggregator of this space method
  /// @return SafePtr to the MethodData
  virtual Common::SafePtr< FluxReconstructionSolverData > getData()
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
  
  /// Defined the strategy list of this Method
  std::vector< Common::SafePtr< Framework::NumericalStrategy > > getStrategyList() const;

protected: // interface implementation functions

  /// Sets up data, commands and strategies of this Method
  virtual void setMethodImpl(){
    SpaceMethod::setMethodImpl();

   // setupCommandsAndStrategies();
    cf_assert(m_setup.isNotNull());
    m_setup->execute(); 
    setupCommandsAndStrategies();
  }

  /// Unsets the data, commands and strategies of this Method
  virtual void unsetMethodImpl(){
    cf_assert(m_unsetup.isNotNull());
    m_unsetup->execute();
    unsetupCommandsAndStrategies();

    SpaceMethod::unsetMethodImpl(); 
  }

  /// Extrapolates the states to the node positions
  virtual void extrapolateStatesToNodesImpl();

  /// Initialize the solution before starting the computation
  void initializeSolutionImpl(bool isRestart);

  /// Set matrix, right hand side and solve system
  virtual void computeSpaceResidualImpl(CFreal factor);

  /// Compute the time contribution to residual
  virtual void computeTimeResidualImpl(CFreal factor);

  /// Apply boundary conditions
  virtual void applyBCImpl();
  
  /// Apply boundary conditions for diffusive terms
  virtual void applyBCDiffImpl();
  
  /// Add source terms
  virtual void addSourceTermsImpl();

  /// Prepare to compute
  virtual void prepareComputationImpl();

  /// Postprocess the solution.
  virtual void postProcessSolutionImpl();

  /// Executed on "CF_ON_MESHADAPTER_BEFOREMESHUPDATE" event
  Common::Signal::return_t beforeMeshUpdateActionImpl(Common::Signal::arg_t eBefore);

  /// Executed on "CF_ON_MESHADAPTER_AFTERMESHUPDATE" event
  Common::Signal::return_t afterMeshUpdateActionImpl(Common::Signal::arg_t eAfter);
  
private: // functions

  void configureSourceTermCommands( Config::ConfigArgs& args );
  
  void configureInitCommands( Config::ConfigArgs& args );
  
  void configureBcCommands( Config::ConfigArgs& args );

private: // data

  ///The Setup command to use
  Common::SelfRegistPtr< FluxReconstructionSolverCom > m_setup;

  ///The UnSetup command to use
  Common::SelfRegistPtr< FluxReconstructionSolverCom > m_unsetup;
  
  /// The extrapolate command
  Common::SelfRegistPtr< FluxReconstructionSolverCom > m_extrapolate;

  /// Command used to prepare the computation
  Common::SelfRegistPtr< FluxReconstructionSolverCom > m_prepare;

  ///The convective solve command
  Common::SelfRegistPtr< FluxReconstructionSolverCom > m_convSolve;
  
  ///The diffusion solve command
  Common::SelfRegistPtr< FluxReconstructionSolverCom > m_diffSolve;
  
  /// The command that computes the contribution of the time discretization to the rhs
  Common::SelfRegistPtr< FluxReconstructionSolverCom > m_timeRHSJacob;
  
  /// Command used to limit a solution
  Common::SelfRegistPtr< FluxReconstructionSolverCom > m_limiter;
  
  /// Command used to preprocess a solution
  Common::SelfRegistPtr< FluxReconstructionSolverCom > m_preprocess;
  
  /// Command used to add artificial viscosity
  Common::SelfRegistPtr< FluxReconstructionSolverCom > m_artificialVisc;
  
  /// Command used to enforce physicality of the solution
  Common::SelfRegistPtr< FluxReconstructionSolverCom > m_physicality;
  
  /// Command used to compute the error of the solution
  Common::SelfRegistPtr< FluxReconstructionSolverCom > m_computeError;
  
  /// Command used to finalize the computation of the RHS
  Common::SelfRegistPtr< FluxReconstructionSolverCom > m_finalizeRHS;
  
  /// The commands to use for initializing the solution.
  std::vector< Common::SelfRegistPtr< FluxReconstructionSolverCom > > m_inits;

  /// The command to use for the action before the mesh is updated
  Common::SelfRegistPtr< FluxReconstructionSolverCom > _beforeMeshUpdate;

  /// The command to use for the action after the mesh has been updated
  Common::SelfRegistPtr< FluxReconstructionSolverCom > _afterMeshUpdate;


  ///The Setup string for configuration
  std::string m_setupStr;

  ///The UnSetup string for configuration
  std::string m_unsetupStr;
  
  /// The string for configuration of the extrapolate command
  std::string m_extrapolateStr;

  /// The string for configuration of the m_prepare command
  std::string m_prepareStr;

  /// The string for configuration of the m_convSolve command
  std::string m_convSolveStr;
  
  /// The string for configuration of the m_diffSolve command
  std::string m_diffSolveStr;
  
  /// The string for configuration of the m_limiter command
  std::string m_limiterStr;
  
  /// The string for configuration of the m_preprocess command
  std::string m_preprocessStr;
  
  /// The string for configuration of the m_artificialVisc command
  std::string m_artificialViscStr;
  
  /// The string for configuration of the m_physicality command
  std::string m_physicalityStr;
  
  /// The string for configuration of the commands that compute
  /// RHS (and optionally the Jacobian)
  std::string m_spaceRHSJacobStr;
  
  /// Boolean flag for using blending approach
  bool m_useBlending;
  
  /// The string for configuration of the m_timeRHSJacob command
  std::string m_timeRHSJacobStr;
  
  ///The computeError string for configuration
  std::string m_computeErrorStr;
  
  ///The finalizeRHS string for configuration
  std::string m_finalizeRHSStr;

  /// The solution initializing command types
  std::vector<std::string> m_initTypeStr;

  /// The solution initializing command names for configuration
  std::vector<std::string> m_initNameStr;
  
  /// The commands to use for applying the boundary conditions for the convective terms,
  /// with ConvBndFaceTermRHSFluxReconstruction as type
  std::vector< Common::SafePtr< ConvBndCorrectionsRHSFluxReconstruction > > m_bcs;

  /// The commands to use for applying the boundary conditions for the convective terms
  std::vector< Common::SelfRegistPtr< FluxReconstructionSolverCom > > m_bcsComs;
  
  /// The source term command types
  std::vector<std::string> m_srcTermTypeStr;

  /// The source term command names for configuration
  std::vector<std::string> m_srcTermNameStr;
  
  /// The commands for the source terms
  std::vector< Common::SelfRegistPtr< FluxReconstructionSolverCom > > m_srcTerms;
  
  /// The boundary condition command names for configuration
  std::vector<std::string> m_bcNameStr;
  
  /// The diffusion boundary condition command names for configuration
  std::vector<std::string> m_bcNameDiffStr;
  
  ///The data to share between FluxReconstructionSolverCom commands
  Common::SharedPtr< FluxReconstructionSolverData > m_data;
  
  /// The commands to use for applying the boundary conditions for the diffusive terms,
  /// with DiffBndCorrectionsRHSFluxReconstruction as type
  std::vector< Common::SafePtr< DiffBndCorrectionsRHSFluxReconstruction > > m_bcsDiff;

  /// The commands to use for applying the boundary conditions for the diffusive terms
  std::vector< Common::SelfRegistPtr< FluxReconstructionSolverCom > > m_bcsDiffComs;

  /// The string for the configuration of the action before mesh update
  std::string _beforeMeshUpdateStr;

  /// The string for the configuration of the action after mesh update
  std::string _afterMeshUpdateStr;

//////////////////////////////////////////////////////////////////////////////

}; // class FluxReconstructionSolver

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionMethod_FluxReconstructionSolver_hh

