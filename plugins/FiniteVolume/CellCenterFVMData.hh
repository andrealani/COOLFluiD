#ifndef COOLFluiD_Numerics_FiniteVolume_CellCenterFVMData_hh
#define COOLFluiD_Numerics_FiniteVolume_CellCenterFVMData_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/SafePtr.hh"
#include "Common/OwnedObject.hh"
#include "Config/ConfigObject.hh"
#include "Framework/MethodCommand.hh"
#include "Framework/ConvectiveVarSet.hh"
#include "Framework/DiffusiveVarSet.hh"
#include "Framework/LinearSystemSolver.hh"
#include "Framework/ConvergenceMethod.hh"
#include "Framework/GeoDataComputer.hh"
#include "Framework/FluxSplitter.hh"
#include "Framework/Limiter.hh"
#include "Framework/VarSetMatrixTransformer.hh"
#include "Framework/MultiMethodHandle.hh"
#include "Framework/GeometricEntityPool.hh"
#include "Framework/TrsGeoWithNodesBuilder.hh"
#include "Framework/FaceTrsGeoBuilder.hh"
#include "Framework/FaceCellTrsGeoBuilder.hh"
#include "Framework/CellTrsGeoBuilder.hh"
#include "Framework/SpaceMethodData.hh"
#include "Framework/JacobianLinearizer.hh"
#include "Framework/VolumeIntegrator.hh"
#include "FiniteVolume/ComputeSourceTermFVMCC.hh"
#include "FiniteVolume/FVMCC_EquationFilter.hh"
#include "Framework/NodalStatesExtrapolator.hh"
#include "FiniteVolume/FVMCC_PolyRec.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  
  namespace Numerics {
    
    namespace FiniteVolume {
      
      class ComputeDiffusiveFlux;
      class DerivativeComputer;
      class FVMCC_BC;
      
//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a Data Object that is accessed by the different
 * CellCenterFVMCom 's that compose the CellCenterFVM.
 *
 * @todo there is missing documentation in this class.
 *
 * @author Andrea Lani
 * @author Mehmet Sarp Yalim
 *
 */
class CellCenterFVMData : public Framework::SpaceMethodData {
public:
    
  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Default constructor without arguments.
   */
  CellCenterFVMData(Common::SafePtr<Framework::Method> owner);

  /**
   * Destructor.
   */
  ~CellCenterFVMData();

  /// Configure the data from the supplied arguments.
  /// @param args configuration arguments
  virtual void configure ( Config::ConfigArgs& args );

  /// Sets up the method data
  virtual void setup();

  /// Unsets the method data
  virtual void unsetup();

  /**
   * @return true if the simulation is axisymmetric
   */
  bool isAxisymmetric() const
  {
    return _isAxisymm;
  }

  /**
   * @return true if there is a source term
   */
  bool hasSourceTerm() const
  {
    return _hasSourceTerm;
  }

  /**
   * Set the residual factor
   */
  void setResFactor(const CFreal& factor)
  {
    _resFactor = factor;
  }

  /**
   * Get the residual factor
   */
  CFreal getResFactor() const
  {
    return _resFactor;
  }
 
  /**
   * Get the computer of geometric data
   * @return a reference to the computer of geometric data
   */
  Common::SafePtr<Framework::GeoDataComputer<CellCenterFVMData> > getGeoDataComputer() const
  {
    cf_assert(_geoDataComputer.isNotNull());
    return _geoDataComputer.getPtr();
  }
  
  /**
   * Get the flux splitter
   * @return a reference to the flux splitter
   */
  Common::SafePtr<Framework::FluxSplitter<CellCenterFVMData> > getFluxSplitter() const
  {
    cf_assert(_fluxSplitter.isNotNull());
    return _fluxSplitter.getPtr();
  }

  /**
   * Get the diffusive flux computer
   * @return a reference to the diffusive flux computer
   */
  Common::SafePtr<ComputeDiffusiveFlux> getDiffusiveFluxComputer() const
  {
    cf_assert(_diffusiveFlux.isNotNull());
    return _diffusiveFlux.getPtr();
  }

  /**
   * Get the polynomial reconstructor
   */
  Common::SafePtr<FVMCC_PolyRec> getPolyReconstructor() const
  {
    cf_assert(_polyRec.isNotNull());
    return _polyRec.getPtr();
  }
  
  /**
   * Get the limiter
   */
  Common::SafePtr<Framework::Limiter<CellCenterFVMData> > getLimiter() const
  {
    cf_assert(_limiter.isNotNull());
    return _limiter.getPtr();
  }
  
  /**
   * Get the jacobian linearizer
   */
  Common::SafePtr<Framework::JacobianLinearizer> getJacobianLinearizer() const
  {
    cf_assert(_linearizer.isNotNull());
    return _linearizer.getPtr();
  }
  
  /**
   * Get the nodal extrapolator
   */
  Common::SafePtr<Framework::NodalStatesExtrapolator<CellCenterFVMData> > 
  getNodalStatesExtrapolator() const
  {
    cf_assert(_nStatesExtrapolator.isNotNull());
    return _nStatesExtrapolator.getPtr();
  }

  /**
   * Get the DerivativeComputer
   */
  Common::SafePtr<DerivativeComputer> getDerivativeComputer() const
  {
    return _derivComputer.getPtr();
  }

  /**
   * Flag telling if the solution variables have to be reconstructed
   */
  bool reconstructSolVars() const
  {
    return _reconstructSolVars;
  }

  /**
   * Use the analytical jacobian for the convective fluxes
   */
  bool useAnalyticalConvJacob() const
  {
    return _useAnalyticalConvJacob;
  }

  /**
   * Set the analytical jacobian for the convective fluxes
   */
  void setUseAnalyticalConvJacob(bool flag)
  {
    _useAnalyticalConvJacob = flag;
  }

  /**
   * Get the linearization variables name
   */
  std::string getLinearVarStr() const
  {
    return _linearVarStr;
  }
  
  /**
   * Get the reconstruction variables name
   */
  std::string getReconstructVarStr() const
  {
    return _reconstructVarStr;
  }
  
  /**
   * Get the matrix transformer from solution to update variables
   * starting from update variables
   */
  Common::SafePtr<Framework::VarSetMatrixTransformer>
  getSolToUpdateInUpdateMatTrans() const
  {
    cf_assert(_solToUpdateInUpdateMatTrans.isNotNull());
    return _solToUpdateInUpdateMatTrans.getPtr();
  }

  /**
   * Get the matrix transformer from update to solution
   * starting from update variables
   */
  Common::SafePtr<Framework::VarSetMatrixTransformer>
  getUpdateToSolutionInUpdateMatTrans() const
  {
    cf_assert(_updateToSolutionInUpdateMatTrans.isNotNull());
    return _updateToSolutionInUpdateMatTrans.getPtr();
  }
  
  /**
   * Get the matrix transformer from solution to linear
   * starting from update variables
   */
  Common::SafePtr<Framework::VarSetTransformer>
  getSolutionToLinearVecTrans() const
  {
    cf_assert(_solutionToLinearVecTrans.isNotNull());
    return _solutionToLinearVecTrans.getPtr();
  }

  /**
   * Get the vector transformer from update to solution variables
   */
  Common::SafePtr<Framework::VarSetTransformer>
  getUpdateToSolutionVecTrans() const
  {
    cf_assert(_updateToSolutionVecTrans.isNotNull());
    return _updateToSolutionVecTrans.getPtr();
  } 
  
  /**
   * Get the vector transformer from update to reconstruction variables
   */
  Common::SafePtr<Framework::VarSetTransformer>
  getUpdateToReconstructionVecTrans() const
  {
    cf_assert(_updateToReconstrVecTrans.isNotNull());
    return _updateToReconstrVecTrans.getPtr();
  }
  
  /**
   * Get the vector transformer from reconstruction to update variables
   */
  Common::SafePtr<Framework::VarSetTransformer>
  getReconstructionToUpdateVecTrans() const
  {
    cf_assert(_reconstrToUpdateVecTrans.isNotNull());
    return _reconstrToUpdateVecTrans.getPtr();
  }
  
  /**
   * Sets the ConvergenceMethod for this SpaceMethod to use
   * @pre the pointer to ConvergenceMethod is not constant to
   *      allow dynamic_casting
   */
  void setConvergenceMethod(Framework::MultiMethodHandle<Framework::ConvergenceMethod> convMtd)
  {
    _convergenceMtd = convMtd;
  }

  /**
   * Get the ConvergenceMethod
   */
  Framework::MultiMethodHandle<Framework::ConvergenceMethod> getConvergenceMethod() const
  {
    cf_assert(_convergenceMtd.isNotNull());
    return _convergenceMtd;
  }
  
  /// get the list of equation filters
  Common::SafePtr< std::vector<Common::SelfRegistPtr<Framework::EquationFilter<CellCenterFVMData> > > >
  getEquationFilters()
  {
    cf_assert(_eqFilters.size() > 0);
    return &_eqFilters;
  }
  
  /**
   * Get the source term computer
   * @return the pointer to the source term computer
   */
  Common::SafePtr< std::vector<Common::SelfRegistPtr<Framework::ComputeSourceTerm<CellCenterFVMData> > > >
  getSourceTermComputer()
  {
    cf_assert(_stComputer.size() > 0);
    return &_stComputer;
  }

  /**
   * Get the source term computer names
   * @return the vector of source term computer names
   */
  std::vector<std::string> getSourceTermComputerNames()
  {
    cf_assert(_stComputer.size() > 0);
    return _stComputerStr;
  }
  
  /**
   * @return the GeometricEntity builder
   */
  Common::SafePtr<Framework::GeometricEntityPool<Framework::FaceTrsGeoBuilder> >
  getFaceTrsGeoBuilder()
  {
    return &_faceTrsGeoBuilder;
  }
  
  /**
   * @return the GeometricEntity builder
   */
  Common::SafePtr<Framework::GeometricEntityPool<Framework::FaceCellTrsGeoBuilder> >
  getFaceCellTrsGeoBuilder()
  {
    return &_faceCellTrsGeoBuilder;
  }
  
  /**
   * @return the GeometricEntity builder for cells
   */
  Common::SafePtr<Framework::GeometricEntityPool<Framework::CellTrsGeoBuilder> >
  getCellTrsGeoBuilder()
  {
    return &_cellTrsGeoBuilder;
  }

 /**
   * @return the GeometricEntity builder
   */
  Common::SafePtr<Framework::GeometricEntityPool<Framework::TrsGeoWithNodesBuilder> >
  getGeoWithNodesBuilder()
  {
    return &_geoWithNodesBuilder;
  }

  /**
   * Tell if a variable has to be applied for the residual
   */
  bool isResidualTransformationNeeded() const
  {
    return (_updateVarStr != _solutionVarStr);
  }

  /**
   * Set the preprocessing phase
   */
  void setIsPreProcessedSolution(bool flag)
  {
    _isPreProcessedSolution = flag;
  }

  /**
   * Tell if the preprocessing phase has been accomplished
   */
  bool isPreProcessedSolution() const
  {
    return _isPreProcessedSolution;
  }
  
  /**
   * Set the initialization phase
   */
  void setIsInitializationPhase(bool flag)
  {
    _isInitializationPhase = flag;
  }

  /**
   * Tell if the initialization phase is ongoing
   */
  bool isInitializationPhase() const
  {
    return _isInitializationPhase;
  }
  
  /**
   * Preprocess BC flag on
   */
  bool& preProcessBCFlag() {return _preProcessBCFlag;}
  
  /**
   * Tell the names of the TRSs on which ghost states should be placed on the face
   */
  const std::vector<std::string>& getTRSsWithGhostsOnFace() const
  {
    return _trssWithGhostsOnFace;
  }

  /**
   * Tell the names of the TRSs for which no BCs must be applied
   */
  const std::vector<std::string>& getTRSsWithNoBC() const
  {
    return _trssWithNoBC;
  }

  /**
   * Gets the Class name
   */
  static std::string getClassName()
  {
    return "CellCenterFVM";
  }

  /**
   * Gets the opposite local face ID (in 2D/3D structured meshes or
   * boundary layers)
   */
  CFuint getOppositeIFace(CFuint iFace, CFuint dim,
			  CFuint nbCellNodes) const;
  
  /// get the face unit normal 
  RealVector& getUnitNormal() {return _unitNormal;}
  
  /// get a pointer to the current face
  Framework::GeometricEntity*& getCurrentFace() {return _currFace;}
  
  /// get flag to tell to use the average flux on the current face
  bool getUseAverageFlux() const {return _useAverageFlux;}
    
  /// set flag to tell to use the average flux on the current face
  void setUseAverageFlux(bool useAverageFlux) {_useAverageFlux = useAverageFlux;}
  
  /// get flag to tell to construct left and right cells
  bool getBuildAllCells() const {return _buildAllCells;}
  
  /// set flag to tell to use the average flux on the current face
  void setBuildAllCells(bool buildAllCells) {_buildAllCells = buildAllCells;}
  
  /*
   * Gets a map for ghost cell states
   */
  Common::SafePtr<Common::CFMap<Framework::State*,CFuint> > getMapGhostStateIDs()
  {
     return &_mapGhostStateIDs;
  }

  /// get the volume integrator
  Common::SafePtr<Framework::VolumeIntegrator> getVolumeIntegrator()
  {
    return &_volumeIntegrator;
  }
  
  /**
   * Set the list of boundary conditions
   */
  void setBCList(const std::vector<Common::SelfRegistPtr<Framework::MethodCommand<CellCenterFVMData> > >& bcList);
  
  /**
   * Get the BC map, associating the TRS index with the corresponding BC command
   * @pre this cannot be inlined for compilation reasons: FVMCC_BC needs to be forward declared
   */
  Common::SafePtr<Common::CFMap<CFuint, FVMCC_BC*> > getMapBC();
  
private:
  
  /**
   * Configures the GeoDataComputer
   */
  void configureGeoDataComputer ( Config::ConfigArgs& args );
  
  /**
   * Configures the FluxSplitter
   */
  void configureFluxSplitter ( Config::ConfigArgs& args );

  /**
   * Configures the ComputeDiffusiveFlux
   */
  void configureDiffusiveFluxComputer ( Config::ConfigArgs& args );

  /**
   * Configures the PolyReconstructor
   */
  void configurePolyReconstructor ( Config::ConfigArgs& args );

  /**
   * Configures the JacobianLinearizer
   */
  void configureJacobianLinearizer ( Config::ConfigArgs& args );
  
  /**
   * Configures the Limiter
   */
  void configureLimiter ( Config::ConfigArgs& args );

  /**
   * Configures the NodalStatesExtrapolator
   */
  void configureNodalStatesExtrapolator ( Config::ConfigArgs& args );

  /**
   * Configures the source term computer
   */
  void configureSourceTermComputer ( Config::ConfigArgs& args );

  /// Configures the equation filters
  void configureEquationFilters ( Config::ConfigArgs& args );
  
  /**
   * Configures the variable set transfomers
   */
  void configureVarSetTransformers ( Config::ConfigArgs& args );

private:
  
  /// flag indicating if solution initialization is ongoing
  bool  _isInitializationPhase;

  /// flag indicating if preprocesing has been accomplished
  bool  _isPreProcessedSolution;
  
  /// volume integrator
  Framework::VolumeIntegrator _volumeIntegrator;
  
  /// Diffusive flux computer
  Common::SelfRegistPtr<ComputeDiffusiveFlux> _diffusiveFlux;

  /// object computing the derivatives
  Common::SelfRegistPtr<DerivativeComputer> _derivComputer;

  /// Geometric data computer
  Common::SelfRegistPtr<Framework::GeoDataComputer<CellCenterFVMData> > _geoDataComputer;
  
  /// Flux splitter
  Common::SelfRegistPtr<Framework::FluxSplitter<CellCenterFVMData> > _fluxSplitter;

  /// Polynomial reconstructor
  Common::SelfRegistPtr<FVMCC_PolyRec> _polyRec;
  
  /// Limiter
  Common::SelfRegistPtr<Framework::Limiter<CellCenterFVMData> > _limiter;

  /// nodal states extrapolator
  Common::SelfRegistPtr<Framework::NodalStatesExtrapolator<CellCenterFVMData> > _nStatesExtrapolator;

  /// Convergence Method
  Framework::MultiMethodHandle<Framework::ConvergenceMethod> _convergenceMtd;

  /// Equation filters
  std::vector<Common::SelfRegistPtr<Framework::EquationFilter<CellCenterFVMData> > > _eqFilters;
  
  /// Source term computer
  std::vector<Common::SelfRegistPtr<Framework::ComputeSourceTerm<CellCenterFVMData> > > _stComputer;

  /// Matrix transformer from solution to update variables
  /// starting from update variables
  Common::SelfRegistPtr<Framework::VarSetMatrixTransformer>
  _solToUpdateInUpdateMatTrans;
  
  /// Matrix transformer from update to solution variables
  /// starting from update variables
  Common::SelfRegistPtr<Framework::VarSetMatrixTransformer>
  _updateToSolutionInUpdateMatTrans;
  
  /// Vector transformer from update to solution variables
  Common::SelfRegistPtr<Framework::VarSetTransformer> _updateToSolutionVecTrans;
  
  /// Vector transformer from solution to linear variables
  Common::SelfRegistPtr<Framework::VarSetTransformer> _solutionToLinearVecTrans;
  
  /// Vector transformer from update to reconstruction variables
  Common::SelfRegistPtr<Framework::VarSetTransformer> _updateToReconstrVecTrans;
  
  /// Vector transformer from reconstruction to update variables
  Common::SelfRegistPtr<Framework::VarSetTransformer> _reconstrToUpdateVecTrans;
  
  /// jacobian linearizer
  Common::SelfRegistPtr<Framework::JacobianLinearizer> _linearizer;
  
  /// builder for faces
  Framework::GeometricEntityPool<Framework::FaceTrsGeoBuilder> _faceTrsGeoBuilder;
  
  /// builder for faces with cells
  Framework::GeometricEntityPool<Framework::FaceCellTrsGeoBuilder> _faceCellTrsGeoBuilder;
  
  /// builder for cells
  Framework::GeometricEntityPool<Framework::CellTrsGeoBuilder> _cellTrsGeoBuilder;

  // builder of GeometricEntity's with Node's
  Framework::GeometricEntityPool<Framework::TrsGeoWithNodesBuilder> _geoWithNodesBuilder;

  /// current face
  Framework::GeometricEntity* _currFace;
  
  ///The command to use for computing the boundary conditions.
  Common::CFMap<CFuint, FVMCC_BC*> _bcMap;
  
  /// adimensional normal
  RealVector _unitNormal;
  
  /// flag that tells to perform preprocessing for BCs (one action per all BCs at once)
  bool _preProcessBCFlag;
  
  /// flag that tells to use the average flux on the current face
  bool _useAverageFlux;

  /// flag that tells if the simulation has a source term
  bool _hasSourceTerm;
  
  /// flag that tells if both left and right cells must be built
  bool _buildAllCells;
    
  /// residual factor
  CFreal _resFactor;

  /// linearizazion variable set name
  std::string _linearVarStr;
  
  /// reconstruction variable set name
  std::string _reconstructVarStr;
  
  /// string for the configuration of the geometric data computer
  std::string _geoDataComputerStr;

  /// string for the configuration of the flux splitter
  std::string _fluxSplitterStr;

  /// string for the configuration of the diffusive flux computer
  std::string _diffusiveFluxStr;

  /// string for the configurqtion of the derivatives computer
  std::string _derivComputerStr;

  /// string for the configuration of the polynomial reconstructor
  std::string _polyRecStr;

  /// string for the configuration of the limiter
  std::string _limiterStr;

  /// string for the configuration of the numerical integrator QuadratureType
  std::string _integratorQuadratureStr;
  
  /// string for the configuration of the numerical integrator Order
  std::string _integratorOrderStr;
  
  /// object computing the solution extrapolation in the nodal states
  std::string _nStatesExtrapolatorStr;

  /// strings for configuration of the equation filters
  std::vector<std::string> _eqFiltersStr;
  
  /// strings for configuration of the source term computer
  std::vector<std::string> _stComputerStr;

  ///  name of the TRSs on which ghost should be placed on the face
  std::vector<std::string> _trssWithGhostsOnFace;
  
  ///  name of TRSs for which a BC doesn't have to be applied
  std::vector<std::string> _trssWithNoBC;
  
  /// flag that tells if the simulation is axisymmetric
  bool _isAxisymm;

  /// Use analytical convective flux jacobian
  bool _useAnalyticalConvJacob;

  /// reconstruct the solution (conservative) variables
  bool _reconstructSolVars;
  
  /// GhostStates / IDs Map
  Common::CFMap<Framework::State*, CFuint> _mapGhostStateIDs;

 }; // end of class CellCenterFVMData

//////////////////////////////////////////////////////////////////////////////

/// Definition of a command for FiniteVolume
typedef Framework::MethodCommand<CellCenterFVMData> CellCenterFVMCom;

/// Definition of a command provider for FiniteVolume
typedef Framework::MethodCommand<CellCenterFVMData>::PROVIDER CellCenterFVMComProvider;

/// Definition of a Strategy for FiniteElement
typedef Framework::MethodStrategy<CellCenterFVMData> CellCenterFVMStrategy;

/// Definition of a Strategy provider for FiniteElement
typedef Framework::MethodStrategy<CellCenterFVMData>::PROVIDER CellCenterFVMStrategyProvider;

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_CellCenterFVMData_hh
