#ifndef COOLFluiD_Numerics_FiniteVolume_FVMCC_ComputeRHS_hh
#define COOLFluiD_Numerics_FiniteVolume_FVMCC_ComputeRHS_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolume/CellCenterFVMData.hh"
#include "Framework/DataSocketSink.hh"
#include "FiniteVolume/ComputeDiffusiveFlux.hh"
#include "FiniteVolume/FVMCC_PolyRec.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  
  namespace MathTools { class MatrixInverter; }
  
  namespace Numerics {

    namespace FiniteVolume {
      class FVMCC_BC;
      
//////////////////////////////////////////////////////////////////////////////

/**
 * This class represent a command that computes the RHS using
 * standard cell center FVM schemes
 *
 * @author Andrea Lani
 *
 */
class FVMCC_ComputeRHS : public CellCenterFVMCom {
public:

  /**
   * Constructor.
   */
  explicit FVMCC_ComputeRHS(const std::string& name);

  /**
   * Destructor.
   */
  virtual ~FVMCC_ComputeRHS();

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  virtual void setup();

  /**
   * Un Setup private data and data of the aggregated classes
   * in this command after processing phase
   */
  virtual void unsetup();

  /**
   * Configures the command.
   */
  virtual void configure ( Config::ConfigArgs& args );

  /**
   * Execute Processing actions
   */
  virtual void execute();

 
  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();
    
protected:
  
  /// Restore the backed up left states
  virtual void restoreState(CFuint iCell) {}
  
  /// Initialize the computation of RHS
  virtual void initializeComputationRHS();
  
  /// Compute the jacobian of the RHS
  virtual void computeRHSJacobian();
  
  /// Finalize the computation of RHS
  virtual void finalizeComputationRHS();
  
  /// Compute between convective and diffusive term
  virtual void computeInterConvDiff();
  
  /**
   * Get the factor multiplying the residual
   */
  CFreal getResFactor() {return getMethodData().getResFactor();}
  
  /**
   * Compute the physical data in the states
   */
  virtual void computePhysicalData()
  {
    computeStatesData();
  }
  
  /**
   * Compute the physical data in the states
   */
  void computeStatesData()
  {
    std::vector<Framework::State*>& states = _polyRec->getExtrapolatedValues();
    std::vector<RealVector>& pdata = _polyRec->getExtrapolatedPhysicaData();
    _reconstrVar->computePhysicalData(*states[0], pdata[0]);
    _reconstrVar->computePhysicalData(*states[1], pdata[1]);
  }

  /**
   * Back up the physical data in the states
   */
  void computeAndBackUpStatesData()
  {
    std::vector<Framework::State*>& states = _polyRec->getExtrapolatedValues();
    std::vector<RealVector>& pdata = _polyRec->getExtrapolatedPhysicaData();
    std::vector<RealVector>& pdataBkp = _polyRec->getBackupPhysicaData();
    
    _reconstrVar->computePhysicalData(*states[0], pdata[0]);
    _reconstrVar->computePhysicalData(*states[1], pdata[1]);
    pdataBkp[0] = pdata[0];
    pdataBkp[1] = pdata[1];
    
    // _reconstrVar->computeAndBackUpStatesData
    //       (_rStates, _rExtraVars, 2*_nbQPointsInFace[_faceIdx]);
  }
   
  /// Compute the update factors for left and right states in non axisymmetric case
  void computeNoAxiUpFactors()
  {  
    _upFactor[0] = getResFactor();
    _upFactor[1] = -1.0;
    _upStFactor[0] = _upStFactor[1] = -getResFactor();
  }
  
  /// Compute the update factors for left and right states in axisymmetric case
  void computeAxiUpFactors()
  {  
    _upFactor[0] = _rMid*_invr[0]*getResFactor();
    _upFactor[1] = (-_invr[1]/_invr[0]);
    _upStFactor[0] = -getResFactor()*_invr[0];
    _upStFactor[1] = -getResFactor()*_invr[1];
  }
  
  /**
   * Set the integrator data
   */
  virtual void setFaceIntegratorData();

  /**
   * Compute the contribution of the current face to the RHS
   */
  virtual void updateRHS();

  /**
   * Compute the source term
   */
  virtual void computeSourceTerm();
  
  /**
   * Transform residual
   */
  void transformResidual();

  /**
   * Get the flux to compute the jacobian
   */
  RealVector& getJacobianFlux()
  {
    return _flux;
  }
   
  /// Tell if source term has to be computed
  bool computeSourceTermJacob(CFuint idx, CFuint iCell, const std::vector<CFuint>& ids) const
  {
    return (idx == iCell) && _sourceJacobOnCell[idx] && (ids.size() > 0);
  }
  
  /// Tell if source term has to be computed
  bool computeSourceTermJacob(CFuint iCell, const std::vector<CFuint>& ids) const
  {
    return _sourceJacobOnCell[iCell] && (ids.size() > 0);
  }
  
  /// Compute the transformation matrix dP/dU analytically
  const RealMatrix& computeAnalyticalTransMatrix(const Framework::State& state) 
  {
    _solutionToUpdateMatTrans->setMatrix(state);
    return *_solutionToUpdateMatTrans->getMatrix();
  }
  
  /// Compute the transformation matrix dP/dU numerically
  RealMatrix& computeNumericalTransMatrix(Framework::State& state);
  
protected:
  
  /// flags for cells
  Framework::DataSocketSink<bool> socket_cellFlag;

  /// storage of the face normals
  Framework::DataSocketSink<CFreal> socket_normals;
  
  /// storage of the face areas
  Framework::DataSocketSink<CFreal> socket_faceAreas;
  
  /// storage of the rhs
  Framework::DataSocketSink<CFreal> socket_rhs;

  /// storage of the update coefficients
  Framework::DataSocketSink<CFreal> socket_updateCoeff;
  
  /// storage of the States
  Framework::DataSocketSink < Framework::State* , Framework::GLOBAL > socket_states;

  /// storage of the cell volumes
  Framework::DataSocketSink<CFreal> socket_volumes;

  /// storage of the cell volumes
  Framework::DataSocketSink<CFint> socket_isOutward;
  
  /// storage of the nodal states
  Framework::DataSocketSink<RealVector> socket_nstates;

  /// storage of the limiter
  Framework::DataSocketSink<CFreal> socket_limiter;
  
  /// storage of the ghost states
  Framework::DataSocketSink<Framework::State*> socket_gstates;

  /// storage of the nodes
  Framework::DataSocketSink < Framework::Node* , Framework::GLOBAL > socket_nodes;
  
  /// flux splitter
  Common::SafePtr<Framework::FluxSplitter<CellCenterFVMData> > _fluxSplitter;  
  /// diffusive flux computer
  Common::SafePtr<ComputeDiffusiveFlux> _diffusiveFlux;

  /// update variable set
  Common::SafePtr<Framework::ConvectiveVarSet> _reconstrVar;

  /// diffusive variable set
  Common::SafePtr<Framework::DiffusiveVarSet> _diffVar;
  
  /// polynomial reconstructor
  Common::SafePtr<FVMCC_PolyRec> _polyRec;
  
  /// nodal extrapolator
  Common::SafePtr<Framework::NodalStatesExtrapolator<CellCenterFVMData> > _nodalExtrapolator;
  
  /// Matrix transformer from solution to update variables
  /// starting from update variables
  Common::SafePtr<Framework::VarSetMatrixTransformer> _solutionToUpdateMatTrans;

  /// Vector transformer from update to solution variables
  Common::SafePtr<Framework::VarSetTransformer> _updateToSolutionVecTrans;
  
  // source term computers
  typedef Framework::ComputeSourceTerm<CellCenterFVMData> SourceTerm;
  Common::SafePtr< std::vector<Common::SelfRegistPtr<SourceTerm> > > _stComputers;
  
  // equation filters
  typedef Framework::EquationFilter<CellCenterFVMData> EqFilter;
  Common::SafePtr< std::vector<Common::SelfRegistPtr<EqFilter> > > _eqFilters;
  
  // pointer to the current face
  Framework::GeometricEntity* _currFace;
  
  // current BC
  FVMCC_BC* _currBC;
  
  // flag telling if there is a diffusive term
  bool _hasDiffusiveTerm;
  
  // flag telling if diffusion is active
  bool _isDiffusionActive;
  
  // face index inside the loop
  CFuint _faceIdx;
  
  // medium radius of the face
  CFreal _rMid;
  
  /// update factor for fluxes LR
  RealVector _upFactor;
  
  /// update factor for source terms LR 
  RealVector _upStFactor;
  
  /// inverse of radius for left and right states
  RealVector _invr;
  
  // temporary flux
  RealVector _flux;
  
  // temporary diffusive flux
  RealVector _dFlux;
  
  // temporary flux
  RealVector _rFlux;
  
  // dummy jacobian
  RealMatrix _jacobDummy;
  
  // inverse of the dummy jacobian
  RealMatrix _invJacobDummy;
  
  // temporary source terms
  std::vector<std::vector<RealVector> > _source;
  
  // temporary source term jacobians
  std::vector<std::vector<RealMatrix> > _sourceJacobian;
  
  // IDs corresponding to source term computers with numerical jacobian
  std::vector<CFuint> _stNumJacobIDs;
  
  // IDs corresponding to source term computers with analytical jacobian
  std::vector<CFuint> _stAnJacobIDs;
  
  /// flags for telling to compute the source term jacobian 
  std::vector<bool> _sourceJacobOnCell;
  
  /// storage for temporary extrapolated solution
  /// in the quadrature points
  Common::SafePtr<Framework::FluxSplitterData> _fluxData;
  
  /// temporary unit normal
  RealVector _tempUnitNormal;
    
  /// reconstructed extra variables
  std::vector<RealVector*> _rExtraVars;
  
  /// temporary data for holding the matrix inverter
  std::auto_ptr<MathTools::MatrixInverter>  _inverter;
  
  /// flag forcing to freeze the diffusive coefficients
  /// while doing numerical perturbation of the jacobians
  bool _freezeDiffCoeff;
  
  /// flag forcing to compute analytically the jacobian of the
  /// diffusive fluxes
  bool _analyticalDiffJacob;
  
  /// flag forcing to extrapolate the solution to nodes consistently
  /// (normally useless)
  bool _extrapolateInNodes;
  
  /// flag telling if to use analytical transformation matrix
  bool _useAnalyticalMatrix;
  
}; // class FVMCC_ComputeRHS

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_FVMCC_ComputeRHS_hh

