#ifndef COOLFluiD_Numerics_FiniteVolume_FVMCC_ComputeRHSSingleState_hh
#define COOLFluiD_Numerics_FiniteVolume_FVMCC_ComputeRHSSingleState_hh

//////////////////////////////////////////////////////////////////////////////

#include "CellCenterFVMData.hh"
#include "Framework/DataSocketSink.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represent a command that computes the RHS using
 * standard cell center FVM schemes
 *
 * @author Andrea Lani
 *
 */
class FVMCC_ComputeRHSSingleState : public CellCenterFVMCom {
public:

  /**
   * Constructor.
   */
  explicit FVMCC_ComputeRHSSingleState(const std::string& name);

  /**
   * Destructor.
   */
  virtual ~FVMCC_ComputeRHSSingleState();

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
   * Execute Processing actions
   */
  virtual void execute();

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();
  
protected: // functions
  
  /**
   * Set the integrator data
   */
  void setFaceIntegratorData();
  
  /**
   * Compute the face areas and other face related data
   */
  void computeFaceRelatedData();
  
  /**
   * Get the factor multiplying the residual
   */
  CFreal getResFactor() // const
  {
    return getMethodData().getResFactor();
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
  
protected:
  
  /// storage of the face normals
  Framework::DataSocketSink<CFreal> socket_normals;
  
  /// storage of the face areas
  Framework::DataSocketSink<CFreal> socket_faceAreas;
  
  /// storage of the States
  Framework::DataSocketSink < Framework::State* , Framework::GLOBAL > socket_states;
  
  /// storage of the cell volumes
  Framework::DataSocketSink<CFint> socket_isOutward;
  
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
  
  // flag telling if there is a diffusive term
  bool _hasDiffusiveTerm;
  
  // pointer to the current face
  Framework::GeometricEntity* _currFace;
  
  /// builder of cells
  Framework::GeometricEntityPool<Framework::CellTrsGeoBuilder> _cellBuilder;
  
  // face index inside the loop
  CFuint _faceIdx;
  
  // medium radius of the face
  CFreal _rMid;
  
  // temporary flux
  RealVector _flux;
  
  // temporary diffusive flux
  RealVector _dFlux;
  
  // temporary flux
  RealVector _rFlux;
  
  /// temporary unit normal
  RealVector _tempUnitNormal;
  
  /// storage for temporary interpolated coordinates
  /// in the (face) quadrature points
  Framework::Node* _faceCoord;
  
  /// reconstructed extra variables
  std::vector<RealVector*> _rExtraVars;
  
}; // class FVMCC_ComputeRHSSingleState

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_FVMCC_ComputeRHSSingleState_hh

