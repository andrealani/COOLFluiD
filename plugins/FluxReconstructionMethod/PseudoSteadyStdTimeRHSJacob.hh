#ifndef COOLFluiD_FluxReconstructionMethod_PseudoSteadyStdTimeRHSJacob_hh
#define COOLFluiD_FluxReconstructionMethod_PseudoSteadyStdTimeRHSJacob_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/DataSocketSink.hh"
#include "FluxReconstructionMethod/FluxReconstructionSolverData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * @author Kris Van den Abeele
 */
class PseudoSteadyStdTimeRHSJacob : public FluxReconstructionSolverCom {
public:

  /**
   * Constructor.
   */
  explicit PseudoSteadyStdTimeRHSJacob(const std::string& name);

  /**
   * Destructor.
   */
  virtual ~PseudoSteadyStdTimeRHSJacob();

  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  virtual void setup();

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

  /**
   * Adds the contribution of the time residual to the rhs and the jacobian
   */
  virtual void addTimeResidual();

protected:

  /// storage of volumes (for SD --> Jacobian determinants)
  Framework::DataSocketSink< CFreal> socket_volumes;

  /// storage of the rhs
  Framework::DataSocketSink< CFreal> socket_rhs;

  /// storage of the update coefficients
  Framework::DataSocketSink<CFreal> socket_updateCoeff;

  /// storage of the past states
  Framework::DataSocketSink< Framework::State*> socket_pastStates;

  /// builder of cells
  Common::SafePtr<Framework::GeometricEntityPool<Framework::StdTrsGeoBuilder> > m_cellBuilder;

  /// pointer to the linear system solver
  Common::SafePtr<Framework::LinearSystemSolver> m_lss;

  /// pointer to Jacobian matrix
  Common::SafePtr<Framework::LSSMatrix> m_jacobMatrix;

  /// number of equations in the physical model
  CFuint m_nbrEqs;

  /// pointer to cellstates
  std::vector< Framework::State* >* m_cellStates;

  /// solution point mapped coordinates
  Common::SafePtr< std::vector< RealVector > > m_solPntsLocalCoords;

  /// diagonal value on the LSS matrix
  RealVector m_diagValues;

  /// boolean telling whether computation is unsteady
  bool m_isUnsteady;

}; // class PseudoSteadyStdTimeRHSJacob

//////////////////////////////////////////////////////////////////////////////

  } // namespace FluxReconstructionMethod

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionMethod_PseudoSteadyStdTimeRHSJacob_hh
