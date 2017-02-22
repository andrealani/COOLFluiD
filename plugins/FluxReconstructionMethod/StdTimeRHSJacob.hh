#ifndef COOLFluiD_FluxReconstructionMethod_StdComputeTimeRhs_hh
#define COOLFluiD_FluxReconstructionMethod_StdComputeTimeRhs_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/DataSocketSink.hh"
#include "FluxReconstructionMethod/FluxReconstructionSolverData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 */
class StdTimeRHSJacob : public FluxReconstructionSolverCom {
public:

  /**
   * Constructor.
   */
  explicit StdTimeRHSJacob(const std::string& name);

  /**
   * Destructor.
   */
  virtual ~StdTimeRHSJacob();

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

  /// builder of cells
  Common::SafePtr<Framework::GeometricEntityPool<Framework::StdTrsGeoBuilder> > m_cellBuilder;

  /// pointer to the linear system solver
  Common::SafePtr<Framework::LinearSystemSolver> m_lss;

  /// storage of the update coefficients
  Framework::DataSocketSink<CFreal> socket_updateCoeff;

  /// solution point mapped coordinates
  Common::SafePtr< std::vector< RealVector > > m_solPntsLocalCoords;

  /// number of equations in the physical model
  CFuint m_nbrEqs;

}; // class StdTimeRHSJacob

//////////////////////////////////////////////////////////////////////////////

  } // namespace FluxReconstructionMethod

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionMethod_StdComputeTimeRhs_hh
