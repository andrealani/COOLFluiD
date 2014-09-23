#ifndef COOLFluiD_SpectralFD_PseudoSteadyStdTimeDiagBlockJacob_hh
#define COOLFluiD_SpectralFD_PseudoSteadyStdTimeDiagBlockJacob_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/DataSocketSink.hh"
#include "SpectralFD/SpectralFDMethodData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

/**
 * @author Kris Van den Abeele
 */
class PseudoSteadyStdTimeDiagBlockJacob : public SpectralFDMethodCom {
public:

  /**
   * Constructor.
   */
  explicit PseudoSteadyStdTimeDiagBlockJacob(const std::string& name);

  /**
   * Destructor.
   */
  virtual ~PseudoSteadyStdTimeDiagBlockJacob();

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

  /// storage of the update coefficients
  Framework::DataSocketSink<CFreal> socket_updateCoeff;

  /// storage of the past states
//   Framework::DataSocketSink< Framework::State*> socket_pastStates;

  /// socket for diagonal block Jacobian matrices
  Framework::DataSocketSink< RealMatrix > socket_diagBlockJacobMatr;

  /// builder of cells
  Common::SafePtr<Framework::GeometricEntityPool<Framework::StdTrsGeoBuilder> > m_cellBuilder;

  /// current diagonal block Jacobian matrix
  RealMatrix* m_currDiagMatrix;

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

  /// element index
  CFuint m_elemIdx;

}; // class PseudoSteadyStdTimeDiagBlockJacob

//////////////////////////////////////////////////////////////////////////////

  } // namespace SpectralFD

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_SpectralFD_PseudoSteadyStdTimeDiagBlockJacob_hh
