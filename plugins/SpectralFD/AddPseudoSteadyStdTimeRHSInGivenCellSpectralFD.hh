#ifndef COOLFluiD_SpectralFD_AddPseudoSteadyStdTimeRHSInGivenCellSpectralFD_hh
#define COOLFluiD_SpectralFD_AddPseudoSteadyStdTimeRHSInGivenCellSpectralFD_hh

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
class AddPseudoSteadyStdTimeRHSInGivenCellSpectralFD : public SpectralFDMethodCom {
public:

  /**
   * Constructor.
   */
  explicit AddPseudoSteadyStdTimeRHSInGivenCellSpectralFD(const std::string& name);

  /**
   * Destructor.
   */
  virtual ~AddPseudoSteadyStdTimeRHSInGivenCellSpectralFD();

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

  /// socket for rhs of current set of states
  Framework::DataSocketSink< CFreal > socket_rhsCurrStatesSet;

  /// socket for update coefficients of current set of states
  Framework::DataSocketSink< CFreal > socket_updateCoeff;

  /// socket for current states set index
  Framework::DataSocketSink< CFint > socket_statesSetIdx;

  /// storage of the past states
  Framework::DataSocketSink< Framework::State*> socket_pastStates;

  /// builder of cells
  Common::SafePtr<Framework::GeometricEntityPool<Framework::StdTrsGeoBuilder> > m_cellBuilder;

  /// variable for cell
  Framework::GeometricEntity* m_cell;

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

  /// index of element type
  CFuint m_iElemType;

  /// current iteration
  CFuint m_currIter;

}; // class AddPseudoSteadyStdTimeRHSInGivenCellSpectralFD

//////////////////////////////////////////////////////////////////////////////

  } // namespace SpectralFD

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_SpectralFD_AddPseudoSteadyStdTimeRHSInGivenCellSpectralFD_hh
