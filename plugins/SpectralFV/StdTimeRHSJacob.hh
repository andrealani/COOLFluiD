#ifndef COOLFluiD_SpectralFV_StdComputeTimeRhs_hh
#define COOLFluiD_SpectralFV_StdComputeTimeRhs_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/DataSocketSink.hh"
#include "SpectralFV/SpectralFVMethodData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFV {

//////////////////////////////////////////////////////////////////////////////

/**
 */
class StdTimeRHSJacob : public SpectralFVMethodCom {
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

  /// inverse volume fractions of CVs
  Common::SafePtr< std::vector< CFreal > > m_volFracCVs;

  /// number of equations in the physical model
  CFuint m_nbrEqs;

}; // class StdTimeRHSJacob

//////////////////////////////////////////////////////////////////////////////

  } // namespace SpectralFV

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_SpectralFV_StdComputeTimeRhs_hh
