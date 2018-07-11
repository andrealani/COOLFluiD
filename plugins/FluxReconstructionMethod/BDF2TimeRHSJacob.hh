#ifndef COOLFluiD_FluxReconstructionMethod_BDF2TimeRHSJacob_hh
#define COOLFluiD_FluxReconstructionMethod_BDF2TimeRHSJacob_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/DataSocketSink.hh"
#include "FluxReconstructionMethod/PseudoSteadyStdTimeRHSJacob.hh"
#include "FluxReconstructionMethod/FluxReconstructionSolverData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * @author Ray Vandenhoeck
 */
class BDF2TimeRHSJacob : public PseudoSteadyStdTimeRHSJacob {
public:

  /**
   * Constructor.
   */
  explicit BDF2TimeRHSJacob(const std::string& name);

  /**
   * Destructor.
   */
  virtual ~BDF2TimeRHSJacob();

  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  virtual void setup();

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

  /// storage of the past time rhs
  Framework::DataSocketSink< CFreal> socket_pastTimeRhs;

}; // class BDF2TimeRHSJacob

//////////////////////////////////////////////////////////////////////////////

  } // namespace FluxReconstructionMethod

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionMethod_BDF2TimeRHSJacob_hh
