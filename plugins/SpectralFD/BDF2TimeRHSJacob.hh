#ifndef COOLFluiD_SpectralFD_BDF2TimeRHSJacob_hh
#define COOLFluiD_SpectralFD_BDF2TimeRHSJacob_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/DataSocketSink.hh"
#include "SpectralFD/PseudoSteadyStdTimeRHSJacob.hh"
#include "SpectralFD/SpectralFDMethodData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

/**
 * @author Kris Van den Abeele
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

  } // namespace SpectralFD

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_SpectralFD_BDF2TimeRHSJacob_hh
