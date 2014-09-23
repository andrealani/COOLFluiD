#ifndef COOLFluiD_SpectralFD_BDF3TimeDiagBlockJacob_hh
#define COOLFluiD_SpectralFD_BDF3TimeDiagBlockJacob_hh

//////////////////////////////////////////////////////////////////////////////

#include "SpectralFD/PseudoSteadyStdTimeDiagBlockJacob.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

/**
 * @author Matteo Parsani
 * @author Kris Van den Abeele
 */
class BDF3TimeDiagBlockJacob : public PseudoSteadyStdTimeDiagBlockJacob {
public:

  /**
   * Constructor.
   */
  explicit BDF3TimeDiagBlockJacob(const std::string& name);

  /**
   * Destructor.
   */
  virtual ~BDF3TimeDiagBlockJacob();

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

  /// pointer to the 3 steps time marching scheme parameters (BDF3 with variable time step)
  Common::SafePtr< RealVector > m_3StepsTMSparams;

}; // class BDF3TimeDiagBlockJacob

//////////////////////////////////////////////////////////////////////////////

  } // namespace SpectralFD

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_SpectralFD_BDF3TimeDiagBlockJacob_hh
