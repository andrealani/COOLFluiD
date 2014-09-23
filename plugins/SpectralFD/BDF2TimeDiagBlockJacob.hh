#ifndef COOLFluiD_SpectralFD_BDF2TimeDiagBlockJacob_hh
#define COOLFluiD_SpectralFD_BDF2TimeDiagBlockJacob_hh

//////////////////////////////////////////////////////////////////////////////

#include "SpectralFD/PseudoSteadyStdTimeDiagBlockJacob.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

/**
 * @author Kris Van den Abeele
 */
class BDF2TimeDiagBlockJacob : public PseudoSteadyStdTimeDiagBlockJacob {
public:

  /**
   * Constructor.
   */
  explicit BDF2TimeDiagBlockJacob(const std::string& name);

  /**
   * Destructor.
   */
  virtual ~BDF2TimeDiagBlockJacob();

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

}; // class BDF2TimeDiagBlockJacob

//////////////////////////////////////////////////////////////////////////////

  } // namespace SpectralFD

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_SpectralFD_BDF2TimeDiagBlockJacob_hh
