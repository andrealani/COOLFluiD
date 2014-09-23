#ifndef COOLFluiD_SpectralFD_AddBDF3TimeRHSInGivenCellSpectralFD_hh
#define COOLFluiD_SpectralFD_AddBDF3TimeRHSInGivenCellSpectralFD_hh

//////////////////////////////////////////////////////////////////////////////

#include "SpectralFD/AddPseudoSteadyStdTimeRHSInGivenCellSpectralFD.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

/**
 * @author Matteo Parsani
 * @author Kris Van den Abeele
 */
class AddBDF3TimeRHSInGivenCellSpectralFD : public AddPseudoSteadyStdTimeRHSInGivenCellSpectralFD {
public:

  /**
   * Constructor.
   */
  explicit AddBDF3TimeRHSInGivenCellSpectralFD(const std::string& name);

  /**
   * Destructor.
   */
  virtual ~AddBDF3TimeRHSInGivenCellSpectralFD();

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

  /// storage of the past past time rhs
  Framework::DataSocketSink< CFreal> socket_pastPastTimeRhs;

  /// pointer to the 3 steps time marching scheme parameters (BDF3 with variable time step)
  Common::SafePtr< RealVector > m_3StepsTMSparams;

}; // class AddBDF3TimeRHSInGivenCellSpectralFD

//////////////////////////////////////////////////////////////////////////////

  } // namespace SpectralFD

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_SpectralFD_AddBDF3TimeRHSInGivenCellSpectralFD_hh
