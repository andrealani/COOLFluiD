#ifndef COOLFluiD_SpectralFD_AddBDF2TimeRHSInGivenCellSpectralFD_hh
#define COOLFluiD_SpectralFD_AddBDF2TimeRHSInGivenCellSpectralFD_hh

//////////////////////////////////////////////////////////////////////////////

#include "SpectralFD/AddPseudoSteadyStdTimeRHSInGivenCellSpectralFD.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

/**
 * @author Kris Van den Abeele
 */
class AddBDF2TimeRHSInGivenCellSpectralFD : public AddPseudoSteadyStdTimeRHSInGivenCellSpectralFD {
public:

  /**
   * Constructor.
   */
  explicit AddBDF2TimeRHSInGivenCellSpectralFD(const std::string& name);

  /**
   * Destructor.
   */
  virtual ~AddBDF2TimeRHSInGivenCellSpectralFD();

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

}; // class AddBDF2TimeRHSInGivenCellSpectralFD

//////////////////////////////////////////////////////////////////////////////

  } // namespace SpectralFD

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_SpectralFD_AddBDF2TimeRHSInGivenCellSpectralFD_hh
