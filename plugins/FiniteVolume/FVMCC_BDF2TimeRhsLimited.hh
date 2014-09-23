#ifndef COOLFluiD_Numerics_FiniteVolume_FVMCC_BDF2TimeRhsLimited_hh
#define COOLFluiD_Numerics_FiniteVolume_FVMCC_BDF2TimeRhsLimited_hh

//////////////////////////////////////////////////////////////////////////////

#include "FVMCC_PseudoSteadyTimeRhs.hh"
#include "Framework/DataSocketSource.hh"
#include "FiniteVolume/TimeLimiter.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
    class BlockAccumulator;
  }

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represent a command that computes the pseudo steady RHS
 * using standard cell center FVM schemes
 *
 * @author Thomas Wuilbaut
 *
 */
class FVMCC_BDF2TimeRhsLimited : public FVMCC_PseudoSteadyTimeRhs {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor.
   */
  explicit FVMCC_BDF2TimeRhsLimited(const std::string& name);

  /**
   * Destructor.
   */
  ~FVMCC_BDF2TimeRhsLimited();

  /**
   * Configures this command with supplied arguments.
   */
  virtual void configure ( Config::ConfigArgs& args );

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

private:

  /**
   * Compute the transformation matrix numerically
   */
  void computeNumericalTransMatrix(const CFuint iState);

  /**
   * Compute the analytical transformation matrix
   */
  void computeAnalyticalTransMatrix(const CFuint iState);


protected: // data

  /// storage of the past time rhs
  Framework::DataSocketSink< CFreal> socket_pastTimeRhs;

  /// time limiter
  Common::SelfRegistPtr<TimeLimiter> _timeLimiter;

  std::string _timeLimiterStr;

}; // class FVMCC_BDF2TimeRhsLimited

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_FVMCC_BDF2TimeRhsLimited_hh
