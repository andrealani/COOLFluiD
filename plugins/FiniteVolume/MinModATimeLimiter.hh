#ifndef COOLFluiD_Numerics_FiniteVolume_MinModATimeLimiter_hh
#define COOLFluiD_Numerics_FiniteVolume_MinModATimeLimiter_hh

//////////////////////////////////////////////////////////////////////////////

#include "TimeLimiter.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class computes the limiting value for time using a MinModA limiter
 *
 * @author Thomas Wuilbaut
 *
 */
class MinModATimeLimiter : public TimeLimiter
{
public: // functions

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor
   */
  MinModATimeLimiter(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~MinModATimeLimiter();

  /**
   * Compute the limiter value
   */
  virtual CFreal computeTimeLimiterValue(const CFreal& slopeOLD, const CFreal& slopeNEW);

private:

  CFreal _slopeRatioThreshold;


}; // end of class MinModATimeLimiter

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_MinModATimeLimiter_hh
