#ifndef COOLFluiD_Numerics_FiniteVolume_MinModBTimeLimiter_hh
#define COOLFluiD_Numerics_FiniteVolume_MinModBTimeLimiter_hh

//////////////////////////////////////////////////////////////////////////////

#include "TimeLimiter.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class computes the limiting value for time using a MinModB limiter
 *
 * @author Thomas Wuilbaut
 *
 */
class MinModBTimeLimiter : public TimeLimiter
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
  MinModBTimeLimiter(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~MinModBTimeLimiter();

  /**
   * Compute the limiter value
   */
  virtual CFreal computeTimeLimiterValue(const CFreal& slopeOLD, const CFreal& slopeNEW);

private:

  CFreal _slopeRatioThreshold;

}; // end of class MinModBTimeLimiter

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_MinModBTimeLimiter_hh
