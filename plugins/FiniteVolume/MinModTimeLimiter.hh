#ifndef COOLFluiD_Numerics_FiniteVolume_MinModTimeLimiter_hh
#define COOLFluiD_Numerics_FiniteVolume_MinModTimeLimiter_hh

//////////////////////////////////////////////////////////////////////////////

#include "TimeLimiter.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class computes the limiting value for time using a Minmod limiter
 *
 * @author Thomas Wuilbaut
 *
 */
class MinModTimeLimiter : public TimeLimiter
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
  MinModTimeLimiter(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~MinModTimeLimiter();

  /**
   * Compute the limiter value
   */
  virtual CFreal computeTimeLimiterValue(const CFreal& slopeOLD, const CFreal& slopeNEW);

private:

  CFreal _slopeRatioThreshold;


}; // end of class MinModTimeLimiter

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_MinModTimeLimiter_hh
