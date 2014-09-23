#ifndef COOLFluiD_Numerics_FiniteVolume_MCTimeLimiter_hh
#define COOLFluiD_Numerics_FiniteVolume_MCTimeLimiter_hh

//////////////////////////////////////////////////////////////////////////////

#include "TimeLimiter.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class computes the limiting value for time using a
 * Monotonized Central limiter
 *
 * @author Thomas Wuilbaut
 *
 */
class MCTimeLimiter : public TimeLimiter
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
  MCTimeLimiter(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~MCTimeLimiter();

  /**
   * Compute the limiter value
   */
  virtual CFreal computeTimeLimiterValue(const CFreal& slopeOLD, const CFreal& slopeNEW);


}; // end of class MCTimeLimiter

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_MCTimeLimiter_hh
