#ifndef COOLFluiD_Numerics_FiniteVolume_MinMod2TimeLimiter_hh
#define COOLFluiD_Numerics_FiniteVolume_MinMod2TimeLimiter_hh

//////////////////////////////////////////////////////////////////////////////

#include "TimeLimiter.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class computes the limiting value for time using a MinMod2 limiter
 *
 * @author Thomas Wuilbaut
 *
 */
class MinMod2TimeLimiter : public TimeLimiter
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
  MinMod2TimeLimiter(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~MinMod2TimeLimiter();

  /**
   * Compute the limiter value
   */
  virtual CFreal computeTimeLimiterValue(const CFreal& slopeOLD, const CFreal& slopeNEW);

private:

  CFreal _slopeRatioThreshold;



}; // end of class MinMod2TimeLimiter

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_MinMod2TimeLimiter_hh
