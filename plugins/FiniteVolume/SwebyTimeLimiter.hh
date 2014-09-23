#ifndef COOLFluiD_Numerics_FiniteVolume_SwebyTimeLimiter_hh
#define COOLFluiD_Numerics_FiniteVolume_SwebyTimeLimiter_hh

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
class SwebyTimeLimiter : public TimeLimiter
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
  SwebyTimeLimiter(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~SwebyTimeLimiter();

  /**
   * Compute the limiter value
   */
  virtual CFreal computeTimeLimiterValue(const CFreal& slopeOLD, const CFreal& slopeNEW);


private: //data

  CFreal m_beta;

}; // end of class SwebyTimeLimiter

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_SwebyTimeLimiter_hh
