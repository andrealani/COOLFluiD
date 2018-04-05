#ifndef COOLFluiD_Numerics_FiniteVolume_TimeLimiter_hh
#define COOLFluiD_Numerics_FiniteVolume_TimeLimiter_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/ConcreteProvider.hh"
#include "Common/OwnedObject.hh"
#include "Common/NonCopyable.hh"
#include "Config/ConfigObject.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class offers a basic interface to limit
 * in time
 *
 * @author Thomas Wuilbaut
 *
 */
class TimeLimiter : public Common::NonCopyable<TimeLimiter>,
                    public Common::OwnedObject,
                    public Config::ConfigObject
{
public: // functions

  /// the provider of this type of classes
  typedef Environment::ConcreteProvider<TimeLimiter,1> PROVIDER;

  /// the first argument in the creation should be the name
  typedef const std::string& ARG1;

  /**
   * Constructor
   */
  explicit TimeLimiter(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~TimeLimiter();

  /**
   * Compute the limiter value
   */
  virtual CFreal computeTimeLimiterValue(const CFreal& slopeOLD, const CFreal& slopeNEW) = 0;

  /**
   * Gets the Class name
   */
  static std::string getClassName()
  {
    return "TimeLimiter";
  }

}; // end of class TimeLimiter

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_TimeLimiter_hh
