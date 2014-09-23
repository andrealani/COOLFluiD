

#include "Environment/ObjectProvider.hh"

#include "FiniteVolume/FiniteVolume.hh"
#include "FiniteVolume/SuperbeeTimeLimiter.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<SuperbeeTimeLimiter,
                            TimeLimiter,
                            FiniteVolumeModule,
                            1>
SuperbeeTimeLimiterProvider("Superbee");

//////////////////////////////////////////////////////////////////////////////

void SuperbeeTimeLimiter::defineConfigOptions(Config::OptionList& options)
{
}

//////////////////////////////////////////////////////////////////////////////

SuperbeeTimeLimiter::SuperbeeTimeLimiter(const std::string& name) :
TimeLimiter(name)
{
  addConfigOptionsTo(this);
}

//////////////////////////////////////////////////////////////////////////////

SuperbeeTimeLimiter::~SuperbeeTimeLimiter()
{
}

//////////////////////////////////////////////////////////////////////////////

CFreal SuperbeeTimeLimiter::computeTimeLimiterValue(const CFreal& slopeOLD, const CFreal& slopeNEW)
{

  const CFreal r = slopeNEW/slopeOLD;
  CFreal limiterValue = 0.5 * max<CFreal>(max<CFreal>((CFreal)0.,min<CFreal>((CFreal)1.,2.*r)),min<CFreal>((CFreal)2.,r));
  
  return limiterValue;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
