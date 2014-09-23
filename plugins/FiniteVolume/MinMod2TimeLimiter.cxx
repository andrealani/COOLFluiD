

#include "Environment/ObjectProvider.hh"

#include "FiniteVolume/FiniteVolume.hh"
#include "FiniteVolume/MinMod2TimeLimiter.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<MinMod2TimeLimiter,
                            TimeLimiter,
                            FiniteVolumeModule,
                            1>
MinMod2TimeLimiterProvider("MinMod2");

//////////////////////////////////////////////////////////////////////////////

void MinMod2TimeLimiter::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal >("SlopeRatio","Threshold for the slope ratio.");
}

//////////////////////////////////////////////////////////////////////////////

MinMod2TimeLimiter::MinMod2TimeLimiter(const std::string& name) :
TimeLimiter(name)
{
  addConfigOptionsTo(this);

  _slopeRatioThreshold = 1.;
  setParameter("SlopeRatio",&_slopeRatioThreshold);

}

//////////////////////////////////////////////////////////////////////////////

MinMod2TimeLimiter::~MinMod2TimeLimiter()
{
}

//////////////////////////////////////////////////////////////////////////////

CFreal MinMod2TimeLimiter::computeTimeLimiterValue(const CFreal& slopeOLD, const CFreal& slopeNEW)
{

  CFreal limiterValue = 1.;

  if(_slopeRatioThreshold*fabs(slopeNEW) < fabs(slopeOLD)) limiterValue = fabs(slopeNEW)/fabs(slopeOLD);

  limiterValue = sqrt(limiterValue);

  return limiterValue;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
