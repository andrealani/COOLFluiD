

#include "Environment/ObjectProvider.hh"

#include "FiniteVolume/FiniteVolume.hh"
#include "FiniteVolume/MinModTimeLimiter.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<MinModTimeLimiter,
                            TimeLimiter,
                            FiniteVolumeModule,
                            1>
MinModTimeLimiterProvider("MinMod");

//////////////////////////////////////////////////////////////////////////////

void MinModTimeLimiter::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal >("SlopeRatio","Threshold for the slope ratio.");
}

//////////////////////////////////////////////////////////////////////////////

MinModTimeLimiter::MinModTimeLimiter(const std::string& name) :
TimeLimiter(name)
{
  addConfigOptionsTo(this);

  _slopeRatioThreshold = 1.;
  setParameter("SlopeRatio",&_slopeRatioThreshold);

}

//////////////////////////////////////////////////////////////////////////////

MinModTimeLimiter::~MinModTimeLimiter()
{
}

//////////////////////////////////////////////////////////////////////////////

CFreal MinModTimeLimiter::computeTimeLimiterValue(const CFreal& slopeOLD, const CFreal& slopeNEW)
{

  CFreal limiterValue = 1.;
  if(slopeNEW*slopeOLD > 0.)
  {
    if(_slopeRatioThreshold*fabs(slopeNEW) < fabs(slopeOLD)) limiterValue = slopeNEW/slopeOLD;
    if(fabs(slopeNEW) > _slopeRatioThreshold*fabs(slopeOLD)) limiterValue = slopeOLD/slopeNEW;
    //else limiterValue = slopeOLD/slopeNEW;
  }
  else{
    if(slopeOLD != 0.) limiterValue = 0.;
  }

  limiterValue = sqrt(limiterValue);

  return limiterValue;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
