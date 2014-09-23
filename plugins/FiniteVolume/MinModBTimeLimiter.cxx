

#include "Environment/ObjectProvider.hh"

#include "FiniteVolume/FiniteVolume.hh"
#include "FiniteVolume/MinModBTimeLimiter.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<MinModBTimeLimiter,
                            TimeLimiter,
                            FiniteVolumeModule,
                            1>
MinModBTimeLimiterProvider("MinModB");

//////////////////////////////////////////////////////////////////////////////

void MinModBTimeLimiter::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal >("SlopeRatio","Threshold for the slope ratio.");
}

//////////////////////////////////////////////////////////////////////////////

MinModBTimeLimiter::MinModBTimeLimiter(const std::string& name) :
TimeLimiter(name)
{
  addConfigOptionsTo(this);

  _slopeRatioThreshold = 1.;
  setParameter("SlopeRatio",&_slopeRatioThreshold);

}

//////////////////////////////////////////////////////////////////////////////

MinModBTimeLimiter::~MinModBTimeLimiter()
{
}

//////////////////////////////////////////////////////////////////////////////

CFreal MinModBTimeLimiter::computeTimeLimiterValue(const CFreal& slopeOLD, const CFreal& slopeNEW)
{

  CFreal limiterValue = 1.;
  if(slopeNEW*slopeOLD > 0.)
  {
    if(_slopeRatioThreshold*fabs(slopeNEW) < fabs(slopeOLD)) limiterValue = slopeNEW/slopeOLD;
  }
  else{
    if(slopeOLD != 0.) limiterValue = 0.;
  }

  limiterValue = sqrt(limiterValue);


//   const CFreal r = slopeNEW/slopeOLD;
//   CFreal limiterValue = std::max(0., std::min(1., _slopeRatioThreshold*r));
//
//   if(slopeOLD == 0.) limiterValue = 1.;
//
//   limiterValue = sqrt(limiterValue);

  return limiterValue;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
