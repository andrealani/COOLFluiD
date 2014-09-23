

#include "Environment/ObjectProvider.hh"

#include "FiniteVolume/FiniteVolume.hh"
#include "FiniteVolume/MinModATimeLimiter.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<MinModATimeLimiter,
                            TimeLimiter,
                            FiniteVolumeModule,
                            1>
MinModATimeLimiterProvider("MinModA");

//////////////////////////////////////////////////////////////////////////////

void MinModATimeLimiter::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal >("SlopeRatio","Threshold for the slope ratio.");
}

//////////////////////////////////////////////////////////////////////////////

MinModATimeLimiter::MinModATimeLimiter(const std::string& name) :
TimeLimiter(name)
{
  addConfigOptionsTo(this);

  _slopeRatioThreshold = 1.;
  setParameter("SlopeRatio",&_slopeRatioThreshold);

}

//////////////////////////////////////////////////////////////////////////////

MinModATimeLimiter::~MinModATimeLimiter()
{
}

//////////////////////////////////////////////////////////////////////////////

CFreal MinModATimeLimiter::computeTimeLimiterValue(const CFreal& slopeOLD, const CFreal& slopeNEW)
{

  CFreal limiterValue = 1.;
  if(slopeNEW*slopeOLD > 0.)
  {
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
