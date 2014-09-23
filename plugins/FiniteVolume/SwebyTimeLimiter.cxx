

#include "Environment/ObjectProvider.hh"

#include "FiniteVolume/FiniteVolume.hh"
#include "FiniteVolume/SwebyTimeLimiter.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<SwebyTimeLimiter,
                            TimeLimiter,
                            FiniteVolumeModule,
                            1>
SwebyTimeLimiterProvider("Sweby");

//////////////////////////////////////////////////////////////////////////////

void SwebyTimeLimiter::defineConfigOptions(Config::OptionList& options)
{
 options.addConfigOption< CFreal >("Beta","Beta value for the limiter.");
}

//////////////////////////////////////////////////////////////////////////////

SwebyTimeLimiter::SwebyTimeLimiter(const std::string& name) :
TimeLimiter(name)
{
  addConfigOptionsTo(this);

  m_beta = 1.5;
  setParameter("Beta",&m_beta);

}

//////////////////////////////////////////////////////////////////////////////

SwebyTimeLimiter::~SwebyTimeLimiter()
{
}

//////////////////////////////////////////////////////////////////////////////

CFreal SwebyTimeLimiter::computeTimeLimiterValue(const CFreal& slopeOLD, const CFreal& slopeNEW)
{

  const CFreal r = slopeOLD/slopeNEW;
  CFreal limiterValue = max((CFreal)0.,max(min(m_beta*r,(CFreal)1.),min(r,m_beta)));

  limiterValue = min(limiterValue, (CFreal)1.);

  return limiterValue;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
