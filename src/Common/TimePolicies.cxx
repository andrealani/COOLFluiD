// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Common/StringOps.hh"
#include "Common/TimePolicies.hh"
#include "Common/HourMinSec.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Common {

//////////////////////////////////////////////////////////////////////////////

#ifdef CF_HAVE_SYS_RESOURCE_H

TimePolicy_rusage::TimePolicy_rusage() : m_start(0.), m_stop(0.) {}

void TimePolicy_rusage::initStartTime()
{
  m_start = seconds();
}

void TimePolicy_rusage::takeStopTime()
{
  m_stop = seconds();
}

void TimePolicy_rusage::accumulateTime(CFreal& accTime)
{
  accTime += (m_stop - m_start);
}

CFdouble TimePolicy_rusage::getDelta() const
{
  return (seconds() - m_start);
}

CFdouble TimePolicy_rusage::seconds() const
{
  rusage usg;
  getrusage(RUSAGE_SELF, &usg);
  timeval usertime = usg.ru_utime;
  return (usertime.tv_sec + (usertime.tv_usec * HourMinSec::usecPerSec));
}

#endif // CF_HAVE_SYS_RESOURCE_H

//////////////////////////////////////////////////////////////////////////////

#if defined (CF_HAVE_TIME_H) || defined (CF_HAVE_SYS_TIME_H)

TimePolicy_cclock::TimePolicy_cclock() : m_start(0.), m_stop(0.) {}

void TimePolicy_cclock::initStartTime()
{
  m_start = seconds();
}

void TimePolicy_cclock::takeStopTime()
{
  m_stop = seconds();
}

void TimePolicy_cclock::accumulateTime(CFreal& accTime)
{
  accTime += (m_stop - m_start);
}

CFdouble TimePolicy_cclock::getDelta() const
{
  return (seconds() - m_start);
}

CFdouble TimePolicy_cclock::seconds() const
{
  // not so accurate but very portable
  const CFdouble secs_per_tick = 1.0 / CLOCKS_PER_SEC;
  return ( static_cast<CFdouble>(clock()) ) * secs_per_tick;
}

#endif // CF_HAVE_TIME_H || CF_HAVE_SYS_TIME_H

//////////////////////////////////////////////////////////////////////////////

#ifdef CF_HAVE_GETTIMEOFDAY

TimePolicy_gettimeofday::TimePolicy_gettimeofday()
{
  gettimeofday(&m_start,0);
  m_stop = m_start;
}

void TimePolicy_gettimeofday::initStartTime()
{
  gettimeofday(&m_start,0);
}

void TimePolicy_gettimeofday::takeStopTime()
{
  gettimeofday(&m_stop,0);
}

void TimePolicy_gettimeofday::accumulateTime(CFreal& accTime)
{
  timeval diff;
  subTimeval(diff, m_stop, m_start);
  accTime += toDouble(diff);
}

CFdouble TimePolicy_gettimeofday::getDelta() const
{
  timeval diff, now;
  gettimeofday (&now, 0);
  bool nega = subTimeval(diff, now, m_start);
  return (nega ? toDouble(diff) * -1 : toDouble(diff));
}

CFdouble TimePolicy_gettimeofday::toDouble(const timeval & result) const
{
  return (result.tv_sec + static_cast<CFdouble>(result.tv_usec / 1000000.0L));
}

bool TimePolicy_gettimeofday::subTimeval (timeval & Res, const timeval & X1, const timeval & Y1) const
{
  timeval X = X1;
  timeval Y = Y1;
  if (X.tv_usec < Y.tv_usec)
  {
      int nsec = (Y.tv_usec - X.tv_usec) / 1000000;
      Y.tv_usec += 1000000 * nsec;
      Y.tv_sec -= nsec;
  }
  if (X.tv_usec - Y.tv_usec > 1000000)
  {
      int nsec = (Y.tv_usec - X.tv_usec) / 1000000;
      Y.tv_usec += 1000000 * nsec;
      Y.tv_sec -= nsec;
  }
  Res.tv_sec = X.tv_sec - Y.tv_sec;
  Res.tv_usec = X.tv_usec - Y.tv_usec;

  return X.tv_sec < Y.tv_sec;
}

#endif // CF_HAVE_GETTIMEOFDAY

//////////////////////////////////////////////////////////////////////////////

  } // namespace Common

} // namespace COOLFluiD

