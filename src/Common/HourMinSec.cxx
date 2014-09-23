// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include <sstream>

#include "Common/HourMinSec.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Common {

//////////////////////////////////////////////////////////////////////////////

const CFdouble HourMinSec::secPerHour = 1. / 3600.;

const CFdouble HourMinSec::secPerMin = 1. / 60.;

const CFdouble HourMinSec::usecPerSec = 1. / 1000000.;

//////////////////////////////////////////////////////////////////////////////

HourMinSec::HourMinSec() : m_h(0.), m_m(0.), m_s(0.)
{
}

//////////////////////////////////////////////////////////////////////////////

HourMinSec::HourMinSec(const CFdouble& sec)
{
  set(sec);
}

//////////////////////////////////////////////////////////////////////////////

HourMinSec::~HourMinSec()
{
}

//////////////////////////////////////////////////////////////////////////////

std::string HourMinSec::str() const
{
  std::ostringstream oss;
  if (m_h > 0) {
    oss << m_h << " h ";
  }
  if (m_m > 0) {
    oss << m_m << " min ";
  }
  oss << m_s << " sec";
  return oss.str();
}

//////////////////////////////////////////////////////////////////////////////

void HourMinSec::set(const CFdouble& total)
{
  m_h = floor(total * secPerHour);
  m_m = floor((total - m_h * 3600.) * secPerMin);
  m_s = total - m_h * 3600. - m_m * 60.;
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Common

} // namespace COOLFluiD

