// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Common_HourMinSec_hh
#define COOLFluiD_Common_HourMinSec_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/COOLFluiD.hh"
#include "Common/Common.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Common {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a place holder for some timing data.
/// @author Andrea Lani
/// @author Tiago Quintino
class Common_API HourMinSec {
public:

  /// Default constructor without arguments
  HourMinSec();

  /// Constructor with nb of seconds.
  /// @param sec number of seconds to initializes
  HourMinSec(const CFdouble& sec);

  /// Default destructor
  ~HourMinSec();

  /// Conversion to std::string
  std::string str() const;

  /// Set the values corresponding to the given time (in sec)
  void set(const CFdouble& total);

public:
  /// number of seconds per hour
  static const CFdouble secPerHour;
  /// number of seconds per minute
  static const CFdouble secPerMin;
  /// number of micro seconds per second
  static const CFdouble usecPerSec;

private:
  /// number of hours
  CFdouble m_h;
  /// number of minutes
  CFdouble m_m;
  /// number of seconds
  CFdouble m_s;

}; // end HourMinSec

//////////////////////////////////////////////////////////////////////////////

  } // namespace Common

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Common_HourMinSec_hh

