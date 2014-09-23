// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Common_Stopwatch_hh
#define COOLFluiD_Common_Stopwatch_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/HourMinSec.hh"
#include "Common/TimePolicies.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Common {

//////////////////////////////////////////////////////////////////////////////

/// Performance Stopwatch  used for benchmarking codes.
/// Measure elapsed seconds
/// @author Tiago Quintino
template <typename TIMEPOLICY = Common::WallTime>
class Stopwatch {
public:

  /// Constructor
  Stopwatch();

  /// Destructor
  ~Stopwatch();

  /// Start timing from 0.00.
  /// @post m_running == true
  /// @post m_total == 0
  void start();

  /// Restart timing from 0.00. Same as calling start()
  /// @post m_running == true
  /// @post m_total == 0
  void restart();

  /// Reset the timer
  /// Clears the elapsed time.
  /// @post m_running == false
  /// @post m_total == 0
  void reset();

  /// Resumes counting. Doesn't clear the elapsed time.
  /// No effect if isRunning()
  /// @post m_running == true
  /// @post m_total >= 0
  void resume();

  /// Stop timing. Doesn't clear the elapsed time.
  /// No effect if isNotRunning()
  /// @post m_running == false
  /// @post m_total >= 0
  void stop();

  /// Read current time.
  /// @return current time in seconds
  CFreal read() const;
  
  /// Converst the current time to CFdouble
  /// @return current time in seconds
  operator CFreal () const { return read(); }
  
  /// Read current time.
  /// @return current time in Hour, Minutes and Seconds
  HourMinSec readTimeHMS();

  /// Checks if the Stopwatch is m_running
  /// @return TRUE if it is running
  bool isRunning() const;

  /// Checks if the Stopwatch is not running
  /// @return TRUE if it is not running
  bool isNotRunning() const;

private:

  /// Initializes the starting time
  void initStartTime();

  /// Takes the stoping time
  void takeStopTime();

  /// Adds the elapsed time to the accumulated time
  /// @post m_total is updated
  void accumulateTime();

private:

  /// flag to record if it is running
  bool   m_running;

  /// the accumulated time
  CFreal m_total;
  
  /// the way to count the time
  TIMEPOLICY impl;

}; // Class Stopwatch

//////////////////////////////////////////////////////////////////////////////

  } // namespace Common

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "Stopwatch.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Common_Stopwatch_hh

