// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_MaxNumberStepsCondition_hh
#define COOLFluiD_MaxNumberStepsCondition_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/StopCondition.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

  /// This class represents the stop condition dictated by a maximum
  /// number of iterations
  /// @author Andrea Lani
class Framework_API MaxNumberStepsCondition : public StopCondition {
public: // methods

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Default constructor without arguments
  /// @see StopCondition()
  MaxNumberStepsCondition(const std::string& name);

  /// Default destructor
  ~MaxNumberStepsCondition();

  /// returns true if we have to evaluate the stopcondition on all CPU domains.
  virtual bool IsGlobal () const ;

  /// Take the combined value from all CPU's and decide if the simulation
  /// should stop
  virtual bool isAchieved (const ConvergenceStatus& status) ;

private: // data

  // maximum number of steps to perform
  CFuint _maxNbSteps;

}; // end of class MaxNumberStepsCondition

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_MaxNumberStepsCondition_hh
