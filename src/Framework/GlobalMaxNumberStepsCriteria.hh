// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_GlobalMaxNumberStepsCriteria_hh
#define COOLFluiD_GlobalMaxNumberStepsCriteria_hh

//////////////////////////////////////////////////////////////////////////////

#include "GlobalStopCriteria.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

  /// This class represents the Global stop criteria dictated by a maximum
  /// number of iterations
  /// @author Tiago Quintino
class Framework_API GlobalMaxNumberStepsCriteria : public GlobalStopCriteria {
public:

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Default constructor without arguments
  /// @see GlobalStopCriteria()
  GlobalMaxNumberStepsCriteria(const std::string& name);

  /// Default destructor
  ~GlobalMaxNumberStepsCriteria();

  /// Take the combined value from all CPU's and decide if the simulation
  /// should stop
  virtual bool isSatisfied();

private:

  // maximum number of steps to perform
  CFuint m_maxNbSteps;

}; // end of class GlobalMaxNumberStepsCriteria

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_GlobalMaxNumberStepsCriteria_hh
