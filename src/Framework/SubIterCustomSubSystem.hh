// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_SubIterCustomSubSystem_hh
#define COOLFluiD_Framework_SubIterCustomSubSystem_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/CustomSubSystem.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// A SubIterCustomSubSystem is a concrete implementation of the SubSystem
/// that allows to run a customized subsystem with subiterations
/// @author Thomas Wuilbaut

class Framework_API SubIterCustomSubSystem : public CustomSubSystem {
public: // functions

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// The default constructor without arguments.
  /// @see ~SubIterCustomSubSystem()
  SubIterCustomSubSystem(const std::string& name);

  /// Destructor
  /// @see SubIterCustomSubSystem()
  virtual ~SubIterCustomSubSystem();

  /// Configures this Simualtion.
  /// Sets up the data for this Object.
  virtual void configure( Config::ConfigArgs& args );

  /// Run (Process) Phase. All the big number crunching work goes here.
  virtual void run();

  /// Run one function of the sequence
  virtual void runSequenceID(const CFuint ID);

protected: // data

  ///flag for stopping the subIter loop
  bool _endSubIteration;

  CFuint m_nbSubIterations;

  /// Type of StopCondition to configure
  std::string m_subIterStopConditionStr;

  /// Name of StopCondition to configure
  std::string m_subIterStopConditionNameStr;

  /// The StopConditionController for the subIteration
  std::auto_ptr<StopConditionController> m_subIterStopCondControler;


}; // class SubIterCustomSubSystem

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_SubIterCustomSubSystem_hh
