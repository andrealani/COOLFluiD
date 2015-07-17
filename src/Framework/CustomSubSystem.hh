// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_CustomSubSystem_hh
#define COOLFluiD_Framework_CustomSubSystem_hh

//////////////////////////////////////////////////////////////////////////////

#include "StandardSubSystem.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// A CustomSubSystem is a concrete implementation of the SubSystem
/// that allows to run a customized subsystem
/// @author Thomas Wuilbaut
class Framework_API CustomSubSystem : public StandardSubSystem {
public: // functions

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// The default constructor without arguments.
  /// @see ~CustomSubSystem()
  CustomSubSystem(const std::string& name);

  /// Destructor
  /// @see CustomSubSystem()
  virtual ~CustomSubSystem();

  /// Run (Process) Phase. All the big number crunching work goes here.
  virtual void run();

  /// Run one function of the sequence
  virtual void runSequenceID(const CFuint ID);
  
 protected: // function
  
  /// Helper function to run the functionwith the method passed
  void executeMethodFunction (const std::string& method_name,
			      const std::string& function_name);
  
protected: // data

  /// collection of strings used to describe the sequence of commands to run
  std::vector<std::string> m_runSequenceStr;

}; // class CustomSubSystem

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_CustomSubSystem_hh
