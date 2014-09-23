// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_Maestro_hh
#define COOLFluiD_Framework_Maestro_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/NonCopyable.hh"
#include "Common/DynamicObject.hh"
#include "Common/OwnedObject.hh"
#include "Common/SelfRegistPtr.hh"

#include "Config/ConfigObject.hh"

#include "Environment/ConcreteProvider.hh"

#include "Framework/Framework.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

      class Simulator;
      class GlobalStopCriteria;

//////////////////////////////////////////////////////////////////////////////

/// Controls the flow of actions in the SubSystem's through signals
/// @author Tiago Quintino
/// @author Thomas Wuilbaut
class Framework_API Maestro :
      public Common::DynamicObject,
      public Common::OwnedObject,
      public Config::ConfigObject,
      public Common::NonCopyable<Maestro>   {


public: // typedefs

  typedef Environment::ConcreteProvider<Maestro,1> PROVIDER;
  typedef const std::string& ARG1;

public: // functions

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Default constructor without arguments.
  Maestro(const std::string& name);

  /// Default destructor.
  virtual ~Maestro();

  /// Configures this Simulator
  virtual void configure ( Config::ConfigArgs& args );

  /// Gets the Class name
  static std::string getClassName() { return "Maestro"; }

  /// Takes control of the Simulator
  void manage ( Common::SharedPtr<Simulator> sim );

protected: // data

  /// Name of StopCriteria to configure
  Common::SharedPtr<Simulator> m_sim;

  /// Save each iteration to different Tecplot file (with suffix _Globaliter#).
  bool m_append_iter;
  /// Restart from previous Global Iteration.
  bool m_restart_from_previous;

  /// Name of StopCriteria to configure
  std::string m_stopcriteria_str;
  /// The StopCriteriaController
  Common::SelfRegistPtr<GlobalStopCriteria> m_stopcriteria;

}; // end of class Maestro

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

Framework_Factory(Maestro) // define the factoty instance

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_Maestro_hh
