// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_GlobalStopCriteria_hh
#define COOLFluiD_Framework_GlobalStopCriteria_hh

//////////////////////////////////////////////////////////////////////////////

#include "Config/ConfigObject.hh"
#include "Common/OwnedObject.hh"
#include "Environment/ConcreteProvider.hh"
#include "Framework/Framework.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This class defines the criteria to be satisfied to stop the Simulation
/// @author Tiago Quintino
/// @author Thomas Wuilbaut
class Framework_API GlobalStopCriteria : public Common::OwnedObject,
                           public Config::ConfigObject {
public:

  typedef Environment::ConcreteProvider<GlobalStopCriteria,1> PROVIDER;
  typedef const std::string& ARG1;

  /// Constructor
  GlobalStopCriteria(const std::string& name);

  /// Default destructor
  virtual ~GlobalStopCriteria();

  /// Gets the Class name
  static std::string getClassName()
  {
    return "GlobalStopCriteria";
  }

  /// Returns true if the simulation should end
  virtual bool isSatisfied() = 0;

}; // end of class GlobalStopCriteria

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

Framework_Factory(GlobalStopCriteria) // define the factoty instance

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_GlobalStopCriteria_hh
