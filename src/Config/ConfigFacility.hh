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

#ifndef COOLFluiD_Config_ConfigFacility_hh
#define COOLFluiD_Config_ConfigFacility_hh

//////////////////////////////////////////////////////////////////////////////

#include <set>

#include "Common/NonCopyable.hh"
#include "Config/ConfigArgs.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Config {

    class ConfigObject;

//////////////////////////////////////////////////////////////////////////////

/// The centralized configuration facility
/// @author Tiago Quintino
class Config_API ConfigFacility : public Common::NonCopyable<ConfigFacility> {

public: // methods

  /// @return the instance of this singleton
  static ConfigFacility& getInstance();

  /// Registers this object in the facility
  /// If object already exists nothing happens
  void registConfigObj( ConfigObject* cobj );

  /// Remove sthe object form the facility
  void unregistConfigObj( ConfigObject* cobj );

  /// Configure the dynamic options with the arguments passed
  void configureDynamicOptions ( ConfigArgs& args );

private: // methods

  /// Default constructor
  ConfigFacility();

  /// Default destructor
  ~ConfigFacility();

private: // data

  /// storage of the config objects that are registered
  std::set<ConfigObject*> m_config_objs;

}; // end of class ConfigFacility

//////////////////////////////////////////////////////////////////////////////

  } // namespace Config

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Config_ConfigFacility_hh
