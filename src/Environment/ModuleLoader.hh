// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Environment_ModuleLoader_hh
#define COOLFluiD_Environment_ModuleLoader_hh

//////////////////////////////////////////////////////////////////////////////

#include <boost/filesystem/path.hpp>

#include "Common/NonCopyable.hh"
#include "Config/ConfigObject.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Environment {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a class that handles
/// dynamic module loader upon request of the user,
/// @author Tiago Quintino
class Environment_API ModuleLoader :
  public Config::ConfigObject,
  public Common::NonCopyable<ModuleLoader>   {

public:

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Default constructor
  ModuleLoader();

  /// Default destructor
  ~ModuleLoader();

  /// Configures the names of the libraries to load.
  /// @param args is the Config::ConfigArgs with the arguments to be parsed.
  virtual void configure ( Config::ConfigArgs& args );

  /// Loads external modules
  /// which were compiled as dynamic loadable libraries.
  /// It searchs the paths set in the DirPath singleton
  /// object.
  void loadExternalModules();

  /// Loads one module with the supplied name
  /// @param mod name of the module to load
  void loadModule(const std::string& mod);

  /// Sets the dir search paths
  /// @param paths vector with all paths to search when adding a module
  void setSearchPaths(std::vector< boost::filesystem::path >& paths);

private: // data

  /// configuration variable for the modules to be loaded
  std::vector<std::string> m_moduleNames;

}; // end of class ModuleLoader

//////////////////////////////////////////////////////////////////////////////

  } // namespace Environment

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Environment_ModuleLoader_hh
