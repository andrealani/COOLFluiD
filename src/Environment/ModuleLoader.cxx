// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Common/MemFunArg.hh"
#include "Common/OSystem.hh"
#include "Common/LibLoader.hh"

#include "Environment/DirPaths.hh"
#include "Environment/CFEnv.hh"
#include "Environment/ModuleRegisterBase.hh"
#include "Environment/ModuleRegistry.hh"
#include "Environment/ModuleLoader.hh"
#include "Environment/ModuleLoadFailedException.hh"
#include "Environment/SingleBehaviorFactory.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;

namespace COOLFluiD {

  namespace Environment {

//////////////////////////////////////////////////////////////////////////////

void ModuleLoader::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< std::vector<std::string> >("Libs","Module libraries to load.");
}

//////////////////////////////////////////////////////////////////////////////

ModuleLoader::ModuleLoader() : ConfigObject("Modules"),
  m_moduleNames()
{
   addConfigOptionsTo(this);
   setParameter("Libs",&m_moduleNames);
}

//////////////////////////////////////////////////////////////////////////////

ModuleLoader::~ModuleLoader()
{
}

//////////////////////////////////////////////////////////////////////////////

void ModuleLoader::setSearchPaths(std::vector< boost::filesystem::path >& paths)
{
  OSystem::getInstance().getLibLoader()->set_search_paths(paths);
}

//////////////////////////////////////////////////////////////////////////////

void ModuleLoader::loadExternalModules()
{
  std::vector< boost::filesystem::path > paths =  Environment::DirPaths::getInstance().getModulesDir();

  setSearchPaths(paths);

  // attempt to load all modules from the list
  if( !m_moduleNames.empty() )
  {
    std::for_each( m_moduleNames.begin(), m_moduleNames.end(), mem_fun_arg(*this,(&ModuleLoader::loadModule)) );
  }
  else
  {
    CFLog(NOTICE,"ModuleLoader: No external modules loaded\n");
  }
}

//////////////////////////////////////////////////////////////////////////////

void ModuleLoader::loadModule(const std::string& mod)
{
  CFAUTOTRACE;

  cf_assert( OSystem::getInstance().getLibLoader().isNotNull() );

  // we assume that library name is the same as the module being loaded
  // with the lib prefix stripped out
  std::string lmod;
  if (Common::StringOps::startsWith(mod,"lib"))
  {
    lmod = std::string(mod,3); // create a copy without the first 3 letters
  }
  else
  {
    lmod = mod;
  }

  // check if the module has already been loaded
  // if not, then try to load it
  if ( !Environment::CFEnv::getInstance().getModuleRegistry()->isRegistered(lmod) )
  {
    try
    {
      OSystem::getInstance().getLibLoader()->load_library(lmod);
    }
    catch ( LibLoaderException& exp )
    {
      throw ModuleLoadFailedException (FromHere(),exp.what());
    }
  }
  else
  {
    CFLog(NOTICE,"ModuleLoader: Module " << mod << " is already loaded.\n");
  }
}

//////////////////////////////////////////////////////////////////////////////

void ModuleLoader::configure ( Config::ConfigArgs& args )
{
  ConfigObject::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Environment

} // namespace COOLFluiD
