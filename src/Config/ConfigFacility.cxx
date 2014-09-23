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

#include "Config/ConfigFacility.hh"
#include "Config/ConfigObject.hh"

#include "Common/CFLog.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;

namespace COOLFluiD {

  namespace Config {

//////////////////////////////////////////////////////////////////////////////

ConfigFacility::ConfigFacility()
{
}

//////////////////////////////////////////////////////////////////////////////

ConfigFacility::~ConfigFacility()
{
}

//////////////////////////////////////////////////////////////////////////////

ConfigFacility& ConfigFacility::getInstance()
{
  static ConfigFacility theConfigFacility;
  return theConfigFacility;
}

//////////////////////////////////////////////////////////////////////////////

void ConfigFacility::registConfigObj( ConfigObject* cobj )
{
  cf_assert( cobj != CFNULL );
  m_config_objs.insert(cobj);
}

//////////////////////////////////////////////////////////////////////////////

void ConfigFacility::unregistConfigObj( ConfigObject* cobj )
{
  cf_assert( cobj != CFNULL );
  m_config_objs.erase(cobj);
}

//////////////////////////////////////////////////////////////////////////////

void ConfigFacility::configureDynamicOptions( ConfigArgs& args )
{
  CFAUTOTRACE;

  for ( std::set<ConfigObject*>::iterator itr = m_config_objs.begin(),
                                         stop = m_config_objs.end();
        itr != stop; ++itr )
  {
    ConfigObject * cobj = *itr;
    cf_assert ( *itr != CFNULL );
    cobj->configureDynamic(args);
  }
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Config

} // namespace COOLFluiD
