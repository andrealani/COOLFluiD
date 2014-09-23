// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include <sstream>

#include "Common/StringOps.hh"
#include "Config/ConfigArgs.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;

namespace COOLFluiD {

  namespace Config {

//////////////////////////////////////////////////////////////////////////////

ConfigArgs::ConfigArgs() {}

//////////////////////////////////////////////////////////////////////////////

ConfigArgs::~ConfigArgs() {}

//////////////////////////////////////////////////////////////////////////////

ConfigArgs& ConfigArgs::operator=(const std::map<std::string,std::string>& value)
{
  std::map<std::string,std::string>::operator=(value);
  return *this;
}

//////////////////////////////////////////////////////////////////////////////

std::string ConfigArgs::str () const
{
  std::ostringstream oss;
  for (ConfigArgs::const_iterator it = begin(); it != end(); ++it)
  {
    oss << it->first << " " << it->second << "\n";
  }
  return oss.str();
}

//////////////////////////////////////////////////////////////////////////////

void ConfigArgs::consume(const std::vector<ConfigKey>& keys)
{
  std::vector<ConfigKey>::const_iterator itr = keys.begin();
  for ( ; itr != keys.end(); ++itr)
    erase( *itr );
}

//////////////////////////////////////////////////////////////////////////////

void ConfigArgs::remove_filter(const ConfigKey& filtkey)
{
  std::vector<ConfigKey> remove_args;
  for (ConfigArgs::iterator it = begin(); it != end(); ++it)
  {
    const ConfigKey& this_key = it->first;
    if (StringOps::startsWith(this_key,filtkey))
      remove_args.push_back(this_key);
  }
  consume(remove_args);
}

//////////////////////////////////////////////////////////////////////////////

void ConfigArgs::pass_filter(const ConfigKey& filtkey)
{
  std::vector<ConfigKey> remove_args;
  for (ConfigArgs::iterator it = begin(); it != end(); ++it)
  {
    const ConfigKey& this_key = it->first;
    if ( ! StringOps::startsWith(this_key,filtkey) )
      remove_args.push_back(this_key);
  }
  consume(remove_args);
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Config

} // namespace COOLFluiD

