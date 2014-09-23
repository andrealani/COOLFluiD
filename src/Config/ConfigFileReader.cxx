// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include <fstream>

#include "Common/StringOps.hh"
#include "Common/FilesystemException.hh"
#include "Config/ConfigArgs.hh"
#include "Config/ConfigFileReader.hh"
#include "Config/XMLConfigFileReader.hh"
#include "Config/NestedConfigFileReader.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;

namespace COOLFluiD {

  namespace Config {

//////////////////////////////////////////////////////////////////////////////

ConfigFileReader::ConfigFileReader()
{
}

//////////////////////////////////////////////////////////////////////////////

ConfigFileReader::~ConfigFileReader()
{
}

//////////////////////////////////////////////////////////////////////////////

void ConfigFileReader::parse(const boost::filesystem::path& name, ConfigArgs& args)
{
    // check the file type
    std::ifstream iss(name.string().c_str());
    if(!iss) throw Common::FilesystemException (FromHere(),"Could not open file: " + name.string());
    std::string line;
    getline(iss, line);
    bool is_xml = Common::StringOps::startsWith(line,"<?xml");
    iss.close();

    if (is_xml)
    {
      XMLConfigFileReader xmlConfig;
      xmlConfig.parse(name.string(), args);
    }
    else // not xml
    {
      NestedConfigFileReader configFile;
      configFile.parse(name.string(), args);
    }
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Config

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
