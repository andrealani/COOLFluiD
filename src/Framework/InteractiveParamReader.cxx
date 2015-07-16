// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Common/PEFunctions.hh"
#include "Common/FilesystemException.hh"

#include "Config/ConfigFileReader.hh"
#include "Config/ConfigFacility.hh"

#include "Environment/FileHandlerInput.hh"
#include "Environment/SingleBehaviorFactory.hh"
#include "Environment/DirPaths.hh"

#include "Framework/SubSystemStatus.hh"
#include "Framework/InteractiveParamReader.hh"
#include "Framework/MeshData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

void InteractiveParamReader::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< std::string >("FileName","Name of the file with the interactive parameters.");
  options.addConfigOption< CFuint >("readRate","The frequency of reading from the file with the interactive parameters.");
}

//////////////////////////////////////////////////////////////////////////////

InteractiveParamReader::InteractiveParamReader(const std::string& name)
  : ConfigObject(name)
{
  addConfigOptionsTo(this);

  m_filename = "";
  setParameter("FileName",&m_filename);

  m_read_rate = 1;
  setParameter("readRate",&m_read_rate);
}

//////////////////////////////////////////////////////////////////////////////

InteractiveParamReader::~InteractiveParamReader()
{
}

//////////////////////////////////////////////////////////////////////////////

void InteractiveParamReader::configure ( Config::ConfigArgs& args )
{
  ConfigObject::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void InteractiveParamReader::readFile()
{
  CFAUTOTRACE;

  // skip if no file is set
  if ( m_filename.empty() ) return;
  
  const CFuint iter = SubSystemStatusStack::getActive()->getNbIter();
  
  if ( !(iter % m_read_rate) ) {
    // to read in parallel: all the processes must stop and read from the master process's file
    // if this is a parallel simulation, only ONE process at a time reads the file
    readFileSequentially();
  }
}
    
//////////////////////////////////////////////////////////////////////////////

void InteractiveParamReader::readFileSequentially()
{
  CFAUTOTRACE;

  boost::filesystem::path fname;
  if (Common::StringOps::startsWith(m_filename,"."))
  {
    fname = boost::filesystem::path(m_filename);
  }
  else
  {
    fname = Environment::DirPaths::getInstance().getBaseDir() / boost::filesystem::path(m_filename);
  }

  // this would test if the file exists
  //   Common::SelfRegistPtr<Environment::FileHandlerInput> fhandle = Environment::SingleBehaviorFactory<Environment::FileHandlerInput>::getInstance().create();
  //   ifstream& inputFile = fhandle->open(fname);
  //   fhandle->close();

  Config::ConfigArgs args;
  Config::ConfigFileReader config_parser;
  config_parser.parse(fname, args);

  Config::ConfigFacility::getInstance().configureDynamicOptions(args);
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
