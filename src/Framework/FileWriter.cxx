// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/FileWriter.hh"
#include "Common/CFLog.hh"
#include "Environment/FileHandlerOutput.hh"
#include "Environment/SingleBehaviorFactory.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

FileWriter::FileWriter()
{
}

//////////////////////////////////////////////////////////////////////////////

FileWriter::~FileWriter()
{
}

//////////////////////////////////////////////////////////////////////////////

void FileWriter::writeToFile(const boost::filesystem::path& filepath)
{
  CFAUTOTRACE;

  Common::SelfRegistPtr<Environment::FileHandlerOutput>* fhandle = 	
	Environment::SingleBehaviorFactory<Environment::FileHandlerOutput>::getInstance().createPtr();
  ofstream& file = (*fhandle)->open(filepath);

  writeToFileStream(file);

  (*fhandle)->close();
  delete fhandle;
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

