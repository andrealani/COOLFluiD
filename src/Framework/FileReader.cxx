// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include <limits>
#include "Common/CFLog.hh"
#include "Common/FactoryRegistry.hh"
#include "Framework/FileReader.hh"
#include "Framework/BadFormatException.hh"
#include "Environment/FileHandlerInput.hh"
#include "Environment/SingleBehaviorFactory.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Environment;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

FileReader::FileReader() :
  m_fr(CFNULL),
  m_readAgain(false),
  m_readCount(0)
{
}

//////////////////////////////////////////////////////////////////////////////

FileReader::~FileReader()
{
}

//////////////////////////////////////////////////////////////////////////////

void FileReader::readFromFile(const boost::filesystem::path& filepath)
{
  CFAUTOTRACE;
  
  CFLog(VERBOSE, "FileReader::readFromFile() => start\n");
  
  const CFuint MAXLINESINFILE = std::numeric_limits<CFuint>::max()-1;

  m_readCount = 0;

  do {

    Common::SelfRegistPtr<Environment::FileHandlerInput>* fhandle =
      Environment::SingleBehaviorFactory<Environment::FileHandlerInput>::getInstance().createPtr();

    std::ifstream& file = (*fhandle)->open(filepath);
    
    m_readAgain = false;
    bool keepOnReading = true;
    CFuint linesRead = 0;

    do {
      keepOnReading = readString(file);

      if (++linesRead > MAXLINESINFILE)
        throw BadFormatException (FromHere(),"File too long, probably misformatted."
				  "\nSee void FileReader::readFromFile for more info\n");

    } while (keepOnReading);

    CFLog(VERBOSE, "FileReader::readFromFile() => 2\n");    
    (*fhandle)->close();
    CFLog(VERBOSE, "FileReader::readFromFile() => 3\n");

    m_readCount++;
    CFLog(VERBOSE, "FileReader::readFromFile() => 4\n");    

    delete fhandle;
  } while (m_readAgain);
 
  CFLog(VERBOSE, "FileReader::readFromFile() => 5\n");
  finish();
  CFLog(VERBOSE, "FileReader::readFromFile() => 6\n");

  CFLog(VERBOSE, "FileReader::readFromFile() => end\n");
}

//////////////////////////////////////////////////////////////////////////////

void FileReader::finish()
{
}

//////////////////////////////////////////////////////////////////////////////
 
void FileReader::setFactoryRegistry(Common::SafePtr<Common::FactoryRegistry> fr)
{
  m_fr = fr;
}

//////////////////////////////////////////////////////////////////////////////

 Common::SafePtr<Common::FactoryRegistry> FileReader::getFactoryRegistry() 
 {
#ifdef CF_HAVE_SINGLE_EXEC
   return m_fr;
#endif
 }

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

