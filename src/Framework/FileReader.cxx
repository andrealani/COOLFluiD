// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include <limits>
#include "Common/CFLog.hh"
#include "Framework/FileReader.hh"
#include "Framework/BadFormatException.hh"
#include "Environment/FileHandlerInput.hh"
#include "Environment/SingleBehaviorFactory.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

FileReader::FileReader() :
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

  const CFuint MAXLINESINFILE = std::numeric_limits<CFuint>::max()-1;

  m_readCount = 0;

  do {
    Common::SelfRegistPtr<Environment::FileHandlerInput> fhandle =
      Environment::SingleBehaviorFactory<Environment::FileHandlerInput>::getInstance().create();
    ifstream& file = fhandle->open(filepath);

    m_readAgain = false;
    bool keepOnReading = true;
    CFuint linesRead = 0;

    do {
      keepOnReading = readString(file);

      if (++linesRead > MAXLINESINFILE)
        throw BadFormatException (FromHere(),"File too long, probably misformatted."
				  "\nSee void FileReader::readFromFile for more info\n");

    } while (keepOnReading);

    fhandle->close();

    m_readCount++;

  } while (m_readAgain);

  finish();
}

//////////////////////////////////////////////////////////////////////////////

void FileReader::finish()
{
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

