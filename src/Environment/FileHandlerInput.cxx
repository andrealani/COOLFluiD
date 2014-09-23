// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Environment/FileHandlerInput.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace Environment {

//////////////////////////////////////////////////////////////////////////////

FileHandlerInput::FileHandlerInput() : m_isopen(false)
{
}

//////////////////////////////////////////////////////////////////////////////

FileHandlerInput::~FileHandlerInput()
{
}

//////////////////////////////////////////////////////////////////////////////

std::ifstream& FileHandlerInput::open(const std::string& filepath,  std::ios_base::openmode mode)
{
   boost::filesystem::path p (filepath);
   return this->open(p,mode);
}

//////////////////////////////////////////////////////////////////////////////

std::fstream& FileHandlerInput::openBinary(const std::string& filepath,  std::ios_base::openmode mode)
{
  boost::filesystem::path p (filepath);
  return this->openBinary(p,mode);
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace Environment

} // namespace COOLFluiD
