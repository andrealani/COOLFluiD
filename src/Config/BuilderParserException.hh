// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Config_BadMatchException_hh
#define COOLFluiD_Config_BadMatchException_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/Exception.hh"
#include "Config/Config.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

 namespace Config {

//////////////////////////////////////////////////////////////////////////////

  /// This exception is thrown by the builder parser on any error
  /// @author Quentin Gasper
  class Config_API BuilderParserException : public Common::Exception {
   public:

  /// Constructor
  BuilderParserException(const Common::CodeLocation& where, const std::string& what);

  /// A copy constructor is necessary for exceptions, for the C++
  /// exception mechanism to work.
  BuilderParserException(const BuilderParserException& e);
  
  }; // class BuilderParserException

 //////////////////////////////////////////////////////////////////////////////

 } // namespace Config
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Config_BadMatchException_hh
