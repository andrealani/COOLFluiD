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

/// This class represents exception that arise
/// when the searched Enviromental Variable does not exist.
/// @author Tiago Quintino
class Config_API BadMatchException : public Common::Exception {
public:

  /// Constructor
  BadMatchException(const Common::CodeLocation& where, const std::string& what) :
    Exception(where, what,"BadMatchException") {}

  /// A copy constructor is necessary for exceptions, for the C++
  /// exception mechanism to work.
  BadMatchException(const BadMatchException& e) throw() : Exception(e) {}

}; // end of class BadMatchException

//////////////////////////////////////////////////////////////////////////////

  } // namespace Config

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Config_BadMatchException_hh
