// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_BadFormatException_hh
#define COOLFluiD_Framework_BadFormatException_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/Exception.hh"
#include "Framework/Framework.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This class represents the exception thrown if a certain path
/// to a directory does not exist
/// @author Tiago Quintino
class Framework_API BadFormatException : public Common::Exception {
public:

  /// Constructor
  BadFormatException(const Common::CodeLocation& where, const std::string& what) :
    Exception(where, what,"BadFormatException") {}

  /// A copy constructor is necessary for exceptions, for the C++
  /// exception mechanism to work.
  BadFormatException(const BadFormatException& e) throw() : Exception(e) {}

}; // end of class BadFormatException

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_BadFormatException_hh
