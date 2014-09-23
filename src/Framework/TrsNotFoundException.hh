// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_TrsNotFoundException_hh
#define COOLFluiD_Framework_TrsNotFoundException_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/Exception.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This class represents an Exception thrown when a Trs with
/// was not found in the available Trs's.
/// @author Tiago Quintino
class Framework_API TrsNotFoundException : public Common::Exception {
public:

  /// Constructor
  /// @see COOLFluiD::Exception()
  TrsNotFoundException (const Common::CodeLocation& where, const std::string& what) :
    Exception(where,what,"TrsNotFoundException")
  {
  }

  /// A copy constructor is necessary for exceptions, for the C++
  /// exception mechanism to work.
  TrsNotFoundException(const TrsNotFoundException& e)
    throw() : Exception(e)
  {
  }

}; // end of class TrsNotFoundException

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_TrsNotFoundException_hh
