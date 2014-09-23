// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_ConsistencyException_hh
#define COOLFluiD_Framework_ConsistencyException_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/Exception.hh"
#include "Framework/Framework.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This class represents an Exception thrown when a consistency check has failed
/// @author Tiago Quintino
class Framework_API ConsistencyException : public Common::Exception {
public:

  /// Constructor
  ConsistencyException ( const Common::CodeLocation& where, const std::string& what) :
    Common::Exception(where, what, "ConsistencyException") {}

  /// Copy constructor
  ConsistencyException ( const ConsistencyException& e) throw () : Exception(e) {}

}; // end of class ConsistencyException

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_ConsistencyException_hh
