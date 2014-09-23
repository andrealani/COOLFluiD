// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Common_NoSuchValueException_hh
#define COOLFluiD_Common_NoSuchValueException_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/Exception.hh"

#include "Common/Common.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Common {

//////////////////////////////////////////////////////////////////////////////

/// This class represents an Exception thrown when a certain
/// value is not found in a storage or container.
/// @author Andrea Lani
/// @author Tiago Quintino
class Common_API NoSuchValueException : public Common::Exception {
public:

  /// Constructor
  NoSuchValueException ( const Common::CodeLocation& where, const std::string& what);

  /// Copy constructor
  NoSuchValueException ( const NoSuchValueException& e) throw ();

}; // end of class NoSuchValueException

//////////////////////////////////////////////////////////////////////////////

  } // namespace Common

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Common_NoSuchValueException_hh
