// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Common_NoSuchStorageException_hh
#define COOLFluiD_Common_NoSuchStorageException_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/Exception.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Common {

//////////////////////////////////////////////////////////////////////////////

  /// This class represents an Exception thrown when a certain
  /// value is not found in a storage.
  /// @author Andrea Lani
  /// @author Tiago Quintino
class NoSuchStorageException : public Common::Exception {
public:

  /// Constructor
  NoSuchStorageException ( const Common::CodeLocation& where, const std::string& what) :
    Common::Exception(where, what, "NoSuchStorageException") {}

  /// Copy constructor
  NoSuchStorageException ( const NoSuchStorageException& e) throw () : Exception(e) {}

}; // end of class NoSuchStorageException

//////////////////////////////////////////////////////////////////////////////

  } // namespace Common

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Common_NoSuchStorageException_hh
