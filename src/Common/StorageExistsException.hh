// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Common_StorageExistsException_hh
#define COOLFluiD_Common_StorageExistsException_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/Exception.hh"

#include "Common/Common.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Common {

//////////////////////////////////////////////////////////////////////////////

/// This class represents an Exception thrown when a certain
/// value already exists in a storage.
/// @author Andrea Lani
/// @author Tiago Quintino
class Common_API StorageExistsException : public Common::Exception {
public:

  /// Constructor
  StorageExistsException ( const Common::CodeLocation& where, const std::string& what) :
    Common::Exception(where, what, "StorageExistsException") {}

  /// Copy constructor
  StorageExistsException ( const StorageExistsException& e) throw () : Exception(e) {}

}; // end of class StorageExistsException

//////////////////////////////////////////////////////////////////////////////

  } // namespace Common

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Common_StorageExistsException_hh
