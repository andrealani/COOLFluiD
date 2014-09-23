// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Common_FileFormatException_hh
#define COOLFluiD_Common_FileFormatException_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/Exception.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Common {

//////////////////////////////////////////////////////////////////////////////

/// This class represents an Exception thrown when a certain
/// value is not found in a storage or container.
/// @author Andrea Lani
/// @author Tiago Quintino
class FileFormatException : public Common::Exception {
public:

  /// Constructor
  /// @param what is the value that has been requested,
  ///             but actually doesn't exist
  /// @see Exception()
  FileFormatException (const Common::CodeLocation& where, const std::string& what) :
    Common::Exception(where,what,"FileFormatException")
  {
  }

  /// A copy constructor is necessary for exceptions, for the C++
  /// exception mechanism to work.
  FileFormatException(const FileFormatException& e)
    throw() : Common::Exception(e)
  {
  }

}; // end of class FileFormatException

//////////////////////////////////////////////////////////////////////////////

  } // namespace Common

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Common_FileFormatException_hh
