// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Common_FilesystemException_hh
#define COOLFluiD_Common_FilesystemException_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/Exception.hh"
#include "Common/Common.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Common {

//////////////////////////////////////////////////////////////////////////////

/// This class represents the exception thrown
/// if an error occurs when accessing the filesystem.
/// It is preferable to using directly the boost::filesystem_error
/// exception. These boost exceptions should be intercepted and
/// recast into this.
/// @author Tiago Quintino
class Common_API FilesystemException : public Common::Exception {
public:

  /// Constructor
  FilesystemException ( const Common::CodeLocation& where, const std::string& what);

  /// Copy constructor
  FilesystemException ( const FilesystemException& e) throw ();

}; // end of class FilesystemException

//////////////////////////////////////////////////////////////////////////////

  } // namespace Common

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Common_FilesystemException_hh
