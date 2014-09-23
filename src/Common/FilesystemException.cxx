// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Common/FilesystemException.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Common {

//////////////////////////////////////////////////////////////////////////////

FilesystemException::FilesystemException ( const Common::CodeLocation& where, const std::string& what)
: Common::Exception(where, what, "FilesystemException") 
{}


//////////////////////////////////////////////////////////////////////////////

FilesystemException::FilesystemException ( const FilesystemException& e) throw ()
: Exception(e) 
{}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Common

} // namespace COOLFluiD
