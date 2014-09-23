// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Common/ParserException.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace Common {

//////////////////////////////////////////////////////////////////////////////

ParserException::ParserException ( const Common::CodeLocation& where, const std::string& what) 
: Common::Exception(where, what, "ParserException") {}

ParserException::ParserException ( const ParserException& e) throw () 
: Exception(e) 
{}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Common
} // namespace COOLFluiD
