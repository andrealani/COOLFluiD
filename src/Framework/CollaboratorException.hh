// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_CollaboratorException_hh
#define COOLFluiD_Framework_CollaboratorException_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/Exception.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This class represents the exception thrown
/// when trying to access a collaborator Method
/// @author Tiago Quintino
class Framework_API CollaboratorException : public Common::Exception {
public:

  /// Constructor
  /// @parameter what  is the name of the file that has been requested, but
  ///                  cannot actually be opened
  /// @see Exception()
  CollaboratorException(const Common::CodeLocation& where, const std::string& what) :
    Exception(where,what,"CollaboratorException")
  {
  }

  /// A copy constructor is necessary for exceptions, for the C++
  /// exception mechanism to work.
  CollaboratorException(const CollaboratorException& e)
    throw() : Exception(e)
  {
  }

}; // end of class CollaboratorException

//////////////////////////////////////////////////////////////////////////////

} // namespace Framework
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_CollaboratorException_hh
