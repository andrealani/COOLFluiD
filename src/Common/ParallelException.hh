// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Common_ParallelException_hh
#define COOLFluiD_Common_ParallelException_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/Exception.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
namespace Common {

//////////////////////////////////////////////////////////////////////////////

/// This exception is thrown in any place of the code which
/// by some conceptual impossibility should not be reached.
/// Typically on a switch-case construction where one of the choices
/// should be taken and the default never reached.
/// @author Tiago Quintino
class Common_API ParallelException : public Common::Exception {
public:

  /// Constructor
  ParallelException(const Common::CodeLocation& where, const std::string& what);

}; // end of class ParallelException

//////////////////////////////////////////////////////////////////////////////

  } // namespace Common
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Common_ParallelException_hh
