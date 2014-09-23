// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_CallWithNoEffectException_hh
#define COOLFluiD_Framework_CallWithNoEffectException_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/Exception.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
namespace Framework {

//////////////////////////////////////////////////////////////////////////////

  /// This exception is thrown by {@see Connection} if a non-pure virtual
  /// method in the base class is called and it does nothing else than
  /// throwing this.
  /// @author Andrea Lani
class Framework_API CallWithNoEffectException : public Common::Exception {
public:

  /// Constructor
  /// @param what string holding the message to print if the exception
  ///             is caught
  /// @see Exception()
  CallWithNoEffectException (const Common::CodeLocation& where, const std::string& what)
  : Exception(where,what,"CallWithNoEffect")
  {
  }

  /// A copy constructor is necessary for exceptions, for the C++
  /// exception mechanism to work.
  CallWithNoEffectException(const CallWithNoEffectException& e) throw() :
    Exception(e)
  {
  }

}; // end of class CallWithNoEffectException

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_CallWithNoEffectException_hh
