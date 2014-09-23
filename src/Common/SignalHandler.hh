// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Common_SignalHandler_hh
#define COOLFluiD_Common_SignalHandler_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/NonCopyable.hh"
#include "Common/OwnedObject.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Common {

//////////////////////////////////////////////////////////////////////////////

/// This class is an interface for handlers of signals from the
/// opwrating system
/// @author Tiago Quintino
class Common_API SignalHandler : public Common::NonCopyable<SignalHandler> {

public: // methods

/// signal function type
typedef void (*sighandler_t)(int);

  /// Constructor is private to allow only the friend classes to build it
  SignalHandler();

  /// Default destructor is private to allow only the friend classes to destroy it
  virtual ~SignalHandler();

  /// Regists the signal handlers that will be handled by this class
  virtual void registSignalHandlers() = 0;

  /// Gets the Class name
  static std::string getClassName() { return "SignalHandler"; }

protected: // methods

private: // data

}; // end of class SignalHandler

//////////////////////////////////////////////////////////////////////////////

  } // namespace Common

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOFluiD_Common_SignalHandler_hh
