// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Common_OSystem_hh
#define COOLFluiD_Common_OSystem_hh

#include "Common/Exception.hh"
#include "Common/NonCopyable.hh"
#include "Common/SelfRegistPtr.hh"
#include "Common/SafePtr.hh"

#include "Common/CommonAPI.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Common {

  class ProcessInfo;
  class SignalHandler;
  class LibLoader;

//////////////////////////////////////////////////////////////////////////////

class Common_API OSystemException : public Common::Exception {

public:

  /// Constructor
  OSystemException ( const Common::CodeLocation& where, const std::string& what) :
    Common::Exception(where, what, "OSystemException") {}

  /// Copy constructor
  OSystemException ( const OSystemException& e) throw () : Exception(e) {}

}; // end class OSystemException

//////////////////////////////////////////////////////////////////////////////

/// Represents the operating system
/// @author Tiago Quintino
class Common_API OSystem : public Common::NonCopyable<OSystem> {

public: // methods

  /// @return the single object that represents the operating system
  static OSystem& getInstance();

  /// @return ProcessInfo object
  Common::SafePtr<Common::ProcessInfo> getProcessInfo();

  /// @return SignalHandler object
  Common::SafePtr<Common::SignalHandler> getSignalHandler();

  /// @return LibLoader object
  Common::SafePtr<Common::LibLoader> getLibLoader();

  /// Executes the command passed in the string
  /// @todo should return the output of the command but not yet implemented.
  void executeCommand (const std::string& call);

  /// Sleeps for a certain number of seconds
  /// @param seconds to sleep ( default is 1 )
  void sleep (const CFuint& seconds = 1);

private: // functions

  /// constructor
  OSystem ();
  /// destructor
  ~OSystem ();

private: // data

  /// memory usage object
  Common::ProcessInfo * m_process_info;
  /// signal handler object
  Common::SignalHandler * m_sig_handler;
  /// libloader object
  Common::LibLoader * m_lib_loader;

}; // class FileHandlerOutput

//////////////////////////////////////////////////////////////////////////////

  } // namespace Common

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Common_OSystem_hh
