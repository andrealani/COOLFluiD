// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Common_CFLog_hh
#define COOLFluiD_Common_CFLog_hh

//////////////////////////////////////////////////////////////////////////////

#if defined(CF_HAVE_CRAYSTATIC) || !defined(CF_HAVE_LOG4CPP)
#include <iostream>
#endif

#include "logcpp/Category.hh"

#include "Common/COOLFluiD.hh"

#if (defined(CF_HAVE_IBMSTATIC) || defined(CF_HAVE_CRAYSTATIC) || !defined(CF_HAVE_LOG4CPP))  && defined(CF_HAVE_MPI)
#include <mpi.h>
#endif

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

//////////////////////////////////////////////////////////////////////////////

/// Output levels for the CFLog facility
/// These values are relative to the logcpp levels, so keep them this way.
/// @author Tiago Quintino
enum CFLogLevel { EMERG  = 0,
                  FATAL  = 0,
                  ALERT  = 100,
                  CRIT   = 200,
                  ERROR  = 300, 
                  WARN   = 400,
                  NOTICE = 500,
                  INFO   = 600,
		  VERBOSE = 650,
                  DEBUG_MIN = 700,
		  DEBUG_MED = 725,
		  DEBUG_MAX = 750,
                  NOTSET = 800
};

//////////////////////////////////////////////////////////////////////////////

class Common_API CFLogger : public Common::NonCopyable<CFLogger> {
public: // functions

  /// @return the instance of this singleton
  static CFLogger& getInstance();

  /// Gets the main stream logger
  /// @returns the main logger
  logcpp::Category& getMainLogger();

  /// Gets the trace stream logger
  /// @returns the logger usefull for tracing messages
  logcpp::Category& getTraceLogger();

  /// Sets if the trace output should go to the stdout
  void setTraceToStdOut(bool traceToStdOut);
  
  /// Set the main logger level of priority
  void setMainLoggerLevel(CFuint level);
  
  /// Get the main logger level of priority
  CFuint getMainLoggerLevel() {return getMainLogger().getPriority();}

private: // methods

  /// Default contructor
  CFLogger();
  /// Default destructor
  ~CFLogger();

};

//////////////////////////////////////////////////////////////////////////////
// Logging macros
//////////////////////////////////////////////////////////////////////////////

/// these are always defined
#if !defined(CF_HAVE_CRAYSTATIC) || defined(CF_HAVE_LOG4CPP)  
#define CFout    CFLogger::getInstance().getMainLogger().noticeStream()
#define CFerr    CFLogger::getInstance().getMainLogger().errorStream()
#define CFlog    CFLogger::getInstance().getMainLogger().infoStream()
#define CFtrace  CFLogger::getInstance().getTraceLogger().debugStream()
#define CFendl   logcpp::CategoryStream::ENDLINE
#else
#define CFout    std::cout
#define CFerr    std::cout
#define CFlog    std::cout
#define CFtrace  std::cout
#define CFendl   std::endl
#endif

//////////////////////////////////////////////////////////////////////////////

// bypass default CFLog when compiling with IBM compiler with static linking
#if (defined(CF_HAVE_IBMSTATIC) || defined(CF_HAVE_CRAYSTATIC)) || !defined(CF_HAVE_LOG4CPP) && defined(CF_HAVE_MPI)
static int getCPURank() 
{
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  return rank;
}
#ifdef CF_HAVE_IBMSTATIC
#define CFLog(n,x) if (n <= CFLogger::getInstance().getMainLoggerLevel() && getCPURank() == 0) CFLogger::getInstance().getMainLogger() << n << x 
#endif

#if defined(CF_HAVE_CRAYSTATIC) || !defined(CF_HAVE_LOG4CPP)
#define CFLog(n,x) if (n <= CFLogger::getInstance().getMainLoggerLevel() && getCPURank() == 0) std::cout << x
#endif

#endif

#if !defined(CF_HAVE_IBMSTATIC) && !defined(CF_HAVE_CRAYSTATIC) && defined(CF_HAVE_LOG4CPP)
#ifndef CF_NO_LOG
  #define CFLog(n,x) if (n <= CFLogger::getInstance().getMainLoggerLevel()) CFLogger::getInstance().getMainLogger() << n << x
 #else
  #define CFLog(n,x)
#endif
#endif

//////////////////////////////////////////////////////////////////////////////

#ifndef CF_NO_DEBUG_LOG

  #define CFLogDebug(x)    CFLog(DEBUG_MIN,x)
  #define CFLogDebugMax(x) CFLog(DEBUG_MAX,x)
  #define CFLogDebugMed(x) CFLog(DEBUG_MED,x)
  #define CFLogDebugMin(x) CFLog(DEBUG_MIN,x)

#else // CF_NO_DEBUG_LOG

  #define CFLogDebug(x)
  #define CFLogDebugMax(x)
  #define CFLogDebugMed(x)
  #define CFLogDebugMin(x)

#endif // CF_NO_DEBUG_LOG

#define CFLogInfo(x)   CFLog(INFO,x)
#define CFLogNotice(x) CFLog(NOTICE,x)
#define CFLogWarn(x)   CFLog(WARN,x)
#define CFLogError(x)  CFerr << x ; CFerr.flush()

//////////////////////////////////////////////////////////////////////////////
// Tracing macros
//////////////////////////////////////////////////////////////////////////////

/// Class to help trace functions which is easy to use and catches all function exits (return,throw,...)
/// @author Dries Kimpe
class Common_API AutoTracer {
public:
    /// constructor
    AutoTracer (const char * Function, const char * File, int Line) : m_Function(Function), m_File(File), m_Line(Line)
    {
      CFtrace << "### BEGIN ### " << m_Function << " : " << m_File << " : " << m_Line << "\n" ; CFtrace.flush();
    }

    /// destructor
    ~AutoTracer () { CFtrace << "### END   ### " << m_Function << " : " << m_File << " : " << m_Line << "\n" ; CFtrace.flush(); }

protected: // data
    const char * m_Function;
    const char * m_File;
    int m_Line;
};

//////////////////////////////////////////////////////////////////////////////

#ifndef CF_NO_TRACE

  #define CFAUTOTRACE_PASTE_(a,b) a##b
  #define CFAUTOTRACE_PASTE(a,b) CFAUTOTRACE_PASTE_(a,b)
  #define CFAUTOTRACE ::COOLFluiD::AutoTracer CFAUTOTRACE_PASTE(AutoTrace_Uniq,__LINE__) (__FUNCTION__,__FILE__,__LINE__)

  #define CFTRACEBEGIN CFtrace << "### BEGIN ### " << __FUNCTION__ << " : " << __FILE__ << " : " << __LINE__ << "\n" ; CFtrace.flush();
  #define CFTRACEEND   CFtrace << "### END ##### " << __FUNCTION__ << " : " << __FILE__ << " : " << __LINE__ << "\n" ; CFtrace.flush();

#else // CF_NO_TRACE

#define CFAUTOTRACE

  #define CFTRACEBEGIN
  #define CFTRACEEND

#endif // CF_NO_TRACE

//////////////////////////////////////////////////////////////////////////////
// Debugging macros
//////////////////////////////////////////////////////////////////////////////

#ifndef CF_NO_DEBUG_MACROS

/// Definition of a macro for placing a debug point in the code
#define CF_DEBUG_POINT  CFerr << "DEBUG : " << __FILE__ << " : " << __LINE__ << " : " << __FUNCTION__ << "\n" ; CFerr.flush()
/// Definition of a macro for outputing an object that implements the output stream operator
#define CF_DEBUG_OBJ(x) CFerr << "DEBUG : OBJECT " << #x << " -> " << x << " : " << __FILE__ << " : " << __LINE__ << " : " << __FUNCTION__ << "\n" ; CFerr.flush()
/// Definition of a macro for outputing a debug string in the code
#define CF_DEBUG_STR(x) CFerr << "DEBUG : STRING : " << x << " : " << __FILE__ << " : " << __LINE__ << " : " << __FUNCTION__ << "\n" ; CFerr.flush()
/// Definition of a macro for debug exit
#define CF_DEBUG_EXIT   CFerr << "DEBUG : EXIT "  << __FILE__ << " : " << __LINE__ << " : " << __FUNCTION__ << "\n" ; CFerr.flush() ; exit(0)
/// Definition of a macro for debug abort
#define CF_DEBUG_ABORT  CFerr << "DEBUG : ABORT " << __FILE__ << " : " << __LINE__ << " : " << __FUNCTION__ << "\n" ; CFerr.flush() ; abort()

#else // CF_NO_DEBUG_MACROS

#define CF_DEBUG_POINT
#define CF_DEBUG_OBJ(x)
#define CF_DEBUG_STR(x)
#define CF_DEBUG_EXIT
#define CF_DEBUG_ABORT

#endif // CF_NO_DEBUG_MACROS

//////////////////////////////////////////////////////////////////////////////

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Common_CFLog_hh
