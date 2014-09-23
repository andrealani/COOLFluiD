// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

/*
 * SyslogAppender.hh
 *
 * Copyright 2001, LifeLine Networks BV (www.lifeline.nl). All rights reserved.
 * Copyright 2001, Walter Stroebel. All rights reserved.
 *
 * See the COPYING file for the terms of usage and distribution.
 */

#ifndef CF_LOGCPP_REMOTESYSLOGAPPENDER_HH
#define CF_LOGCPP_REMOTESYSLOGAPPENDER_HH

#include <logcpp/Portability.hh>
#include <string>
#include <stdarg.h>
#include <logcpp/LayoutAppender.hh>
#include <logcpp/Priority.hh>

#ifdef LOGCPP_HAVE_SYSLOG
  #include <syslog.h>
#else
/// from syslog.h
typedef enum {
    LOG_EMERG   = 0,       ///< system is unusable
    LOG_ALERT   = 1,       ///< action must be taken immediately
    LOG_CRIT    = 2,       ///< critical conditions
    LOG_ERR     = 3,       ///< error conditions
    LOG_WARNING = 4,       ///< warning conditions
    LOG_NOTICE  = 5,       ///< normal but significant condition
    LOG_INFO    = 6,       ///< informational
    LOG_DEBUG   = 7,       ///< debug-level messages
} SyslogLevel;

typedef enum {
    LOG_KERN     = (0<<3),  ///< kernel messages
    LOG_USER     = (1<<3),  ///< random user-level messages
    LOG_MAIL     = (2<<3),  ///< mail system
    LOG_DAEMON   = (3<<3),  ///< system daemons
    LOG_AUTH     = (4<<3),  ///< security/authorization messages
    LOG_SYSLOG   = (5<<3),  ///< messages generated internally by syslogd
    LOG_LPR      = (6<<3),  ///< line printer subsystem
    LOG_NEWS     = (7<<3),  ///< network news subsystem
    LOG_UUCP     = (8<<3),  ///< UUCP subsystem
    LOG_CRON     = (9<<3),  ///< clock daemon
    LOG_AUTHPRIV = (10<<3), ///< security/authorization messages (private)
    LOG_FTP      = (11<<3), ///< ftp daemon

    /* other codes through 15 reserved for system use */
    LOG_LOCAL0   = (16<<3), ///< reserved for local use
    LOG_LOCAL1   = (17<<3), ///< reserved for local use
    LOG_LOCAL2   = (18<<3), ///< reserved for local use
    LOG_LOCAL3   = (19<<3), ///< reserved for local use
    LOG_LOCAL4   = (20<<3), ///< reserved for local use
    LOG_LOCAL5   = (21<<3), ///< reserved for local use
    LOG_LOCAL6   = (22<<3), ///< reserved for local use
    LOG_LOCAL7   = (23<<3), ///< reserved for local use
} SyslogFacility;
#endif

namespace logcpp {

    /**
     * RemoteSyslogAppender sends LoggingEvents to a remote syslog system.
     *
     * Also see: draft-ietf-syslog-syslog-12.txt
     **/
    class logcpp_API RemoteSyslogAppender : public LayoutAppender {
        public:

        /**
         * Translates a logcpp priority to a syslog priority
         * @param priority The logcpp priority.
         * @returns the syslog priority.
         **/
        static int toSyslogPriority(Priority::Value priority);

        /**
         * Instantiate a RemoteSyslogAppender with given name and name and
         * facility for syslog.
         * @param name The name of the Appender
         * @param syslogName The ident parameter in the openlog(3) call.
         * @param relayer The IP address or hostname of a standard syslog host.
         * @param facility The syslog facility to log to. Defaults to LOG_USER.
         * Value '-1' implies to use the default.
         * @param portNumber An alternative port number. Defaults to the
         * standard syslog port number (514).
         * Value '-1' implies to use the default.
         **/
        RemoteSyslogAppender(const std::string& name,
                             const std::string& syslogName,
                             const std::string& relayer,
                             int facility = LOG_USER,
                             int portNumber = 514);
        virtual ~RemoteSyslogAppender();

        /**
         * Closes and reopens the socket.
         **/
        virtual bool reopen();

        /**
         * Closes the socket
         **/
        virtual void close();

        protected:

        /**
         * Just creates the socket.
         **/
        virtual void open();

        /**
         * Sends a LoggingEvent to the remote syslog.
         * @param event the LoggingEvent to log.
         **/
        virtual void _append(const LoggingEvent& event);

        const std::string _syslogName;
        const std::string _relayer;
        int _facility;
        int _portNumber;
        int _socket;
        unsigned long _ipAddr;
        private:
        int _cludge;
    };
}

#endif // CF_LOGCPP_REMOTESYSLOGAPPENDER_HH
