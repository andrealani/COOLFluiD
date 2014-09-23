// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

/*
 * RollingFileAppender.hh
 *
 * See the COPYING file for the terms of usage and distribution.
 */

#ifndef CF_LOGCPP_ROLLINGFILEAPPENDER_HH
#define CF_LOGCPP_ROLLINGFILEAPPENDER_HH

#include <logcpp/Portability.hh>
#include <logcpp/FileAppender.hh>
#include <string>
#include <stdarg.h>

namespace logcpp {

    /**
       RollingFileAppender is a FileAppender that rolls over the logfile once
       it has reached a certain size limit.
       @since 0.3.1
    **/
    class logcpp_API RollingFileAppender : public FileAppender {
        public:
        RollingFileAppender(const std::string& name, 
                            const std::string& fileName,
                            size_t maxFileSize = 10*1024*1024, 
                            unsigned int maxBackupIndex = 1,
                            bool append = true,
                            mode_t mode = 00644);

        virtual void setMaxBackupIndex(unsigned int maxBackups);
        virtual unsigned int getMaxBackupIndex() const;
        virtual void setMaximumFileSize(size_t maxFileSize);
        virtual size_t getMaxFileSize() const;

        virtual void rollOver();

        protected:
        virtual void _append(const LoggingEvent& event);

        unsigned int _maxBackupIndex;
        size_t _maxFileSize;
    };
}

#endif // CF_LOGCPP_ROLLINGFILEAPPENDER_HH
