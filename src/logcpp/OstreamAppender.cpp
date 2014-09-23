/*
 * OstreamAppender.cpp
 *
 * Copyright 2000, LifeLine Networks BV (www.lifeline.nl). All rights reserved.
 * Copyright 2000, Bastiaan Bakker. All rights reserved.
 *
 * See the COPYING file for the terms of usage and distribution.
 */

#include "PortabilityImpl.hh"
#ifdef CF_HAVE_UNISTD_H
#    include <unistd.h>
#endif
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <logcpp/OstreamAppender.hh>

namespace logcpp {

    OstreamAppender::OstreamAppender(const std::string& name, std::ostream* stream) : 
        LayoutAppender(name),
        _stream(stream) {
    }
    
    OstreamAppender::~OstreamAppender() {
        close();
    }

    void OstreamAppender::close() {
        // empty
    }

    void OstreamAppender::_append(const LoggingEvent& event) {
        (*_stream) << _getLayout().format(event);
        if (!_stream->good()) {
            // XXX help! help!
        }
    }

    bool OstreamAppender::reopen() {
        return true;
    }      
}
