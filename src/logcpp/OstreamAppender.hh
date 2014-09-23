// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

/*
 * OstreamAppender.hh
 *
 * Copyright 2000, LifeLine Networks BV (www.lifeline.nl). All rights reserved.
 * Copyright 2000, Bastiaan Bakker. All rights reserved.
 *
 * See the COPYING file for the terms of usage and distribution.
 */

#ifndef CF_LOGCPP_OSTREAMAPPENDER_HH
#define CF_LOGCPP_OSTREAMAPPENDER_HH

#include <logcpp/Portability.hh>
#include <string>
#include <iostream>
#include <logcpp/LayoutAppender.hh>

namespace logcpp {

    /**
     * OstreamAppender appends LoggingEvents to ostreams.
     **/
    class logcpp_API OstreamAppender : public LayoutAppender {
        public:
        OstreamAppender(const std::string& name, std::ostream* stream);
        virtual ~OstreamAppender();
        
        virtual bool reopen();
        virtual void close();

        protected:
        virtual void _append(const LoggingEvent& event);

        std::ostream* _stream;
    };
}

#endif // CF_LOGCPP_OSTREAMAPPENDER_HH
