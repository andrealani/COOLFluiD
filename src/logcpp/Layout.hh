// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

/*
 * Layout.hh
 *
 * Copyright 2000, LifeLine Networks BV (www.lifeline.nl). All rights reserved.
 * Copyright 2000, Bastiaan Bakker. All rights reserved.
 *
 * See the COPYING file for the terms of usage and distribution.
 */

#ifndef CF_LOGCPP_LAYOUT_HH
#define CF_LOGCPP_LAYOUT_HH

#include <logcpp/Portability.hh>
#include <logcpp/LoggingEvent.hh>
#include <string>

namespace logcpp {

/**
 * Extend this abstract class to create your own log layout format.
 **/
    class logcpp_API Layout {
        public:
        /**
         * Destructor for Layout.
         **/
        virtual ~Layout() { };

        /**
         * Formats the LoggingEvent data to a string that appenders can log.
         * Implement this method to create your own layout format.
         * @param event The LoggingEvent.
         * @returns an appendable string.
         **/
        virtual std::string format(const LoggingEvent& event) = 0;
    };        
}

#endif // CF_LOGCPP_LAYOUT_HH
