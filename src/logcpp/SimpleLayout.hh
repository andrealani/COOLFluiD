// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

/*
 * SimpleLayout.hh
 *
 * Copyright 2000, LifeLine Networks BV (www.lifeline.nl). All rights reserved.
 * Copyright 2000, Bastiaan Bakker. All rights reserved.
 *
 * See the COPYING file for the terms of usage and distribution.
 */

#ifndef CF_LOGCPP_SIMPLELAYOUT_HH
#define CF_LOGCPP_SIMPLELAYOUT_HH

#include <logcpp/Portability.hh>
#include <logcpp/Layout.hh>

namespace logcpp {
    
    /**
     * BasicLayout is a simple fixed format Layout implementation. 
     **/
    class logcpp_API SimpleLayout : public Layout {
        public:
        SimpleLayout();
        virtual ~SimpleLayout();

        /**
         * Formats the LoggingEvent in SimpleLayout style:<br>
         * "priority - message"
         **/
        virtual std::string format(const LoggingEvent& event);
    };        
}

#endif // CF_LOGCPP_SIMPLELAYOUT_HH
