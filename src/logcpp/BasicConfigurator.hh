// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

/*
 * BasicConfigurator.hh
 *
 * Copyright 2002, Log4cpp Project. All rights reserved.
 *
 * See the COPYING file for the terms of usage and distribution.
 */
#ifndef CF_LOGCPP_BASICCONFIGURATOR_HH
#define CF_LOGCPP_BASICCONFIGURATOR_HH

#include <logcpp/Portability.hh>

namespace logcpp {

    /**
       This class implements a trivial default configuration for logcpp:
       it adds a FileAppender that logs to stdout and uses a BasicLayout to
       the root Category.
       @since 0.3.2
     **/
    class logcpp_API BasicConfigurator {
    public:

        /**
           Performs a minimal configuration of logcpp.  
         **/
        static void configure();
 };
}

#endif
