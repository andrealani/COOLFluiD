// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

/*
 * LoggingEvent.hh
 *
 * Copyright 2000, LifeLine Networks BV (www.lifeline.nl). All rights reserved.
 * Copyright 2000, Bastiaan Bakker. All rights reserved.
 *
 * See the COPYING file for the terms of usage and distribution.
 */

#ifndef CF_LOGCPP_LOGGINGEVENT_HH
#define CF_LOGCPP_LOGGINGEVENT_HH

#include <logcpp/Portability.hh>
#include <string>

#include <logcpp/Priority.hh>
#include <logcpp/TimeStamp.hh>

/**
 * The top level namespace for all 'Log for C++' types and classes.
 **/
namespace logcpp {

    /**
     * The internal representation of logging events. When a affirmative
     * logging decision is made a <code>LoggingEvent</code> instance is
     * created. This instance is passed around the different logcpp
     * components.
     *
     * <p>This class is of concern to those wishing to extend logcpp. 
     **/
    struct logcpp_API LoggingEvent {
    public:
        /**
         * Instantiate a LoggingEvent from the supplied parameters.
         *
         * <p>Except <code>timeStamp</code> all the other fields of
         * <code>LoggingEvent</code> are filled when actually needed.
         * <p>
         * @param category The category of this event.
         * @param message  The message of this event.
         * @param ndc The nested diagnostic context of this event. 
         * @param priority The priority of this event.
         **/
        LoggingEvent(const std::string& category, const std::string& message, 
                     const std::string& ndc, Priority::Value priority);


        /** The category name. */
        const std::string categoryName;

        /** The application supplied message of logging event. */
        const std::string message;

        /** The nested diagnostic context (NDC) of logging event. */
        const std::string ndc;

        /** Priority of logging event. */
        Priority::Value priority;

        /** The name of thread in which this logging event was generated,
            e.g. the PID. 
        */
        const std::string threadName;

        /** The number of seconds elapsed since the epoch 
            (1/1/1970 00:00:00 UTC) until logging event was created. */
        TimeStamp timeStamp;
    };
}

#endif // CF_LOGCPP_LOGGINGEVENT_HH

