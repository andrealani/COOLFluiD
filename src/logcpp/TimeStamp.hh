// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

/*
 * TimeStamp.hh
 *
 * Copyright 2001, LifeLine Networks BV (www.lifeline.nl). All rights reserved.
 * Copyright 2001, Bastiaan Bakker. All rights reserved.
 *
 * See the COPYING file for the terms of usage and distribution.
 */

#ifndef CF_LOGCPP_TIMESTAMP_HH
#define CF_LOGCPP_TIMESTAMP_HH

#include <logcpp/Portability.hh>

namespace logcpp {

    /**
     * A simple TimeStamp abstraction
     **/
    class logcpp_API TimeStamp {
        public:
        /**
           Constructs a TimeStamp representing 'now'.
        **/
        TimeStamp();

        /**
           Constructs a TimeStamp representing the given offset since the
           epoch ( 00:00:00 1970/1/1 UTC).
        **/
        TimeStamp(unsigned int seconds, unsigned int microSeconds = 0);

        /**
           Returns the 'seconds' part of the TimeStamp.
        **/
        inline int getSeconds() const {
            return _seconds;
        };

        /** 
           Returns the 'subseconds' part of the TimeStamp in milliseconds,
           getMilliSeconds() == getMicroSeconds() / 1000. 
        **/
        inline int getMilliSeconds() const {
            return _microSeconds / 1000;
        };

        /**
           Returns the subsecond part of the TimeStamp in microseconds.
           The actual precision of this value depends on the platform and
           may be in the order of milliseconds rather than microseconds.
         **/
        inline int getMicroSeconds() const {
            return _microSeconds;
        };

        /**
           Returns a TimeStamp representing the time at which the application
           started.
        **/
        static inline const TimeStamp& getStartTime() {
            return _startStamp;
        };

        protected:
        static TimeStamp _startStamp;

        int _seconds;
        int _microSeconds;
    };
}

#endif // CF_LOGCPP_TIMESTAMP_HH

