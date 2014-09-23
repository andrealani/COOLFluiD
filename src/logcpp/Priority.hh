// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

/*
 * Priority.hh
 *
 * Copyright 2000, LifeLine Networks BV (www.lifeline.nl). All rights reserved.
 * Copyright 2000, Bastiaan Bakker. All rights reserved.
 *
 * See the COPYING file for the terms of usage and distribution.
 */

#ifndef CF_LOGCPP_PRIORITY_HH
#define CF_LOGCPP_PRIORITY_HH

#include <logcpp/Portability.hh>
#include <string>
#include <stdexcept>

/*
 * Optionally work around rudeness in windows.h on Win32.
 */
#ifdef ERROR
#ifdef LOGCPP_FIX_ERROR_COLLISION

namespace logcpp {
    static const int _tmpERRORValue = ERROR;
}

#undef ERROR
    static const int ERROR = logcpp::_tmpERRORValue;
#define ERROR ERROR

#else  // LOGCPP_FIX_ERROR_COLLISION
#error Naming collision for 'ERROR' detected. Please read the FAQ for a \
       workaround. 
#endif // LOGCPP_FIX_ERROR_COLLISION 

#endif // ERROR

/*
 * Other Win32 rudeness in EDK.h
 */
#ifdef DEBUG

#ifdef LOGCPP_FIX_ERROR_COLLISION

#undef DEBUG
#define DEBUG DEBUG

#else  // LOGCPP_FIX_ERROR_COLLISION
#error Naming collision for 'DEBUG' detected. Please read the FAQ for a \
       workaround. 
#endif // LOGCPP_FIX_ERROR_COLLISION 

#endif // DEBUG


namespace logcpp {

    /**
     * The Priority class provides importance levels with which one
     * can categorize log messages.
     **/
    class Priority {
        public:

        /**
         * Predefined Levels of Priorities. These correspond to the
         * priority levels used by syslog(3).
         **/
        typedef enum {EMERG  = 0, 
		      FATAL  = 0,
                      ALERT  = 100,
                      CRIT   = 200,
                      ERROR  = 300, 
                      WARN   = 400,
                      NOTICE = 500,
                      INFO   = 600,
                      DEBUG  = 700,
                      NOTSET = 800
        } PriorityLevel;

        /**
         * The type of Priority Values
         **/
        typedef int Value;

        /**
         * Returns the name of the given priority value.
         * Currently, if the value is not one of the PriorityLevel values,
         * the method returns the name of the largest priority smaller 
         * the given value.
         * @param priority the numeric value of the priority.
         * @returns a string representing the name of the priority.
         **/
        static const std::string& getPriorityName(int priority) throw();
	
	/**
	 * Returns the value of the given priority name. 
	 * This can be either one of EMERG ... NOTSET or a 
	 * decimal string representation of the value, e.g. '700' for DEBUG.
	 * @param priorityName the string containing the the of the priority
	 * @return the value corresponding with the priority name
	 * @throw std::invalid_argument if the priorityName does not 
	 * correspond with a known Priority name or a number
	 **/
        static Value getPriorityValue(const std::string& priorityName)
	throw(std::invalid_argument);
    };
}

#endif // CF_LOGCPP_PRIORITY_HH
