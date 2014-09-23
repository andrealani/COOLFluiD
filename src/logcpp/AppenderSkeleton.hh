// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

/*
 * AppenderSkeleton.hh
 *
 * Copyright 2001, LifeLine Networks BV (www.lifeline.nl). All rights reserved.
 * Copyright 2001, Bastiaan Bakker. All rights reserved.
 *
 * See the COPYING file for the terms of usage and distribution.
 */

#ifndef CF_LOGCPP_APPENDERSKELETON_HH
#define CF_LOGCPP_APPENDERSKELETON_HH

#include <logcpp/Portability.hh>
#include <logcpp/Appender.hh>
#include <logcpp/Filter.hh>

namespace logcpp {

    /**
     * AppenderSkeleton is a helper class, simplifying implementation of
     * Appenders: it already takes care of handling of Thresholds and 
     * Filters. 
     **/
    class logcpp_API AppenderSkeleton : public Appender {
        protected:
        /**
         * Constructor for AppenderSkeleton. Will only be used in 
         * getAppender() (and in derived classes of course).
         * @param name The name of this Appender.
         **/
        AppenderSkeleton(const std::string& name);

        public:
        /**
         * Destructor for AppenderSkeleton.
         **/
        virtual ~AppenderSkeleton();
        
        /**
         * Log in Appender specific way. 
         * @param event  The LoggingEvent to log. 
         **/
        virtual void doAppend(const LoggingEvent& event);

        /**
         * Reopens the output destination of this Appender, e.g. the logfile 
         * or TCP socket.
         * @returns false if an error occured during reopening, true otherwise.
         **/
        virtual bool reopen();

        /**
         * Release any resources allocated within the appender such as file
         * handles, network connections, etc.
         **/
        virtual void close() = 0;

        /**
         * Check if the appender uses a layout.
         * 
         * @returns true if the appender implementation requires a layout.
         **/
        virtual bool requiresLayout() const = 0;

        /**
         * Set the Layout for this appender.
         * @param layout The layout to use.
         **/
        virtual void setLayout(Layout* layout) = 0;

        /**
         * Set the threshold priority of this Appender. The Appender will not
         * appender LoggingEvents with a priority lower than the threshold.
         * Use Priority::NOTSET to disable threshold checking.
         * @param priority The priority to set.
         **/
        virtual void setThreshold(Priority::Value priority);

        /**
         * Get the threshold priority of this Appender.
         * @returns the threshold
         **/
        virtual Priority::Value getThreshold();

        /**
         * Set a Filter for this appender. 
         **/
        virtual void setFilter(Filter* filter);

        /**
         * Get the Filter for this appender.
         * @returns the filter, or NULL if no filter has been set.
         **/
        virtual Filter* getFilter();

        protected:
        /**
         * Log in Appender specific way. Subclasses of Appender should
         * implement this method to perform actual logging.
         * @param event  The LoggingEvent to log. 
         **/
        virtual void _append(const LoggingEvent& event) = 0;

        
        private:
        Priority::Value _threshold;
        Filter* _filter;
    };
}

#endif // CF_LOGCPP_APPENDERSKELETON_HH
