// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

/*
 * StringQueueAppender.hh
 *
 * Copyright 2000, LifeLine Networks BV (www.lifeline.nl). All rights reserved.
 * Copyright 2000, Bastiaan Bakker. All rights reserved.
 *
 * See the COPYING file for the terms of usage and distribution.
 */

#ifndef CF_LOGCPP_STRINGQUEUEAPPENDER_HH
#define CF_LOGCPP_STRINGQUEUEAPPENDER_HH

#include <logcpp/Portability.hh>
#include <string>
#include <queue>
#include <logcpp/LayoutAppender.hh>

namespace logcpp {

    /**
     * This class puts log messages in an in-memory queue. Its primary use 
     * is in test cases, but it may be useful elsewhere as well.
     *
     * @since 0.2.4
     **/
    class logcpp_API StringQueueAppender : public LayoutAppender {
        public:

        StringQueueAppender(const std::string& name);
        virtual ~StringQueueAppender();
        
        virtual bool reopen();
        virtual void close();

        /**
         * Return the current size of the message queue.
         * Shorthand for getQueue().size().
         * @returns the queue size
         **/
        virtual size_t queueSize() const;

        /**
         * Return the queue to which the Appends adds messages.
         * @returns the message queue
         **/
        virtual std::queue<std::string>& getQueue();

        /**
         * Return the queue to which the Appends adds messages.
         * @returns the message queue
         **/
        virtual const std::queue<std::string>& getQueue() const;

        /**
         * Pop the oldest log message from the front of the queue.
         * @returns the oldest log message
         **/
        virtual std::string popMessage();

        protected:
        
        /**
         * Appends the LoggingEvent to the queue.
         * @param event the LoggingEvent to layout and append to the queue.
         **/
        virtual void _append(const LoggingEvent& event);

        std::queue<std::string> _queue;
    };
}

#endif // CF_LOGCPP_STRINGQUEUEAPPENDER_HH
