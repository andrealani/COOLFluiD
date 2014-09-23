// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

/*
 * CategoryStream.hh
 *
 * Copyright 2001, LifeLine Networks BV (www.lifeline.nl). All rights reserved.
 * Copyright 2001, Bastiaan Bakker. All rights reserved.
 *
 * See the COPYING file for the terms of usage and distribution.
 */

#ifndef CF_LOGCPP_CATEGORYSTREAM_HH
#define CF_LOGCPP_CATEGORYSTREAM_HH

#include <logcpp/Portability.hh>
#include <logcpp/Priority.hh>
#include <sstream>

namespace logcpp {

    class Category; // forward declaration

    /**
     * This class enables streaming simple types and objects to a category.
     * Use category.errorStream(), etc. to obtain a CategoryStream class.
     **/
    class logcpp_API CategoryStream {
        public:

        /**
         * Enumeration of special 'Separators'. Currently only contains the
         * 'ENDLINE' separator, which separates two log messages.
         **/
        typedef enum {
            ENDLINE
        } Separator;

        /**
         * Construct a CategoryStream for given Category with given priority.
         * @param category The category this stream will send log messages to.
         * @param priority The priority the log messages will get or
         * Priority::NOTSET to silently discard any streamed in messages.
         **/
        CategoryStream(Category& category, Priority::Value priority);

        /**
         * Destructor for CategoryStream
         **/
        ~CategoryStream();

        /**
         * Returns the destination Category for this stream.
         * @returns The Category.
         **/
        inline Category& getCategory() const { return _category; };

        /**
         * Returns the priority for this stream.
         * @returns The priority.
         **/
        inline Priority::Value getPriority() const throw() {
            return _priority;
        };

        /**
         * Streams in a Separator. If the separator equals
         * CategoryStream::ENDLINE it sends the contents of the stream buffer
         * to the Category with set priority and empties the buffer.
         * @param separator The Separator
         * @returns A reference to itself.
         **/
        CategoryStream& operator<<(Separator separator);

        /**
         * Flush the contents of the stream buffer to the Category and
         * empties the buffer.
         **/
        void flush();

        /**
         * Stream in arbitrary types and objects.
         * @param t The value or object to stream in.
         * @returns A reference to itself.
         **/
        template<typename T> CategoryStream& operator<<(const T& t) {
            if (getPriority() != Priority::NOTSET) {
                if (!_buffer) {
                    if (!(_buffer = new std::ostringstream)) {
                        // XXX help help help
                    }
                }
                (*_buffer) << t;
            }
            return *this;
        }

        private:
        Category& _category;
        Priority::Value _priority;
        std::ostringstream* _buffer;
    };

}
#endif // CF_LOGCPP_CATEGORYSTREAM_HH
