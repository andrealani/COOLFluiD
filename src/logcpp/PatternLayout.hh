// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

/*
 * PatternLayout.hh
 *
 * Copyright 2002, Bastiaan Bakker. All rights reserved.
 *
 * See the COPYING file for the terms of usage and distribution.
 */

#ifndef CF_LOGCPP_PATTERNLAYOUT_HH
#define CF_LOGCPP_PATTERNLAYOUT_HH

#include <logcpp/Portability.hh>
#include <logcpp/Layout.hh>
#include <logcpp/Configurator.hh>
#include <vector>
#include <sstream>

namespace logcpp {

    /**
     * PatternLayout is a simple fixed format Layout implementation.
     **/
    class logcpp_API PatternLayout : public Layout {
        public:
        /**
           The default conversion pattern
        **/
        static const char* DEFAULT_CONVERSION_PATTERN;

        /**
           A conversion pattern equivalent to the SimpleLayout.
        **/
        static const char* SIMPLE_CONVERSION_PATTERN;

        /**
           A conversion pattern equivalent to the BasicLayout.
        **/
        static const char* BASIC_CONVERSION_PATTERN;

        /**
           A conversion pattern equivalent to the TTCCLayout.
           Note: TTCCLayout is in log4j but not logcpp.
        **/
        static const char* TTCC_CONVERSION_PATTERN;

        PatternLayout();
        virtual ~PatternLayout();

        // NOTE: All double percentage signs ('%%') followed by a character
        //       in the following comments should actually be a single char.
        //       The doubles are included so that doxygen will print them correctly.
        /**
         * Formats the LoggingEvent in the style set by
		 * the setConversionPattern call. By default, set
		 * to "%%m%%n"
         **/
        virtual std::string format(const LoggingEvent& event);

        /**
         * Sets the format of log lines handled by this
         * PatternLayout. By default, set to "%%m%%n".<br>
         * Format characters are as follows:<br>
         * <li><b>%%%%</b> - a single percent sign</li>
         * <li><b>%%c</b> - the category</li>
         * <li><b>%%d</b> - the date\n
         *  Date format: The date format character may be followed by a date format
         *  specifier enclosed between braces. For example, %%d{%%H:%%M:%%S,%%l} or %%d{%%d %%m %%Y %%H:%%M:%%S,%%l}.
         *  If no date format specifier is given then the following format is used:
         *  "Wed Jan 02 02:03:55 1980". The date format specifier admits the same syntax
         *  as the ANSI C function strftime, with 1 addition. The addition is the specifier
         *  %%l for milliseconds, padded with zeros to make 3 digits.</li>
         * <li><b>%%m</b> - the message</li>
         * <li><b>%%n</b> - the platform specific line separator</li>
         * <li><b>%%p</b> - the priority</li>
         * <li><b>%%r</b> - milliseconds since this layout was created.</li>
         * <li><b>%%R</b> - seconds since Jan 1, 1970</li>
         * <li><b>%%u</b> - clock ticks since process start</li>
         * <li><b>%%x</b> - the NDC</li>
         * @param conversionPattern the conversion pattern
         * @exception ConfigureFailure if the pattern is invalid
         **/
        virtual void setConversionPattern(const std::string& conversionPattern)
            throw(ConfigureFailure);

        virtual std::string getConversionPattern() const;

        virtual void clearConversionPattern();

        class logcpp_API PatternComponent {
            public:
            inline virtual ~PatternComponent() {};
            virtual void append(std::ostringstream& out, const LoggingEvent& event) = 0;
        };

        private:
        typedef std::vector<PatternComponent*> ComponentVector;
        ComponentVector _components;

        std::string _conversionPattern;
    };
}

#endif // CF_LOGCPP_PATTERNLAYOUT_HH
