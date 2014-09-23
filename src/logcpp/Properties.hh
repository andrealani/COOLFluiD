// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

/*
 * Properties.hh
 *
 * Copyright 2002, Log4cpp Project. All rights reserved.
 *
 * See the COPYING file for the terms of usage and distribution.
 */

#ifndef CF_LOGCPP_PROPERTIES_HH
#define CF_LOGCPP_PROPERTIES_HH

#include "PortabilityImpl.hh"
#include <string>
#include <iostream>
#include <map>

namespace logcpp {
    
    class Properties : public std::map<std::string, std::string> {
        public:
        Properties();
        virtual ~Properties();

        virtual void load(std::istream& in);
        virtual void save(std::ostream& out);

        virtual int getInt(const std::string& property, int defaultValue);
        virtual bool getBool(const std::string& property, bool defaultValue);
        virtual std::string getString(const std::string& property,
                                      const char* defaultValue);

        protected:
        virtual void _substituteVariables(std::string& value);
    };
}

#endif // CF_LOGCPP_PROPERTIES_HH

