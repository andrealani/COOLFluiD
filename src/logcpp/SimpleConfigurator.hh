// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

/*
 * SimpleConfigurator.hh
 *
 * Copyright 2001, Glen Scott. All rights reserved.
 *
 * See the COPYING file for the terms of usage and distribution.
 */
#ifndef CF_LOGCPP_SIMPLECONFIGURATOR_HH
#define CF_LOGCPP_SIMPLECONFIGURATOR_HH

#include <logcpp/Portability.hh>
#include <iostream>
#include <string>
#include <logcpp/Configurator.hh>

namespace logcpp {

    /**
     * This class implements a simple Configurator for logcpp.
     * It is a temporary hack with an undocumented configuration format.
     * @deprecated As of version 0.3.2 logcpp includes a log4j format 
     * compatible PropertyConfigurator, removing the need for 
     * SimpleConfigurator. This class will be removed in 0.4.0.
     **/
    class logcpp_API SimpleConfigurator {
    public:

        /**
         * Configure logcpp with the configuration in the given file.
         * NB. The configuration file format is undocumented and may change
         * without notice.
         * @since 0.2.6
         * @param initFileName name of the configuration file
         * @exception ConfigureFailure if the method encountered a read or 
         * syntax error.
         **/
        static void configure(const std::string& initFileName) throw (ConfigureFailure);

        /**
         * Configure logcpp with the configuration in the given file.
         * NB. The configuration file format is undocumented and may change
         * without notice.
         * @since 0.3.1
         * @param initFile an input stream to the configuration file
         * @exception ConfigureFailure if the method encountered a read or 
         * syntax error.
         **/
        static void configure(std::istream& initFile) throw (ConfigureFailure);    };
}

#endif
