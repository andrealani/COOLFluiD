// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

/*
 * Configurator.hh
 *
 * Copyright 2001, Glen Scott. All rights reserved.
 *
 * See the COPYING file for the terms of usage and distribution.
 */
#ifndef CF_LOGCPP_CONFIGURATOR_HH
#define CF_LOGCPP_CONFIGURATOR_HH

#include <logcpp/Portability.hh>
#include <logcpp/Export.hh>
#include <string>
#include <stdexcept>

namespace logcpp {

    /**
     * Exception class for configuration.
     */
    class logcpp_API ConfigureFailure : public std::runtime_error {
        public:
        /**
         * Constructor.
         * @param reason String containing the description of the exception.
         */
        ConfigureFailure(const std::string& reason);
    };
}

#endif // CF_LOGCPP_CONFIGURATOR_HH
