/*
 * SimpleLayout.cpp
 *
 * Copyright 2000, LifeLine Networks BV (www.lifeline.nl). All rights reserved.
 * Copyright 2000, Bastiaan Bakker. All rights reserved.
 *
 * See the COPYING file for the terms of usage and distribution.
 */

#include "PortabilityImpl.hh"
#include <logcpp/SimpleLayout.hh>
#include <logcpp/Priority.hh>
#include <sstream>

namespace logcpp {

    SimpleLayout::SimpleLayout() {
    }

    SimpleLayout::~SimpleLayout() {
    }

    std::string SimpleLayout::format(const LoggingEvent& event) {
        std::ostringstream message;

        const std::string& priorityName = Priority::getPriorityName(event.priority);
        message << priorityName << " - " << event.message << std::endl;

        return message.str();
    }
}
