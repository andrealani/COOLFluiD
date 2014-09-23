// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

/*
 * BoostThreads.hh
 *
 * Copyright 2002, LifeLine Networks BV (www.lifeline.nl). All rights reserved.
 * Copyright 2002, Bastiaan Bakker. All rights reserved.
 *
 * See the COPYING file for the terms of usage and distribution.
 */

#ifndef CF_LOGCPP_THREADING_BOOSTTHREADS_HH
#define CF_LOGCPP_THREADING_BOOSTTHREADS_HH

#include <logcpp/Portability.hh>
#include <boost/thread.hpp>
#include <stdio.h>
#include <string>

namespace logcpp {
    namespace threading {

		std::string getThreadId();

        typedef boost::mutex Mutex;
        typedef boost::mutex::scoped_lock ScopedLock;

//        typedef template <typename T> class boost::thread_specific_ptr<T> ThreadLocalDataHolder;
    }
}
#endif
