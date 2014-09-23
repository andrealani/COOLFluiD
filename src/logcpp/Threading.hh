// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

/*
 * Threading.hh
 *
 * Copyright 2002, LifeLine Networks BV (www.lifeline.nl). All rights reserved.
 * Copyright 2002, Bastiaan Bakker. All rights reserved.
 *
 * See the COPYING file for the terms of usage and distribution.
 */

#ifndef CF_LOGCPP_THREADING_THREADING_HH
#define CF_LOGCPP_THREADING_THREADING_HH

#include <logcpp/Portability.hh>

#define CF_LOGCPP_HAVE_THREADING
#define CF_LOGCPP_USE_BOOSTTHREADS

#ifdef CF_LOGCPP_HAVE_THREADING

  #ifdef CF_LOGCPP_USE_BOOSTTHREADS
  #include <logcpp/BoostThreads.hh>
  #endif

#else /* CF_LOGCPP_HAVE_THREADING */
  #include <logcpp/DummyThreads.hh>
#endif /* CF_LOGCPP_HAVE_THREADING */

#endif // CF_LOGCPP_THREADING_THREADING_HH
