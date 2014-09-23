// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

/*
 * Portability.hh
 *
 * Copyright 2001, LifeLine Networks BV (www.lifeline.nl). All rights reserved.
 * Copyright 2001, Bastiaan Bakker. All rights reserved.
 *
 * See the COPYING file for the terms of usage and distribution.
 */

#ifndef CF_LOGCPP_PORTABILITY_HH
#define CF_LOGCPP_PORTABILITY_HH

#if defined (_MSC_VER) || defined(__BORLANDC__)
#    include <logcpp/config-win32.h>
#else
#    include <logcpp/config.h>
#endif

#include <logcpp/Export.hh>

#if defined(_MSC_VER)
#    pragma warning( disable : 4786 ) // 255 char debug symbol limit
#    pragma warning( disable : 4290 ) // throw specifier not implemented
#    pragma warning( disable : 4251 ) // "class XXX should be exported"
#endif

#endif
