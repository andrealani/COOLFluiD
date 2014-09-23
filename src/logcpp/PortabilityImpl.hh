// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

/*
 * PortabilityImpl.hh
 *
 * Copyright 2002, Log4cpp Project. All rights reserved.
 *
 * See the COPYING file for the terms of usage and distribution.
 */

#ifndef CF_LOGCPP_PORTABILITYIMPL_HH
#define CF_LOGCPP_PORTABILITYIMPL_HH

#include <logcpp/Portability.hh>

#ifdef LOGCPP_CSTDLIB_NOT_IN_STD
#include <cstdlib>
namespace std {
    static inline char *getenv(const char *name) { return ::getenv(name); };
    static inline int atoi(const char *nptr) { return ::atoi(nptr); };
    static inline unsigned long int
        strtoul(const char *nptr, char **endptr, int base) { 
        return ::strtol(nptr, endptr, base); 
    };
}
#endif
    
#ifdef LOGCPP_CSTRING_NOT_IN_STD
#include <cstring>
namespace std {
    static inline void *memmove(void *dest, const void *src, size_t n) {
        return ::memmove(dest, src, n);
    };
}
#endif

#ifdef LOGCPP_CTIME_NOT_IN_STD
#include <ctime>
namespace std {
    static inline size_t strftime(char *strDest, size_t maxsize, const char *format, const struct tm *timeptr ) {
        return ::strftime(strDest,maxsize,format,timeptr);
    }
    static inline struct tm *localtime( const time_t *timer ) { return ::localtime(timer); }
}
#endif

#ifdef LOGCPP_CMATH_NOT_IN_STD
#include <cmath>
namespace std {
    static inline int abs(int i) { return ::abs(i); }
}
#endif

#endif // CF_LOGCPP_PORTABILITYIMPL_HH
