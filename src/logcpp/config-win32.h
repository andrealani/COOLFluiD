#ifndef _INCLUDECF_LOGCPP_CONFIG_WIN32_H
#define _INCLUDECF_LOGCPP_CONFIG_WIN32_H 1

/* Define if you have the `ftime' function. */
#ifndef LOGCPP_HAVE_FTIME
#define LOGCPP_HAVE_FTIME  1
#endif

/* Define if you have the <io.h> header file.  */
#ifndef LOGCPP_HAVE_IO_H
#define LOGCPP_HAVE_IO_H 1
#endif

/* define to get around problems with ERROR in windows.h */
#ifndef LOGCPP_FIX_ERROR_COLLISION
#define LOGCPP_FIX_ERROR_COLLISION 1
#endif

/* define WIN32 for Borland */
#ifndef WIN32
#define WIN32
#endif

/* MSVCs <cstdlib> and <cstring> headers are broken in the sense that they
   put functions in the global namespace instead of std::
   The #defines below enable a workaround for MSVC 6 and lower. If MSVC 7
   is still broken please adjust the _MSC_VER version check and report it.
   See also bug report #628211.
*/
#if defined(_MSC_VER) && _MSC_VER < 1300

  #ifndef LOGCPP_CSTDLIB_NOT_IN_STD
  #define LOGCPP_CSTDLIB_NOT_IN_STD
  #endif

  #ifndef LOGCPP_CSTRING_NOT_IN_STD
  #define LOGCPP_CSTRING_NOT_IN_STD
  #endif

  #ifndef LOGCPP_CTIME_NOT_IN_STD
  #define LOGCPP_CTIME_NOT_IN_STD
  #endif

  #ifndef LOGCPP_CMATH_NOT_IN_STD
  #define LOGCPP_CMATH_NOT_IN_STD
  #endif

#endif

/* define mode_t. Move to Portability.hh if more platforms need it */
typedef int mode_t;

#endif /* _INCLUDECF_LOGCPP_CONFIG_WIN32_H */
