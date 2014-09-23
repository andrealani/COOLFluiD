LOG ( "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++" )
LOG ( "Checking compiler features:" )
LOG ( "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++" )

#######################################################################################

  LOG ( "+++++  Checking sizes of types" )
  CHECK_TYPE_SIZE("void *"       CF_SIZEOF_PTR)
  CHECK_TYPE_SIZE(char           CF_SIZEOF_CHAR)
  CHECK_TYPE_SIZE(short          CF_SIZEOF_SHORT)
  CHECK_TYPE_SIZE(int            CF_SIZEOF_INT)
  CHECK_TYPE_SIZE(long           CF_SIZEOF_LONG)
  CHECK_TYPE_SIZE(float          CF_SIZEOF_FLOAT)
  CHECK_TYPE_SIZE(double         CF_SIZEOF_DOUBLE)
  CHECK_TYPE_SIZE("long double"  CF_SIZEOF_LONG_DOUBLE)
#   CHECK_TYPE_SIZE(int8_t         CF_SIZEOF_INT8_T)
#   CHECK_TYPE_SIZE(uint8_t        CF_SIZEOF_UINT8_T)
#   CHECK_TYPE_SIZE(int_least8_t   CF_SIZEOF_INT_LEAST8_T)
#   CHECK_TYPE_SIZE(uint_least8_t  CF_SIZEOF_UINT_LEAST8_T)
#   CHECK_TYPE_SIZE(int_fast8_t    CF_SIZEOF_INT_FAST8_T)
#   CHECK_TYPE_SIZE(uint_fast8_t   CF_SIZEOF_UINT_FAST8_T)
#   CHECK_TYPE_SIZE(int16_t        CF_SIZEOF_INT16_T)
#   CHECK_TYPE_SIZE(uint16_t       CF_SIZEOF_UINT16_T)
#   CHECK_TYPE_SIZE(int_least16_t  CF_SIZEOF_INT_LEAST16_T)
#   CHECK_TYPE_SIZE(uint_least16_t CF_SIZEOF_UINT_LEAST16_T)
#   CHECK_TYPE_SIZE(int_fast16_t   CF_SIZEOF_INT_FAST16_T)
#   CHECK_TYPE_SIZE(uint_fast16_t  CF_SIZEOF_UINT_FAST16_T)
#   CHECK_TYPE_SIZE(int32_t        CF_SIZEOF_INT32_T)
#   CHECK_TYPE_SIZE(uint32_t       CF_SIZEOF_UINT32_T)
#   CHECK_TYPE_SIZE(int_least32_t  CF_SIZEOF_INT_LEAST32_T)
#   CHECK_TYPE_SIZE(uint_least32_t CF_SIZEOF_UINT_LEAST32_T)
#   CHECK_TYPE_SIZE(int_fast32_t   CF_SIZEOF_INT_FAST32_T)
#   CHECK_TYPE_SIZE(uint_fast32_t  CF_SIZEOF_UINT_FAST32_T)
#   CHECK_TYPE_SIZE(int64_t        CF_SIZEOF_INT64_T)
#   CHECK_TYPE_SIZE(uint64_t       CF_SIZEOF_UINT64_T)
#   CHECK_TYPE_SIZE(int_least64_t  CF_SIZEOF_INT_LEAST64_T)
#   CHECK_TYPE_SIZE(uint_least64_t CF_SIZEOF_UINT_LEAST64_T)
#   CHECK_TYPE_SIZE(int_fast64_t   CF_SIZEOF_INT_FAST64_T)
#   CHECK_TYPE_SIZE(uint_fast64_t  CF_SIZEOF_UINT_FAST64_T)
  CHECK_TYPE_SIZE(size_t         CF_SIZEOF_SIZE_T)
#   CHECK_TYPE_SIZE(ssize_t        CF_SIZEOF_SSIZE_T)
#   CHECK_TYPE_SIZE(off_t          CF_SIZEOF_OFF_T)
#   CHECK_TYPE_SIZE(__int64        CF_SIZEOF___INT64)
  CHECK_TYPE_SIZE("long long"    CF_SIZEOF_LONG_LONG)

  MATH ( EXPR CF_OS_BITS "${CF_SIZEOF_PTR} * 8")

#######################################################################################

  LOG ( "+++++  Checking MPI support" )
  INCLUDE(CheckMPI)

#######################################################################################

  # LOG ( "+++++  Checking for pre compiled header support" )
  # INCLUDE(CheckPreCompiledHeaderSupport)

#######################################################################################

  LOG ( "+++++  Checking C++ compiler has namespaces" )
  CHECK_CXX_SOURCE_COMPILES (
  " namespace lolo { struct popo { int i; };  }
    using namespace lolo;
    int main(int argc, char* argv[])
    {
      lolo::popo p;
      popo pp;
      return 0;
    }
  "
  CF_CXX_HAVE_NAMESPACES)
  IF( NOT CF_CXX_HAVE_NAMESPACES )
    MESSAGE ( FATAL_ERROR "C++ compiler does not support namespaces" )
  ENDIF()

#######################################################################################

  LOG ( "+++++  Checking for __FUNCTION__ support" )
  INCLUDE(CheckFunctionDef)

#######################################################################################

  LOG ( "+++++  Checking for restrict keyword" )
  INCLUDE(CheckRestrictKeyword)

#######################################################################################

  # LOG ( "+++++  Checking for explicit template support" ) # output in macro
  INCLUDE ( CheckExplicitInstantiation )

#######################################################################################

  LOG ( "+++++  Checking for mmap support" ) # check memory mmap functions
  # check memory mmap functions
  CHECK_FUNCTION_EXISTS(mmap   CF_HAVE_MMAP)
  CHECK_FUNCTION_EXISTS(munmap CF_HAVE_MUNMAP)
  CHECK_FUNCTION_EXISTS(mremap CF_HAVE_MREMAP)
  IF(CF_HAVE_MMAP AND CF_HAVE_MUNMAP AND CF_HAVE_MREMAP)
    SET( CF_HAVE_ALLOC_MMAP 1 CACHE BOOL "MemoryAllocator_MMAP can be built")
  ENDIF()

#######################################################################################

  LOG ( "+++++  Checking for vsnprintf function" ) # check memory mmap functions
  CHECK_FUNCTION_EXISTS(vsnprintf   CF_HAVE_VSNPRINTF)

#######################################################################################

  LOG ( "+++++  Checking for erfc function" )
  CHECK_CXX_SOURCE_COMPILES (
  " #include <cmath>
    int main (int argc, char* argv[]) { erfc (0.); }
  "
  CF_HAVE_MATH_ERFC )

#######################################################################################

  LOG ( "+++++  Checking for asinh function" )
  CHECK_CXX_SOURCE_COMPILES (
  " #include <cmath>
    int main (int argc, char* argv[]) { asinh (0.); }
  "
  CF_HAVE_MATH_ASINH )

#######################################################################################

  LOG ( "+++++  Checking for acosh function" )
  CHECK_CXX_SOURCE_COMPILES (
  " #include <cmath>
    int main (int argc, char* argv[]) { acosh (0.); }
  "
  CF_HAVE_MATH_ACOSH )

#######################################################################################

  LOG ( "+++++  Checking for atanh function" )
  CHECK_CXX_SOURCE_COMPILES (
  " #include <cmath>
    int main (int argc, char* argv[]) { atanh (0.); }
  "
  CF_HAVE_MATH_ATANH)

#######################################################################################

  LOG ( "+++++  Checking for the POSIX unistd.h header" )    # check unistd.h
  CHECK_INCLUDE_FILE( unistd.h  CF_HAVE_UNISTD_H )

#######################################################################################

  # check for c++ abi, ussually present in GNU compilers
  LOG ( "+++++  Checking for cxxabi" )

  CHECK_CXX_SOURCE_COMPILES (
    "#include <cxxabi.h>
    int main(int argc, char* argv[])
    { char * type; int status;
      char * r = abi::__cxa_demangle(type, 0, 0, &status);
      return 0;
    }"
    CF_HAVE_CXXABI_H)

#######################################################################################

  # check for time headers
  LOG ( "+++++  Checking for headers with time information" )
  CHECK_INCLUDE_FILE(sys/time.h       CF_HAVE_SYS_TIME_H)
  CHECK_INCLUDE_FILE(time.h           CF_HAVE_TIME_H)
  CHECK_INCLUDE_FILE(sys/resource.h   CF_HAVE_SYS_RESOURCE_H)
  CHECK_FUNCTION_EXISTS(gettimeofday  CF_HAVE_GETTIMEOFDAY)

#######################################################################################
# Win32 specific
#######################################################################################

IF(WIN32)

  LOG ( "+++++  Checking for the Win32 windows.h header" )    # check windows.hfor Windows API
  CHECK_INCLUDE_FILE_CXX(windows.h CF_HAVE_WINDOWSH)
  LOG ( "+++++  Checking for the Win32 dbghelp.h header" )    # check dbghelp.h for call stack
  CHECK_INCLUDE_FILE_CXX(dbghelp.h CF_HAVE_DBGHELPH)
  LOG ( "+++++  Checking for the Win32 psapi.h header" )      # check psapi.h for memory info
  CHECK_INCLUDE_FILE_CXX(psapi.h CF_HAVE_PSAPIH)

ENDIF()

#######################################################################################
# UNIX specific
#######################################################################################

IF(UNIX)

  LOG ( "+++++  Checking for the dlfcn.h header" )  # for dlopen
  CHECK_INCLUDE_FILE_CXX(dlfcn.h CF_HAVE_DLOPEN)

ENDIF()

