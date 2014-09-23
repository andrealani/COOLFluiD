MACRO(LOGSYSINFO line)
   LOGVERBOSE (STATUS ${line})
   FILE(APPEND ${CMAKE_CURRENT_BINARY_DIR}/SystemInfo.txt "${line}\n")
ENDMACRO(LOGSYSINFO line)

# recreate the file
FILE( WRITE ${CMAKE_CURRENT_BINARY_DIR}/SystemInfo.txt "COOLFluiD system information file\n")

LOGSYSINFO("---------------------------------------------------------")
LOGSYSINFO("System information: ")
LOGSYSINFO("---------------------------------------------------------")

LOGSYSINFO("C compiler: ${CMAKE_C_COMPILER}")
LOGSYSINFO("CXX compiler: ${CMAKE_CXX_COMPILER}")
LOGSYSINFO("Fortran compiler: ${CMAKE_Fortran_COMPILER}")
LOGSYSINFO("C compiler is gcc: ${CMAKE_COMPILER_IS_GNUCC}")
LOGSYSINFO("CXX compiler is g++: ${CMAKE_COMPILER_IS_GNUCXX}")
LOGSYSINFO("System: ${CMAKE_SYSTEM}")
LOGSYSINFO("System name: ${CMAKE_SYSTEM_NAME}")
LOGSYSINFO("System version: ${CMAKE_SYSTEM_VERSION}")
LOGSYSINFO("Shared Library Link: ${CMAKE_SHARED_LIBRARY_LINK}")

IF(UNIX)
   LOGSYSINFO("This is some kind of UNIX environment")
ENDIF(UNIX)

IF(WIN32)
   LOGSYSINFO("This is Windows")
ENDIF(WIN32)

IF(APPLE)
   LOGSYSINFO("This is an Apple")
ENDIF(APPLE)

IF(MINGW)
   LOGSYSINFO("This is MinGW")
ENDIF(MINGW)

IF(CYGWIN)
   LOGSYSINFO("This is Cygwin")
ENDIF(CYGWIN)

IF(BORLAND)
   LOGSYSINFO("This is for a Borland compiler")
ENDIF(BORLAND)

LOGSYSINFO("---------------------------------------------------------")

