#
# this module look for METIS (http://hdf.ncsa.uiuc.edu) support
# it will define the following values
#
# METIS_INCLUDE_DIR  = where metis.h can be found
# METIS_LIBRARY      = the library to link against (metis etc)
# CF_HAVE_METIS        = set to true after finding the library
#

OPTION ( CF_SKIP_METIS "Skip search for Metis library" OFF )
IF ( NOT CF_SKIP_METIS )

  SET_TRIAL_INCLUDE_PATH ("") # clear include search path
  SET_TRIAL_LIBRARY_PATH ("") # clear library search path

  ADD_TRIAL_INCLUDE_PATH( ${METIS_HOME}/include )
  ADD_TRIAL_INCLUDE_PATH( $ENV{METIS_HOME}/include )

  FIND_PATH(METIS_INCLUDE_DIR metis.h ${TRIAL_INCLUDE_PATHS}  NO_DEFAULT_PATH)
  FIND_PATH(METIS_INCLUDE_DIR metis.h)

  ADD_TRIAL_LIBRARY_PATH(${METIS_HOME}/lib )
  ADD_TRIAL_LIBRARY_PATH($ENV{METIS_HOME}/lib )

  FIND_LIBRARY(METIS_LIBRARY metis ${TRIAL_LIBRARY_PATHS} NO_DEFAULT_PATH)
  FIND_LIBRARY(METIS_LIBRARY metis )

  IF(METIS_INCLUDE_DIR AND METIS_LIBRARY)
    SET(CF_HAVE_METIS 1 CACHE BOOL "Found metis library")
  ELSE()
    SET(CF_HAVE_METIS 0 CACHE BOOL "Not fount metis library")
  ENDIF()

ELSE()
    SET(CF_HAVE_METIS 0 CACHE BOOL "Skipped METIS library")
ENDIF()

MARK_AS_ADVANCED(
  METIS_INCLUDE_DIR
  METIS_LIBRARY
  CF_HAVE_METIS
)

LOG ( "CF_HAVE_METIS: [${CF_HAVE_METIS}]" )
IF(CF_HAVE_METIS)
  LOGFILE ( "  METIS_INCLUDE_DIR: [${METIS_INCLUDE_DIR}]" )
  LOGFILE ( "  METIS_LIBRARY: [${METIS_LIBRARY}]" )
ENDIF()
