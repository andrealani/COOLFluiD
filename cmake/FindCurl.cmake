#Sets:
# CURL_INCLUDE_DIR  = where curl.h can be found
# CURL_LIBRARY      = the library to link against (curl etc)
# CF_HAVE_CURL        = set to true after finding the library

OPTION ( CF_SKIP_CURL "Skip search for Curl library" OFF )

IF ( NOT CF_SKIP_CURL )

  SET_TRIAL_INCLUDE_PATH ("") # clear include search path
  SET_TRIAL_LIBRARY_PATH ("") # clear library search path

  ADD_TRIAL_INCLUDE_PATH( ${CURL_HOME}/include )
  ADD_TRIAL_INCLUDE_PATH( $ENV{CURL_HOME}/include )
  ADD_TRIAL_INCLUDE_PATH( /usr/include )

  FIND_PATH(CURL_INCLUDE_DIR curl/curl.h ${TRIAL_INCLUDE_PATHS} NO_DEFAULT_PATH)

  IF ( NOT CURL_INCLUDE_DIR)
    FIND_PATH(CURL_INCLUDE_DIR curl/curl.h)
  ENDIF()

  ADD_TRIAL_LIBRARY_PATH( ${CURL_HOME}/bin )
  ADD_TRIAL_LIBRARY_PATH( ${CURL_HOME}/lib )
  ADD_TRIAL_LIBRARY_PATH( $ENV{CURL_HOME}/bin )
  ADD_TRIAL_LIBRARY_PATH( $ENV{CURL_HOME}/lib )
  ADD_TRIAL_LIBRARY_PATH( /usr/lib64 )

  FIND_LIBRARY(CURL_LIBRARY curl ${TRIAL_LIBRARY_PATHS} NO_DEFAULT_PATH)

  IF ( NOT CURL_LIBRARY)
    FIND_PATH(CURL_LIBRARY curl)
  ENDIF()

  IF(CURL_INCLUDE_DIR AND CURL_LIBRARY)
    SET(CF_HAVE_CURL 1 CACHE BOOL "Found curl library")
  ELSE()
    SET(CF_HAVE_CURL 0 CACHE BOOL "Not found curl library")
  ENDIF()
  
ELSE()
    SET(CF_HAVE_CURL 0 CACHE BOOL "Skipped Curl library")
ENDIF()

  MARK_AS_ADVANCED(
    CURL_INCLUDE_DIR
    CURL_LIBRARY
    CF_HAVE_CURL
  )

  LOG ( "CF_HAVE_CURL: [${CF_HAVE_CURL}]" )
  IF(CF_HAVE_CURL)
    LOG ( "  CURL_INCLUDE_DIR:  [${CURL_INCLUDE_DIR}]" )
    LOG ( "  CURL_LIBRARY:      [${CURL_LIBRARY}]" )
  ENDIF(CF_HAVE_CURL)
