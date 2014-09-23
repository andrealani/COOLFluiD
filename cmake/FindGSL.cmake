#Sets:
# GSL_INCLUDE_DIR  = where gsl_const.h can be found
# GSL_LIBRARY      = the library to link against (gsl etc)
# CF_HAVE_GSL        = set to true after finding the library

OPTION ( CF_SKIP_GSL "Skip search for gsl library" OFF )

IF ( NOT CF_SKIP_GSL )

  SET_TRIAL_INCLUDE_PATH ("") # clear include search path
  SET_TRIAL_LIBRARY_PATH ("") # clear library search path

  ADD_TRIAL_INCLUDE_PATH( "${GSL_INCLUDEDIR}" )
  ADD_TRIAL_INCLUDE_PATH( "${GSL_HOME}/include" )
  ADD_TRIAL_INCLUDE_PATH( "$ENV{GSL_HOME}/include" )
  ADD_TRIAL_INCLUDE_PATH( "${GSL_HOME}/include/gsl" )
  ADD_TRIAL_INCLUDE_PATH( "$ENV{GSL_HOME}/include/gsl" )

  FIND_PATH(GSL_INCLUDE_DIR gsl/gsl_const.h ${TRIAL_INCLUDE_PATHS}  NO_DEFAULT_PATH)
  FIND_PATH(GSL_INCLUDE_DIR gsl/gsl_const.h)

  ADD_TRIAL_LIBRARY_PATH( "${GSL_LIBRARYDIR}" )
  ADD_TRIAL_LIBRARY_PATH( ${GSL_HOME}/bin ${GSL_HOME}/lib $ENV{GSL_HOME}/bin ${GSL_HOME}/lib ${GSL_HOME}/lib64 )

  FIND_LIBRARY(GSL_LIBRARY gsl ${TRIAL_LIBRARY_PATHS} NO_DEFAULT_PATH)
  FIND_LIBRARY(GSL_LIBRARY gsl )

  IF(GSL_INCLUDE_DIR AND GSL_LIBRARY)
    SET(CF_HAVE_GSL 1 CACHE BOOL "Found gsl library")
  ELSE()
    SET(CF_HAVE_GSL 0 CACHE BOOL "Not found gsl library")
  ENDIF()
  
ELSE()
    SET(CF_HAVE_GSL 0 CACHE BOOL "Skipped gsl library")
ENDIF()

  MARK_AS_ADVANCED(
    GSL_INCLUDE_DIR
    GSL_LIBRARY
    CF_HAVE_GSL
  )

  LOG ( "CF_HAVE_GSL: [${CF_HAVE_GSL}]" )
  IF(CF_HAVE_GSL)
    LOG ( "  GSL_INCLUDE_DIR:  [${GSL_INCLUDE_DIR}]" )
    LOG ( "  GSL_LIBRARY:      [${GSL_LIBRARY}]" )
  ENDIF(CF_HAVE_GSL)
