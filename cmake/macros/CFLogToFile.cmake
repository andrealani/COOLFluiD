##############################################################################
# this macro logs simultaneously to screen and to file
##############################################################################
MACRO ( LOG line )
   MESSAGE(STATUS ${line})
   FILE(APPEND ${PROJECT_LOG_FILE} "${line}\n")
ENDMACRO ( LOG )

##############################################################################
# this macro logs just to file
##############################################################################
MACRO ( LOGFILE line )
   FILE(APPEND ${PROJECT_LOG_FILE} "${line}\n")
ENDMACRO ( LOGFILE )
##############################################################################

##############################################################################
# this macro logs to screen if we are verbose
##############################################################################
MACRO ( LOGVERBOSE line )
    IF ( CF_CMAKE_VERBOSE )
      LOG ( ${line} )
    ELSE ()
      LOGFILE ( ${line} )
    ENDIF()
ENDMACRO ( LOGVERBOSE )
##############################################################################


