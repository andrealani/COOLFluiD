##############################################################################
# this macro separates the sources form the headers
##############################################################################
MACRO( CF_CHECK_FILE_LENGTH SOURCE_FILE )

  STRING(LENGTH ${SOURCE_FILE} FILENAME_LENGTH)
  IF(FILENAME_LENGTH GREATER 90)
    MESSAGE( FATAL_ERROR "The filename [${SOURCE_FILE}] is too long (${FILENAME_LENGTH} chars).\nFilenames longer than 90 characters can cause portability issues with some versions of tar.")
  ENDIF(FILENAME_LENGTH GREATER 90)

ENDMACRO( CF_CHECK_FILE_LENGTH )
##############################################################################

