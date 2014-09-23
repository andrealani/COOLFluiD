##############################################################################
# this macro separates the sources form the headers
##############################################################################
MACRO( CF_SEPARATE_SOURCES FILELIST TGTNAME )

  FOREACH( AFILE ${FILELIST} )

    # check for existance of all declared files
    IF (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${AFILE})

      # check that file lenght is not too big
      CF_CHECK_FILE_LENGTH(${AFILE})

      # separate headers and sources
      IF(${AFILE} MATCHES "(\\.hh|\\.ci|\\.h|\\.hpp)")
        LIST ( APPEND ${TGTNAME}_headers ${AFILE})
      ENDIF()
      IF(${AFILE} MATCHES "(\\.cxx|\\.cpp|\\.cc|\\.c|\\.f|\\.f77|\\.f90)")
        LIST ( APPEND ${TGTNAME}_sources ${AFILE})
      ENDIF()

    ELSE (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${AFILE})
      MESSAGE ( FATAL_ERROR "In directory ${CMAKE_CURRENT_SOURCE_DIR} file ${AFILE} was declared in CMakeLists.txt but not found" )
    ENDIF (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${AFILE})

  ENDFOREACH(AFILE)

ENDMACRO( CF_SEPARATE_SOURCES )
##############################################################################

