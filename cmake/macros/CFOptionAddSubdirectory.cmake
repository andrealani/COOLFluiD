##############################################################################
# only adds the subdirectory if it is present
##############################################################################
MACRO(OPTIONAL_ADD_SUBDIRECTORY OPT_ADIR )

  # IS_DIRECTORY works with full paths
  IF( IS_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/${OPT_ADIR} )
    ADD_SUBDIRECTORY (${OPT_ADIR})
  ENDIF( IS_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/${OPT_ADIR} )

ENDMACRO(OPTIONAL_ADD_SUBDIRECTORY)
##############################################################################
