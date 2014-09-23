##############################################################################
# macro for warning about missing files in the listing of files
# usefull for headers  that are forgotten and give errors compiling
# from installed sources.
# it warns about files that are present in the directory but are neither
# in the headers or the sources of all the libraries and applications
# defined in this directory
##############################################################################
MACRO( CF_WARN_ORPHAN_FILES )

  # first find all the files in the directory
  FOREACH( CFEXT ${CF_FILE_EXTENSIONS} )
    FILE ( GLOB THISEXT_FILES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.${CFEXT})
    LIST ( LENGTH THISEXT_FILES NELEMS )
#     CF_DEBUG_VAR(NELEMS)
    IF ( NELEMS GREATER 0 )
      LIST ( APPEND alldirfiles "${THISEXT_FILES}")
    ENDIF ( NELEMS GREATER 0 )
  ENDFOREACH(CFEXT)

  # remove all found files from orphan list to avoid duplicates
  FOREACH( AFILE ${alldirfiles} )
    CF_CACHE_LIST_REMOVE ( CF_ORPHAN_FILES ${CMAKE_CURRENT_SOURCE_DIR}/${AFILE} )
  ENDFOREACH(AFILE)

  # remove files present in optional list
  FOREACH( AFILE ${OPTIONAL_dirfiles} )
    LIST ( REMOVE_ITEM alldirfiles ${AFILE})
  ENDFOREACH(AFILE)

  # remove files present in libs
  FOREACH( LOCALLIB ${CF_LOCAL_LIBNAMES} )
    FOREACH( AFILE ${${LOCALLIB}_files} )
      LIST ( REMOVE_ITEM alldirfiles ${AFILE})
    ENDFOREACH(AFILE)
  ENDFOREACH( LOCALLIB ${CF_LOCAL_LIBNAMES} )

  # remove files present in apps
  FOREACH( LOCALAPP ${CF_LOCAL_APPNAMES} )
    FOREACH( AFILE ${${LOCALAPP}_files} )
      LIST ( REMOVE_ITEM alldirfiles ${AFILE})
    ENDFOREACH(AFILE)
  ENDFOREACH( LOCALAPP ${CF_LOCAL_APPNAMES} )

  # remove files present in tests
  FOREACH( LOCALTEST ${CF_LOCAL_TESTNAMES} )
    FOREACH( AFILE ${${LOCALTEST}_files} )
      LIST ( REMOVE_ITEM alldirfiles ${AFILE})
    ENDFOREACH(AFILE)
  ENDFOREACH( LOCALTEST ${CF_LOCAL_TESTNAMES} )

  # warn about the other files
  FOREACH( AFILE ${alldirfiles} )
    # temporarily ignore Test files while the unit tests are not added
    IF(${AFILE} MATCHES "^Test" )
    ELSE(${AFILE} MATCHES "^Test" )
      LOGFILE(" +++ WARNING : orphan file ${AFILE}")
      CF_CACHE_LIST_APPEND ( CF_ORPHAN_FILES ${CMAKE_CURRENT_SOURCE_DIR}/${AFILE} )
    ENDIF(${AFILE} MATCHES "^Test" )
  ENDFOREACH(AFILE)

ENDMACRO( CF_WARN_ORPHAN_FILES )
##############################################################################
