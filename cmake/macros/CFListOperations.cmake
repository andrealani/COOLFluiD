##############################################################################
# this macro adds to a cached list if element not yet present
##############################################################################
MACRO( CF_CACHE_LIST_REMOVE THELIST THEVAR )
  IF ( DEFINED ${THELIST} )
    LIST( REMOVE_ITEM ${THELIST} ${THEVAR} )
  ENDIF ( DEFINED ${THELIST} )
ENDMACRO( CF_CACHE_LIST_REMOVE THELIST THEVAR )
##############################################################################

##############################################################################
# this macro adds to a cached list if element not yet present
##############################################################################
MACRO( CF_CACHE_LIST_APPEND THELIST THEVAR )
  CF_CACHE_LIST_REMOVE(${THELIST} ${THEVAR})
  SET ( ${THELIST} ${${THELIST}} ${THEVAR} CACHE INTERNAL "" FORCE )
ENDMACRO( CF_CACHE_LIST_APPEND THELIST THEVAR )
##############################################################################