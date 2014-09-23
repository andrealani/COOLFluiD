##############################################################################
# adds a path to search for when searching for include files
##############################################################################
MACRO( ADD_TRIAL_INCLUDE_PATH ICPATH )
 IF( EXISTS ${ICPATH})
  LIST ( APPEND TRIAL_INCLUDE_PATHS ${ICPATH} )
 ENDIF( EXISTS ${ICPATH})
ENDMACRO( ADD_TRIAL_INCLUDE_PATH )
##############################################################################

##############################################################################
# adds a path to search for when searching for library files
##############################################################################
MACRO( ADD_TRIAL_LIBRARY_PATH ICPATH )
 IF( EXISTS ${ICPATH})
   LIST ( APPEND TRIAL_LIBRARY_PATHS ${ICPATH} )
 ENDIF( EXISTS ${ICPATH})
ENDMACRO( ADD_TRIAL_LIBRARY_PATH )
##############################################################################

##############################################################################
# sets a path to search for when searching for INCLUDE files
##############################################################################
MACRO( SET_TRIAL_INCLUDE_PATH IPATHS )
 SET ( TRIAL_INCLUDE_PATHS "" )
 FOREACH ( path ${IPATHS} )
   IF( EXISTS ${path})
     LIST ( APPEND TRIAL_INCLUDE_PATHS ${path} )
   ENDIF( EXISTS ${path})
 ENDFOREACH ( path ${IPATHS} )
ENDMACRO( SET_TRIAL_INCLUDE_PATH )
##############################################################################

##############################################################################
# sets a path to search for when searching for library files
##############################################################################
MACRO( SET_TRIAL_LIBRARY_PATH LPATHS )
 SET ( TRIAL_LIBRARY_PATHS "" )
 FOREACH ( path ${LPATHS} )
   IF( EXISTS ${path})
     LIST ( APPEND TRIAL_LIBRARY_PATHS ${path} )
   ENDIF( EXISTS ${path})
 ENDFOREACH ( path ${LPATHS} )
ENDMACRO( SET_TRIAL_LIBRARY_PATH )
##############################################################################