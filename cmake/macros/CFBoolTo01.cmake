##############################################################################
# sets all the given variables either to 0 or 1
##############################################################################
MACRO(CF_BOOL_TO_01 FOUND_VAR )
   FOREACH (CURRENT_VAR ${ARGN})
      IF(${FOUND_VAR})
         SET(${CURRENT_VAR} 1)
      ELSE(${FOUND_VAR})
         SET(${CURRENT_VAR} 0)
      ENDIF(${FOUND_VAR})
   ENDFOREACH(CURRENT_VAR)
ENDMACRO(CF_BOOL_TO_01)
##############################################################################
