##############################################################################
# prints the variable name and value
##############################################################################
MACRO( CF_DEBUG_VAR THE_VARIABLE )
  LOG ( " +++++ DEBUG +++++ ${THE_VARIABLE} : [${${THE_VARIABLE}}]")
ENDMACRO( CF_DEBUG_VAR THE_VARIABLE )

##############################################################################
# sets a variable if not yet defined
##############################################################################
MACRO    ( SET_VAR_IF_NOT_DEFINED variable value )
  IF ( NOT DEFINED ${variable} )
    SET ( ${variable} ${value} )
  ENDIF ( NOT DEFINED ${variable} )
ENDMACRO ( SET_VAR_IF_NOT_DEFINED )
