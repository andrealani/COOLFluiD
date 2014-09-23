# find modules
FILE ( GLOB PLUGIN_MODULES "*/CMakeLists.txt" )
LIST ( SORT PLUGIN_MODULES )
FOREACH( ADIR ${PLUGIN_MODULES} )
  STRING (REGEX REPLACE "(/)*CMakeLists.txt" "" MODDIR ${ADIR})
  IF (IS_DIRECTORY ${MODDIR} )
    STRING (REGEX REPLACE "${CMAKE_CURRENT_SOURCE_DIR}/" "" FINALDIR ${MODDIR})
    LIST ( APPEND CF_MODULES_LIST ${FINALDIR} )
  ENDIF()
ENDFOREACH()

# CF_DEBUG_VAR(Subversion_FOUND)
# CF_DEBUG_VAR(CF_ENABLE_AUTOMATIC_UPDATE_MODULES)
# CF_DEBUG_VAR(COOLFluiD_WC_REVISION)

# in case cmake should use subversion to update to correct version
# typically from the CTest automatic testing procedure
IF ( Subversion_FOUND AND CF_ENABLE_AUTOMATIC_UPDATE_MODULES AND DEFINED COOLFluiD_WC_REVISION )

    EXECUTE_PROCESS(COMMAND ${Subversion_SVN_EXECUTABLE} update -r${COOLFluiD_WC_REVISION} ${CF_MODULES_LIST}
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
      OUTPUT_VARIABLE svn_update_output
      ERROR_VARIABLE  svn_update_error )

#     MESSAGE ( "SVN output: ${svn_update_output}" )
#     MESSAGE ( "SVN errors: ${svn_update_error}"  )

ENDIF()

# include all plugin configuration
FILE ( GLOB_RECURSE PLUGIN_CMAKE_CFGS "*/*.cmake" )
LIST ( SORT PLUGIN_CMAKE_CFGS )
FOREACH( ACFG ${PLUGIN_CMAKE_CFGS} )
  INCLUDE(${ACFG})
ENDFOREACH()

# include modules directories
FOREACH( amodule ${CF_MODULES_LIST} )
    LOGVERBOSE ( "\#  MODULE [${amodule}]")
    ADD_SUBDIRECTORY (${amodule})
ENDFOREACH()


