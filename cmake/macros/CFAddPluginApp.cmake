##############################################################################
# macro for adding a application in the project
##############################################################################
MACRO( CF_ADD_PLUGIN_APP APPNAME )

  # option to build it or not
  OPTION ( CF_BUILD_${APPNAME} "Build the ${APPNAME} application" ON )

  # add to list of local apps
  LIST ( APPEND CF_LOCAL_APPNAMES ${APPNAME} )

#   CF_DEBUG_VAR(CF_MODULES_LIST)

  # check if all required modules are present
  SET ( ${APPNAME}_all_mods_pres ON )
  FOREACH ( reqmod ${${APPNAME}_requires_mods} )
    LIST ( FIND CF_MODULES_LIST ${reqmod} pos )
    IF ( ${pos} EQUAL -1 )
      SET ( ${APPNAME}_all_mods_pres OFF )
      IF ( CF_BUILD_${APPNAME} )
          LOGVERBOSE ( "     \# app [${APPNAME}] requires module [${reqmod}] which is not present")
      ENDIF()
    ENDIF()
  ENDFOREACH ( reqmod ${${APPNAME}_requires_mods} )

  SET ( ${APPNAME}_dir ${CMAKE_CURRENT_SOURCE_DIR} )

  IF ( CF_BUILD_${APPNAME} AND ${APPNAME}_all_mods_pres)
    SET ( ${APPNAME}_will_compile ON )
  ELSE()
    SET ( ${APPNAME}_will_compile OFF )
  ENDIF()

  LOGVERBOSE ("app_${APPNAME} = ${${APPNAME}_will_compile}")

  # compile if selected and all required modules are present
  IF (${APPNAME}_will_compile)

    IF( DEFINED ${APPNAME}_includedirs )
      INCLUDE_DIRECTORIES(${${APPNAME}_includedirs})
    ENDIF()

    CF_SEPARATE_SOURCES("${${APPNAME}_files}" ${APPNAME})

    SOURCE_GROUP ( Headers FILES ${${APPNAME}_headers} )
    SOURCE_GROUP ( Sources FILES ${${APPNAME}_sources} )

    LOG ( "    APP : [${APPNAME}]" )

    ADD_EXECUTABLE ( ${APPNAME} ${${APPNAME}_sources} ${${APPNAME}_headers} )

    # add installation paths
    INSTALL( TARGETS ${APPNAME}
      RUNTIME DESTINATION ${CF_INSTALL_BIN_DIR}
      LIBRARY DESTINATION ${CF_INSTALL_LIB_DIR}
      ARCHIVE DESTINATION ${CF_INSTALL_LIB_DIR}
      )

    # if mpi was found add it to the libraries
    IF   (CF_HAVE_MPI AND NOT CF_HAVE_MPI_COMPILER)
#           MESSAGE ( STATUS "${APPNAME} links to ${MPI_LIBRARIES}" )
          TARGET_LINK_LIBRARIES ( ${APPNAME} ${MPI_LIBRARIES} )
    ENDIF()

    # add external dependency libraries if defined
    IF( DEFINED ${APPNAME}_libs )
      TARGET_LINK_LIBRARIES ( ${APPNAME} ${${APPNAME}_libs} )
    ENDIF()

    # internal dependencies
    IF ( CF_ENABLE_INTERNAL_DEPS )
      IF( DEFINED ${APPNAME}_cflibs )
        TARGET_LINK_LIBRARIES ( ${APPNAME} ${${APPNAME}_cflibs} )
      ENDIF()
    ELSE()
      TARGET_LINK_LIBRARIES ( ${APPNAME} ${CF_KERNEL_LIBS} ${CF_KERNEL_STATIC_LIBS} ${CF_Boost_LIBRARIES} )
    ENDIF()

  ENDIF()

  GET_TARGET_PROPERTY ( ${APPNAME}_LINK_LIBRARIES  ${APPNAME} LINK_LIBRARIES )

  # log some info about the app
  LOGFILE("${APPNAME} : [${CF_BUILD_${APPNAME}}]")
  LOGFILE("${APPNAME} : [${${APPNAME}_will_compile}]")
  LOGFILE("${APPNAME}_dir             : [${${APPNAME}_dir}]")
  LOGFILE("${APPNAME}_includedirs     : [${${APPNAME}_includedirs}]")
  LOGFILE("${APPNAME}_libs            : [${${APPNAME}_libs}]")
  LOGFILE("${APPNAME}_all_mods_pres   : [${${APPNAME}_all_mods_pres}]")
  LOGFILE("${APPNAME}_requires_mods   : [${${APPNAME}_requires_mods}]")
  LOGFILE("${APPNAME}_sources         : [${${APPNAME}_sources}]")
  LOGFILE("${APPNAME}_LINK_LIBRARIES  : [${${APPNAME}_LINK_LIBRARIES}]")


ENDMACRO( CF_ADD_PLUGIN_APP )
##############################################################################
