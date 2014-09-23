##############################################################################
# package instrucitons
# if CPack is available and this is a build of just COOLFluiD
# as opposed to a build of from an external project
##############################################################################

IF(EXISTS "${CMAKE_ROOT}/Modules/CPack.cmake")
  IF("${COOLFluiD_BINARY_DIR}" STREQUAL "${CMAKE_BINARY_DIR}")

      SET (CPACK_GENERATOR TGZ)
      SET (CPACK_PACKAGE_DESCRIPTION_SUMMARY "COOLFluiD - An environment for computational fluid dynamics")
      SET (CPACK_PACKAGE_VENDOR "von Karman Institute for Fluid Dynamics")
#       SET(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_CURRENT_SOURCE_DIR}/Copyright.txt")
#       SET(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/Copyright.txt")
      SET (CPACK_PACKAGE_VERSION_MAJOR ${COOLFLUID_VERSION_MAJOR})
      SET (CPACK_PACKAGE_VERSION_MINOR ${COOLFLUID_VERSION_MINOR})
      SET (CPACK_PACKAGE_VERSION_PATCH ${COOLFLUID_VERSION_RELEASE})
      SET(CPACK_PACKAGE_INSTALL_DIRECTORY "coolfluid-${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}")
      SET(CPACK_SOURCE_PACKAGE_FILE_NAME  "coolfluid-${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}")
#       SET(CPACK_PACKAGE_EXECUTABLES "vtk" "VTK" )

      INCLUDE(CPack)

  ENDIF("${COOLFluiD_BINARY_DIR}" STREQUAL "${CMAKE_BINARY_DIR}")
ENDIF(EXISTS "${CMAKE_ROOT}/Modules/CPack.cmake")
