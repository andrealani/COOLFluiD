# sunversion version check
FIND_PACKAGE(Subversion)

IF (Subversion_FOUND)

  Subversion_WC_INFO(${COOLFluiD_SOURCE_DIR} COOLFluiD)
#   MESSAGE("Current revision is ${COOLFluiD_WC_REVISION}")
#   MESSAGE("svn info : ${COOLFluiD_WC_INFO}")

  FIND_PROGRAM(Subversion_SVNVERSION_EXECUTABLE svnversion DOC "subversion svnversion command line client")
  MARK_AS_ADVANCED(Subversion_SVNVERSION_EXECUTABLE)

  IF(Subversion_SVNVERSION_EXECUTABLE)
  EXECUTE_PROCESS(COMMAND ${Subversion_SVNVERSION_EXECUTABLE} -n ${COOLFluiD_SOURCE_DIR}
      WORKING_DIRECTORY ${COOLFluiD_SOURCE_DIR}
      OUTPUT_VARIABLE coolfluid_svnversion
      OUTPUT_STRIP_TRAILING_WHITESPACE)
  ELSE()
    LOG("Subversion svn command was found, but not svnversion.")
    SET(coolfluid_svnversion "NOVERSION-FOUND")
  ENDIF()

ELSE (Subversion_FOUND)
  SET(coolfluid_svnversion "NOVERSION-FOUND")
ENDIF (Subversion_FOUND)
