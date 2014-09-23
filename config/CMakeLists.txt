# file for general configuration for the whole framework
#
# config options that affect all the framework libraries and components
# must be handled in this directory

project(config)

ADD_SUBDIRECTORY ( sysinfo )


string ( TOUPPER "${CMAKE_BUILD_TYPE}" CF_BUILD_TYPE )

# explicit template support
IF ( CF_ENABLE_EXPLICIT_TEMPLATES AND CF_CXX_SUPPORTS_EXPLICIT_TEMPLATES )
  SET ( CF_HAVE_CXX_EXPLICIT_TEMPLATES ON CACHE BOOL "Support for Explicit templates activated" )
ELSE ()
  SET ( CF_HAVE_CXX_EXPLICIT_TEMPLATES OFF CACHE BOOL "Support for Explicit templates deactivated" )
ENDIF()

IF (CF_HAVE_SYS_TIME_H AND CF_HAVE_TIME_H)
  SET ( CF_TIME_WITH_SYS_TIME 1 CACHE BOOL "Have time.h and sys/time.h together")
ENDIF (CF_HAVE_SYS_TIME_H AND CF_HAVE_TIME_H)

#coolfluid_config.h
configure_file(${config_SOURCE_DIR}/coolfluid_config.h.in ${COOLFluiD_BINARY_DIR}/${CF_CONFIG_HH_FILENAME})

# create the header with the svn version
configure_file(${config_SOURCE_DIR}/coolfluid_svnversion.hh.in ${COOLFluiD_BINARY_DIR}/${CF_SVN_HH_FILENAME})
