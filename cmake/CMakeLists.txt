# install the cmake files
FILE(GLOB cmakeFiles "${CMAKE_CURRENT_SOURCE_DIR}/*.cmake")

INSTALL( FILES ${cmakeFiles} DESTINATION cmake )

