IF ( NOT CF_HAVE_SINGLE_EXEC )
LIST ( APPEND AnalyticalModel_files
AnalyticalModel.hh
AnalyticalDM.hh
AnalyticalDM.cxx
)

LIST ( APPEND AnalyticalModel_cflibs Framework )
CF_ADD_PLUGIN_LIBRARY ( AnalyticalModel )

LIST ( APPEND test_analyticaldm_files test_analytical_model.cxx )
LIST ( APPEND test_analyticaldm_libs ${CF_KERNEL_LIBS}  AnalyticalModel  ${CF_Boost_LIBRARIES} )

CONFIGURE_FILE ( config.conf.in ${CMAKE_CURRENT_BINARY_DIR}/config.conf @ONLY )

CF_ADD_PLUGIN_APP ( test_analyticaldm )
ENDIF()
CF_WARN_ORPHAN_FILES()
