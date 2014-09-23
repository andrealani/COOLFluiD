IF ( CF_HAVE_GET )

  LIST ( APPEND GETModel_files
  GETHeaders.hh
  GETModel.hh
  GET_DM.hh
  GET_DM.cxx
  )

CF_DEBUG_VAR( GET_INCLUDE_DIRS )

  LIST ( APPEND GETModel_includedirs ${GET_INCLUDE_DIRS} )
  LIST ( APPEND GETModel_libs ${GET_LIBRARIES})

  CF_ADD_PLUGIN_LIBRARY ( GETModel )

  LIST ( APPEND test_get_dm_files test_get_dm.cxx )
  LIST ( APPEND test_get_dm_libs ${CF_KERNEL_LIBS} GETModel ${GET_LIBRARIES} )

CF_DEBUG_VAR( test_get_dm_libs )

  CF_ADD_PLUGIN_APP ( test_get_dm )

  CF_WARN_ORPHAN_FILES()

ENDIF ( CF_HAVE_GET )

