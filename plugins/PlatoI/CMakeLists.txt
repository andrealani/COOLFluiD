IF ( CF_ENABLE_PLATO AND NOT CF_ENABLE_LONG )
   LIST ( APPEND PlatoI_files
   Plato.hh
   PlatoLibrary.hh
   PlatoLibrary.cxx
   )
   LIST ( APPEND PlatoI_includedirs ${PLATO_INCLUDE_DIR} )
   LIST ( APPEND PlatoI_libs ${PLATO_LIBRARY} )
   LIST ( APPEND PlatoI_cflibs Framework )		
   CF_ADD_PLUGIN_LIBRARY ( PlatoI )

   CF_WARN_ORPHAN_FILES()
ENDIF()
