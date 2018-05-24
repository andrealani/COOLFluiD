IF ( CF_ENABLE_MUTATIONPP )
   LIST ( APPEND MutationppI_files
   MutationLibrarypp.cxx
   MutationLibrarypp.hh
   Mutationpp.hh
   )

IF ( CF_ENABLE_MUTATIONPP_DEBUG )
LIST ( APPEND MutationppI_files
   MutationLibraryppDebug.hh
   MutationLibraryppDebug.cxx
   )
ENDIF()

   LIST ( APPEND MutationppI_includedirs ${MUTATIONPP_INCLUDE_DIR} )
   LIST ( APPEND MutationppI_libs ${MUTATIONPP_LIBRARY} )
   LIST ( APPEND MutationppI_cflibs Framework )

IF ( CF_ENABLE_MUTATIONPP_DEBUG )
LIST ( APPEND MutationppI_cflibs Mutation2OLD Mutation2OLDI )
ENDIF()
		
   CF_ADD_PLUGIN_LIBRARY ( MutationppI )

   CF_WARN_ORPHAN_FILES()
ENDIF()
