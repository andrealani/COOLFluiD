IF ( CF_ENABLE_MUTATIONPP )
   LIST ( APPEND MutationppI_files
   MutationLibrarypp.hh
   Mutationpp.hh
   MutationLibrarypp.cxx
   )
   LIST ( APPEND MutationppI_includedirs ${MUTATIONPP_INCLUDE_DIR} )
   LIST ( APPEND MutationppI_libs ${MUTATIONPP_LIBRARY} )
   LIST ( APPEND MutationppI_cflibs Framework )		
   CF_ADD_PLUGIN_LIBRARY ( MutationppI )

   CF_WARN_ORPHAN_FILES()
ENDIF()
