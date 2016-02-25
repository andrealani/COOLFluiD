IF( NOT CF_SKIP_FORTRAN AND CF_HAVE_MUTATION2 )

LIST ( APPEND Mutation2I_files
MutationLibrary2.hh
Mutation2.hh
MutationLibrary2.cxx
)

LIST ( APPEND Mutation2I_cflibs Framework Mutation2 )
LIST ( APPEND Mutation2I_requires_mods Mutation2.0 )
CF_ADD_PLUGIN_LIBRARY ( Mutation2I )

CF_WARN_ORPHAN_FILES()

ENDIF( NOT CF_SKIP_FORTRAN AND CF_HAVE_MUTATION2 )
