IF( NOT CF_SKIP_FORTRAN )

LIST ( APPEND Mutation2OLDI_files
MutationLibrary2OLD.hh
Mutation2OLD.hh
MutationLibrary2OLD.cxx
)

LIST ( APPEND Mutation2OLDI_cflibs Framework Mutation2OLD )
LIST ( APPEND Mutation2OLDI_requires_mods Mutation2.0.0 )
CF_ADD_PLUGIN_LIBRARY ( Mutation2OLDI )
CF_WARN_ORPHAN_FILES()

ENDIF( NOT CF_SKIP_FORTRAN )
