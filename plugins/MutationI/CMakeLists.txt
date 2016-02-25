IF( NOT CF_SKIP_FORTRAN AND CF_HAVE_MUTATION1 )

LIST ( APPEND MutationI_files
MutationLibrary.hh
Mutation.hh
MutationLibrary.cxx
#OnlySigma.hh
#OnlySigma.cxx
)

LIST ( APPEND MutationI_cflibs Framework Common Mutation )
LIST ( APPEND MutationI_requires_mods Mutation )
CF_ADD_PLUGIN_LIBRARY ( MutationI )

CF_WARN_ORPHAN_FILES()

ENDIF( NOT CF_SKIP_FORTRAN AND CF_HAVE_MUTATION1 )
