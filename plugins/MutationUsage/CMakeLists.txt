IF( NOT CF_SKIP_FORTRAN )

LIST ( APPEND MutationUsage_files
MutationUsage.hh
main.cxx
ComputeLTEComposition.hh
ComputeTransportCoefs.hh
ComputeRhoHTOverAt.cxx
ComputeRhoHTOverAt.hh
ComputeWithMutation.hh
ComputeTransportCoefs.cxx
ComputeLTEComposition.cxx
)

LIST ( APPEND MutationUsage_libs MutationI Mutation ${CF_KERNEL_LIBS} ${CF_Boost_LIBRARIES}   )
LIST ( APPEND MutationUsage_clibs Framework Common MutationI Mutation  )
LIST ( APPEND MutationUsage_requires_mods Common MutationI Mutation )
CF_ADD_PLUGIN_APP ( MutationUsage )

CF_WARN_ORPHAN_FILES()

ENDIF( NOT CF_SKIP_FORTRAN )
