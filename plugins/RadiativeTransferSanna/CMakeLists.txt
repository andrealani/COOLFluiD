LIST ( APPEND RadiativeTransferSanna_files
#MersenneTwister.hh
#MersenneTwister.cxx
ParticleTracking.hh
ParticleTracking.cxx
RadiativeTransferSanna.hh
RadiativeTransferMonteCarlo.hh
RadiativeTransferMonteCarlo.cxx
Slab1DFVMCC.cxx
Slab1DFVMCC.hh
RandomNumberGenerator.hh
RandomNumberGenerator.cxx
)

# Option on MPI
##############################################################################
# TODO: this dependency on MPI should somehow be removed

LIST ( APPEND OPTIONAL_dirfiles
       Slab1DFVMCCMPI.hh
       Slab1DFVMCCMPI.cxx
     )

IF ( CF_HAVE_MPI ) # TEST: works now without PARMETIS ?
  LIST ( APPEND RadiativeTransferSanna_files 
       Slab1DFVMCCMPI.hh
       Slab1DFVMCCMPI.cxx
  )
ENDIF()

LIST ( APPEND RadiativeTransferSanna_requires_mods FiniteVolume FiniteVolumeNavierStokes NavierStokes )
LIST ( APPEND RadiativeTransferSanna_cflibs Framework FiniteVolume FiniteVolumeNavierStokes NavierStokes)

CF_ADD_PLUGIN_LIBRARY ( RadiativeTransferSanna )

##################################################################

ADD_SUBDIRECTORY ( testcases )

##################################################################

CF_WARN_ORPHAN_FILES()
