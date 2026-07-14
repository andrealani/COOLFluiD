IF ( NOT CF_HAVE_SINGLE_EXEC )
SET ( MYLIBNAME "EntropyMHD" )
ELSE ()
SET ( MYLIBNAME "coolfluid-solver" )
ENDIF()

LIST ( APPEND EntropyMHD_files
EntropyMHD.hh
EntropyMHDTerm.hh
EntropyMHDTerm.cxx
EntropyMHD3D.hh
EntropyMHD3D.ci
EntropyMHD3D.cxx
EntropyMHD3DVarSet.hh
EntropyMHD3DVarSet.cxx
EntropyMHD3DCons.hh
EntropyMHD3DCons.cxx
EntropyMHD3DPrim.hh
EntropyMHD3DPrim.cxx
EntropyMHD3DConsToPrim.hh
EntropyMHD3DConsToPrim.cxx
EntropyMHD3DPrimToCons.hh
EntropyMHD3DPrimToCons.cxx
EntropyMHD3DConsToPrimInPrim.hh
EntropyMHD3DConsToPrimInPrim.cxx
)

IF ( NOT CF_HAVE_SINGLE_EXEC )
LIST ( APPEND EntropyMHD_cflibs Framework )
CF_ADD_PLUGIN_LIBRARY ( EntropyMHD )
ELSE()
FOREACH (AFILE ${EntropyMHD_files})
LIST(APPEND coolfluid-solver_files ../../plugins/EntropyMHD/${AFILE})
ENDFOREACH(AFILE)
LIST ( APPEND coolfluid-solver_cflibs Framework )
ENDIF()

CF_WARN_ORPHAN_FILES()
