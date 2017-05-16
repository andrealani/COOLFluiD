IF ( NOT CF_HAVE_SINGLE_EXEC )
SET ( MYLIBNAME "Maxwell" )
ELSE ()
SET ( MYLIBNAME "coolfluid-solver" )
ENDIF()

LIST ( APPEND Maxwell_files
Maxwell.hh
ConvMaxwellTerm.cxx
ConvMaxwellTerm.hh
DiffMaxwellTerm.cxx
DiffMaxwellTerm.hh
Maxwell2DAdimCons.cxx
Maxwell2DAdimCons.hh
Maxwell2DCons.cxx
Maxwell2DCons.hh
Maxwell2DProjectionAdimVarSet.cxx
Maxwell2DProjectionAdimVarSet.hh
Maxwell2DProjectionVarSet.cxx
Maxwell2DProjectionVarSet.hh
Maxwell2DProjectionAdimCons.cxx
Maxwell2DProjectionAdimCons.hh
Maxwell2DProjectionCons.cxx
Maxwell2DProjectionCons.hh
Maxwell2DVarSet.cxx
Maxwell2DVarSet.hh
Maxwell2DVarSetAdim.cxx
Maxwell2DVarSetAdim.hh
Maxwell3DAdimCons.cxx
Maxwell3DAdimCons.hh
Maxwell3DCons.cxx
Maxwell3DCons.hh
Maxwell3DProjectionAdimVarSet.cxx
Maxwell3DProjectionAdimVarSet.hh
Maxwell3DProjectionCons.cxx
Maxwell3DProjectionCons.hh
Maxwell3DProjectionVarSet.cxx
Maxwell3DProjectionVarSet.hh
Maxwell3DVarSet.cxx
Maxwell3DVarSet.hh
Maxwell3DVarSetAdim.cxx
Maxwell3DVarSetAdim.hh
MaxwellModel.ci
MaxwellModel.cxx
MaxwellModel.hh
MaxwellModelAdim.ci
MaxwellModelAdim.cxx
MaxwellModelAdim.hh
MaxwellProjectionAdim.ci
MaxwellProjectionAdim.cxx
MaxwellProjectionAdim.hh
MaxwellProjection.ci
MaxwellProjection.cxx
MaxwellProjection.hh
MaxwellProjectionAdimTerm.cxx
MaxwellProjectionAdimTerm.hh
MaxwellProjectionTerm.cxx
MaxwellProjectionTerm.hh
MaxwellAdimTerm.cxx
MaxwellAdimTerm.hh
MaxwellVarSet.hh
)

IF (CF_HAVE_CUDA)
   LIST ( APPEND Maxwell_files
Maxwell2DProjectionConsT.hh	
Maxwell2DProjectionVarSetT.hh
Maxwell3DProjectionConsT.hh
Maxwell3DProjectionVarSetT.hh
)
ENDIF()

IF ( NOT CF_HAVE_SINGLE_EXEC )
LIST ( APPEND Maxwell_cflibs Framework )
CF_ADD_PLUGIN_LIBRARY ( Maxwell )
ELSE()
 FOREACH (AFILE ${Maxwell_files} )
 LIST(APPEND coolfluid-solver_files ../../plugins/Maxwell/${AFILE} )
 ENDFOREACH(AFILE)
 CF_CACHE_LIST_APPEND ( ${MYLIBNAME}_files  ${coolfluid-solver_files} )
ENDIF()

CF_WARN_ORPHAN_FILES()
ADD_SUBDIRECTORY ( testcases )
