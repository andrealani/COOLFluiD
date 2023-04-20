IF ( NOT CF_HAVE_SINGLE_EXEC )
SET ( MYLIBNAME "MHD" )
ELSE ()
SET ( MYLIBNAME "coolfluid-solver" )
ENDIF()

LIST ( APPEND MHD_files
MHD.hh
MHD2DCons.cxx
MHD2DCons.hh
MHD2DConsToPrim.cxx
MHD2DConsToPrim.hh
MHD2DConsToPrimInRef.cxx
MHD2DConsToPrimInRef.hh
MHD2DLinearCons.cxx
MHD2DLinearCons.hh
MHD2DLinearPrim.cxx
MHD2DLinearPrim.hh
MHD2DPrim.cxx
MHD2DPrim.hh
MHD2DPrimToCons.cxx
MHD2DPrimToCons.hh
MHD2DPrimToConsInRef.cxx
MHD2DPrimToConsInRef.hh
MHD2DProjectionCons.cxx
MHD2DProjectionCons.hh
MHD2DProjectionConsToPrim.cxx
MHD2DProjectionConsToPrim.hh
MHD2DProjectionConsToPrimInRef.cxx
MHD2DProjectionConsToPrimInRef.hh
MHD2DProjectionLinearCons.cxx
MHD2DProjectionLinearCons.hh
MHD2DProjectionLinearPrim.cxx
MHD2DProjectionLinearPrim.hh
MHD2DProjectionPrim.cxx
MHD2DProjectionPrim.hh
MHD2DProjectionPrimToCons.cxx
MHD2DProjectionPrimToCons.hh
MHD2DProjectionPrimToConsInRef.cxx
MHD2DProjectionPrimToConsInRef.hh
MHD2DProjectionVarSet.cxx
MHD2DProjectionVarSet.hh
MHD2DVarSet.cxx
MHD2DVarSet.hh
MHD3DCons.cxx
MHD3DCons.hh
MHD3DConsToPrim.cxx
MHD3DConsToPrim.hh
MHD3DConsToPrimInRef.cxx
MHD3DConsToPrimInRef.hh
MHD3DLinearCons.cxx
MHD3DLinearCons.hh
MHD3DLinearPrim.cxx
MHD3DLinearPrim.hh
MHD3DPrim.cxx
MHD3DPrim.hh
MHD3DPrimToCons.cxx
MHD3DPrimToCons.hh
MHD3DPrimToConsInRef.cxx
MHD3DPrimToConsInRef.hh
MHD3DProjectionCons.cxx
MHD3DProjectionCons.hh
MHD3DProjectionConsToPrim.cxx
MHD3DProjectionConsToPrim.hh
MHD3DProjectionDiffVarSet.cxx
MHD3DProjectionDiffVarSet.hh
MHD3DProjectionDiffPrim.cxx
MHD3DProjectionDiffPrim.hh
MHD3DProjectionPolytropicConsToPrim.cxx
MHD3DProjectionPolytropicConsToPrim.hh
MHD3DProjectionConsToPrimInPrim.cxx
MHD3DProjectionConsToPrimInPrim.hh
MHD3DProjectionConsToPrimInRef.cxx
MHD3DProjectionConsToPrimInRef.hh
MHD3DProjectionLinearCons.cxx
MHD3DProjectionLinearCons.hh
MHD3DProjectionLinearPrim.cxx
MHD3DProjectionLinearPrim.hh
MHD3DProjectionPrim.cxx
MHD3DProjectionPrim.hh
MHD3DProjectionPrimToCons.cxx
MHD3DProjectionPrimToCons.hh
MHD3DProjectionPolytropicCons.cxx
MHD3DProjectionPolytropicCons.hh
MHD3DProjectionPolytropicPrim.cxx
MHD3DProjectionPolytropicPrim.hh
MHD3DProjectionPolytropicPrimToCons.cxx
MHD3DProjectionPolytropicPrimToCons.hh
MHD3DProjectionPolytropicVarSet.hh
MHD3DProjectionPolytropicVarSet.cxx
MHD3DProjectionPrimToConsInRef.cxx
MHD3DProjectionPrimToConsInRef.hh
MHD3DProjectionVarSet.hh
MHD3DProjectionVarSet.cxx
MHD3DVarSet.cxx
MHD3DVarSet.hh
MHDPhysicalModel.ci
MHDPhysicalModel.cxx
MHDPhysicalModel.hh
MHDProjection.ci
MHDProjection.cxx
MHDProjection.hh
MHDProjectionDiffPhysicalModel.ci
MHDProjectionDiffPhysicalModel.cxx
MHDProjectionDiffPhysicalModel.hh
MHDProjectionDiffTerm.cxx
MHDProjectionDiffTerm.hh
MHDProjectionDiffVarSet.cxx
MHDProjectionDiffVarSet.hh
MHDProjectionPolytropic.ci
MHDProjectionPolytropic.cxx
MHDProjectionPolytropic.hh
MHDProjectionTerm.cxx
MHDProjectionTerm.hh
MHDProjectionPolytropicTerm.cxx
MHDProjectionPolytropicTerm.hh
MHDTerm.cxx
MHDTerm.hh
MHD3DProjectionEpsCons.cxx
MHD3DProjectionEpsCons.hh
MHD3DProjectionEpsConsToPrim.cxx
MHD3DProjectionEpsConsToPrim.hh
MHD3DProjectionEpsPrim.cxx
MHD3DProjectionEpsPrim.hh
MHD3DProjectionEpsPrimToCons.cxx
MHD3DProjectionEpsPrimToCons.hh
MHD3DProjectionEpsVarSet.cxx
MHD3DProjectionEpsVarSet.hh
MHDProjectionEpsDiffPhysicalModel.ci
MHDProjectionEpsDiffPhysicalModel.cxx
MHDProjectionEpsDiffPhysicalModel.hh
MHD3DProjectionEpsDiffPrim.cxx
MHD3DProjectionEpsDiffPrim.hh
MHDProjectionEpsTerm.cxx
MHDProjectionEpsTerm.hh
)

IF (CF_HAVE_CUDA)
   LIST ( APPEND MHD_files
MHD2DProjectionConsT.hh
MHD2DProjectionPrimT.hh
MHD2DProjectionVarSetT.hh
MHD3DProjectionConsT.hh
MHD3DProjectionPrimT.hh
MHD3DProjectionVarSetT.hh
)
ENDIF()

IF (CF_HAVE_GSL)
   LIST ( APPEND ${MYLIBNAME}_includedirs ${GSL_INCLUDE_DIR} )
   LIST ( APPEND ${MYLIBNAME}_libs ${GSL_LIBRARY} )
ENDIF()

IF ( NOT CF_HAVE_SINGLE_EXEC )
LIST ( APPEND MHD_cflibs Framework )
CF_ADD_PLUGIN_LIBRARY ( MHD )
ELSE()
IF ( CF_BUILD_MHD )
 FOREACH (AFILE ${MHD_files} )
 LIST(APPEND coolfluid-solver_files ../../plugins/MHD/${AFILE} )
 ENDFOREACH(AFILE)
 CF_CACHE_LIST_APPEND ( ${MYLIBNAME}_files  ${coolfluid-solver_files} )
ENDIF()
ENDIF()

CF_WARN_ORPHAN_FILES()
ADD_SUBDIRECTORY ( testcases )
