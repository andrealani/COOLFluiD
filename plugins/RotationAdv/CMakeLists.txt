LIST ( APPEND RotationAdv_files
RotationAdv.hh
RotationAdv2DLinearPrim.cxx
RotationAdv2DLinearPrim.hh
RotationAdv2DPrim.cxx
RotationAdv2DPrim.hh
RotationAdv2DVarSet.cxx
RotationAdv2DVarSet.hh
RotationAdv3DLinearPrim.cxx
RotationAdv3DLinearPrim.hh
RotationAdv3DPrim.cxx
RotationAdv3DPrim.hh
RotationAdv3DVarSet.cxx
RotationAdv3DVarSet.hh
RotationAdvPhysicalModel.ci
RotationAdvPhysicalModel.cxx
RotationAdvPhysicalModel.hh
RotationAdvTerm.cxx
RotationAdvTerm.hh
RDTerm.cxx
RDTerm.hh
RotationDiffusion2DPrim.cxx
RotationDiffusion2DPrim.hh
RotationDiffusion2DVarSet.cxx
RotationDiffusion2DVarSet.hh
RotationDiffusionData.hh
RotationDiffusionPhysicalModel.ci
RotationDiffusionPhysicalModel.cxx
RotationDiffusionPhysicalModel.hh
RotationDiffusionVarSet.hh
RotationAdvData.hh
)

LIST ( APPEND RotationAdv_cflibs Framework )

CF_ADD_PLUGIN_LIBRARY ( RotationAdv )
CF_WARN_ORPHAN_FILES()
ADD_SUBDIRECTORY ( testcases )
