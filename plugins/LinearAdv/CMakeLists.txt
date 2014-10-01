LIST ( APPEND LinearAdv_files
LinearAdv.hh
LinearAdv2DLinearPrim.cxx
LinearAdv2DLinearPrim.hh
LinearAdv2DPrim.cxx
LinearAdv2DPrim.hh
LinearAdv2DVarSet.cxx
LinearAdv2DVarSet.hh
LinearAdv3DLinearPrim.cxx
LinearAdv3DLinearPrim.hh
LinearAdv3DPrim.cxx
LinearAdv3DPrim.hh
LinearAdv3DVarSet.cxx
LinearAdv3DVarSet.hh
LinearAdvData.hh
LinearAdvPhysicalModel.ci
LinearAdvPhysicalModel.cxx
LinearAdvPhysicalModel.hh
LinearAdvTerm.cxx
LinearAdvTerm.hh
AdvectionDiffusionData.hh
ADTerm.cxx
ADTerm.hh
AdvectionDiffusion2DPrim.cxx
AdvectionDiffusion2DPrim.hh
AdvectionDiffusion2DVarSet.cxx
AdvectionDiffusion2DVarSet.hh
AdvectionDiffusion3DPrim.cxx
AdvectionDiffusion3DPrim.hh
AdvectionDiffusion3DVarSet.cxx
AdvectionDiffusion3DVarSet.hh
AdvectionDiffusionPhysicalModel.ci
AdvectionDiffusionPhysicalModel.cxx
AdvectionDiffusionPhysicalModel.hh
AdvectionDiffusionVarSet.hh
)

LIST ( APPEND LinearAdv_libs Framework )

CF_ADD_PLUGIN_LIBRARY ( LinearAdv )
CF_WARN_ORPHAN_FILES()
#ADD_SUBDIRECTORY ( testcases )
