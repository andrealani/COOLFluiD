LIST ( APPEND NonLinearAdv_files
NonLinearAdv.hh
NonLinearAdv2DLinearPrim.cxx
NonLinearAdv2DLinearPrim.hh
NonLinearAdv2DPrim.cxx
NonLinearAdv2DPrim.hh
NonLinearAdv2DVarSet.cxx
NonLinearAdv2DVarSet.hh
NonLinearAdvPhysicalModel.ci
NonLinearAdvPhysicalModel.cxx
NonLinearAdvPhysicalModel.hh
NonLinearAdvTerm.cxx
NonLinearAdvTerm.hh
)

LIST ( APPEND NonLinearAdv_cflibs Framework )
CF_ADD_PLUGIN_LIBRARY ( NonLinearAdv )
CF_WARN_ORPHAN_FILES()
ADD_SUBDIRECTORY ( testcases )
