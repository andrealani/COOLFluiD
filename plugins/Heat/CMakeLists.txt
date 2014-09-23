LIST ( APPEND Heat_files
ConvHeatTerm.cxx
ConvHeatTerm.hh
DiffHeatTerm.cxx
DiffHeatTerm.hh
Heat.hh
Heat2D.cxx
Heat2D.hh
Heat2DDiffusivePrim.cxx
Heat2DDiffusivePrim.hh
Heat2DDiffusiveVarSet.hh
Heat2DPrim.cxx
Heat2DPrim.hh
Heat2DSourceTConst.cxx
Heat2DSourceTConst.hh
Heat2DSourceTDep.cxx
Heat2DSourceTDep.hh
Heat2DSourceVarSet.hh
Heat3D.cxx
Heat3D.hh
Heat3DDiffusivePrim.cxx
Heat3DDiffusivePrim.hh
Heat3DDiffusiveVarSet.hh
Heat3DPrim.cxx
Heat3DPrim.hh
Heat3DSourceTConst.cxx
Heat3DSourceTConst.hh
Heat3DSourceTDep.cxx
Heat3DSourceTDep.hh
Heat3DSourceVarSet.hh
HeatPhysicalModel.cxx
HeatPhysicalModel.hh
)

LIST ( APPEND Heat_cflibs Framework )

CF_ADD_PLUGIN_LIBRARY ( Heat )
CF_WARN_ORPHAN_FILES()
ADD_SUBDIRECTORY ( testcases )
