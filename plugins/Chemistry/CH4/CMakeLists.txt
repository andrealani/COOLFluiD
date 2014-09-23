LIST ( APPEND CH4_files
CH4.hh
ChemCH43D.cxx
ChemCH43D.hh
ChemCH43DDiffusiveVarSet.hh
ChemCH43DPrim.cxx
ChemCH43DPrim.hh
ChemCH43DPrimArrhenius.cxx
ChemCH43DPrimArrhenius.hh
ChemCH43DPrimDiffusive.cxx
ChemCH43DPrimDiffusive.hh
ChemCH43DSourceVarSet.hh
ChemCH4PhysicalModel.cxx
ChemCH4PhysicalModel.hh
)

LIST ( APPEND CH4_cflibs Framework )
CF_ADD_PLUGIN_LIBRARY ( CH4 )

CF_WARN_ORPHAN_FILES()
