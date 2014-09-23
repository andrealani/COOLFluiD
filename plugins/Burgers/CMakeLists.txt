LIST ( APPEND Burgers_files
Burgers.hh
Burgers2DLinearPrim.cxx
Burgers2DLinearPrim.hh
Burgers2DPrim.cxx
Burgers2DPrim.hh
Burgers2DVarSet.cxx
Burgers2DVarSet.hh
BurgersPhysicalModel.ci
BurgersPhysicalModel.cxx
BurgersPhysicalModel.hh
BurgersTerm.cxx
BurgersTerm.hh
)

LIST ( APPEND Burgers_cflibs Framework )

CF_ADD_PLUGIN_LIBRARY ( Burgers )
CF_WARN_ORPHAN_FILES()
ADD_SUBDIRECTORY ( testcases )
