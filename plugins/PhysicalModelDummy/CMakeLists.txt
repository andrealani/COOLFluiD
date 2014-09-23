LIST ( APPEND PhysicalModelDummy_files
Dummy.hh
DummyPrim.cxx
DummyPrim.hh
DummyTerm.cxx
DummyTerm.hh
DummyVarSet.hh
PhysicalModelDummy.cxx
PhysicalModelDummy.hh
)

LIST ( APPEND PhysicalModelDummy_cflibs Framework )

CF_ADD_PLUGIN_LIBRARY ( PhysicalModelDummy )
CF_WARN_ORPHAN_FILES()
