LIST ( APPEND PhysicalModelDummy_files
CouplingModelDummy.cxx
CouplingModelDummy.hh
CouplingModelDummySendToRecv.cxx
CouplingModelDummySendToRecv.hh
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
