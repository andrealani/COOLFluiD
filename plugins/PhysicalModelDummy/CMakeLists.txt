IF ( NOT CF_HAVE_SINGLE_EXEC )
SET ( MYLIBNAME "PhysicalModelDummy" )
ELSE ()     
SET ( MYLIBNAME "coolfluid-solver" )
ENDIF()

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

IF ( NOT CF_HAVE_SINGLE_EXEC )
LIST ( APPEND PhysicalModelDummy_cflibs Framework )
CF_ADD_PLUGIN_LIBRARY ( PhysicalModelDummy )
ELSE()
 FOREACH (AFILE ${PhysicalModelDummy_files} )
 LIST(APPEND coolfluid-solver_files ../../plugins/PhysicalModelDummy/${AFILE} )
 ENDFOREACH(AFILE)
 CF_CACHE_LIST_APPEND ( ${MYLIBNAME}_files  ${coolfluid-solver_files} )
ENDIF()

CF_WARN_ORPHAN_FILES()

#ADD_SUBDIRECTORY ( testcases )

