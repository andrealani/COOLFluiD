#############################################################################
### UFEM library

LIST(APPEND UFEM_files
UFEM.hh
UFEMSolver.hh
UFEMSolver.cxx
UFEMSolverData.hh
UFEMSolverData.cxx
Solve.hh
Solve.cxx
StdSetup.hh
StdSetup.cxx
StdUnSetup.hh
StdUnSetup.cxx
StdClean.hh
StdClean.cxx
InitState.hh
InitState.cxx
RestartState.hh
RestartState.cxx
UFEM_MeshDataBuilder.hh
UFEM_MeshDataBuilder.cxx
BaseBC.hh
BaseBC.cxx
DirichletBC.hh
DirichletBC.cxx
NeumannBC.hh
NeumannBC.cxx
CopyFromTrsBC.hh
CopyFromTrsBC.cxx
#StandardKEpsilonWallBC.hh
#StandardKEpsilonWallBC.cxx
UFEMTerm.hh
UFEMTerm.cxx
ElemProps.hh
ElemProps.cxx
ElemAssembler.hh
AssemblyData.hh
)

LIST ( APPEND UFEM_cflibs Framework ShapeFunctions )

CF_ADD_PLUGIN_LIBRARY ( UFEM )

#############################################################################
LIST ( APPEND OPTIONAL_dirfiles 
  SetupPLaS.hh 
  SetupPLaS.cxx 
  StandardKEpsilonWallBC.cxx
  StandardKEpsilonWallBC.hh
)

### element dependent libs

ADD_SUBDIRECTORY ( QuadP1P1Cell )
ADD_SUBDIRECTORY ( TriagP1P1Cell )
ADD_SUBDIRECTORY ( TetraP1P1Cell )

#############################################################################
### UFEM PLaS library

LIST ( APPEND OPTIONAL_dirfiles 
  SetupPLaS.hh 
  SetupPLaS.cxx 
)

IF (CF_HAVE_PLAS)
  
  LIST ( APPEND UFEMPLaS_files SetupPLaS.hh SetupPLaS.cxx )
  LIST ( APPEND UFEMPLaS_includedirs ${PLAS_INCLUDE_DIR} )
  
  LIST ( APPEND UFEMPLaS_libs ${PLAS_LIBRARY} )
  LIST ( APPEND UFEMPLaS_cflibs UFEM UFEMTetraP1P1Cell UFEMTriagP1P1Cell UFEMQuadP1P1Cell)
  
  CF_ADD_PLUGIN_LIBRARY ( UFEMPLaS )

ENDIF (CF_HAVE_PLAS)

#############################################################################
### testcases

#ADD_SUBDIRECTORY ( testcases )

OPTION ( CF_BUILD_UFEM_UTESTS "Build the UFEM unit tests, used to compare with CF3" OFF )
mark_as_advanced(CF_BUILD_UFEM_UTESTS)
IF(CF_BUILD_UFEM_UTESTS)
  add_subdirectory(uTests)
ENDIF()

CF_WARN_ORPHAN_FILES()
