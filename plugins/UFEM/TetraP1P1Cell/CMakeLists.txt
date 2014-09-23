LIST ( APPEND UFEMTetraP1P1Cell_files
CellProps.hh
CellProps.cxx
UFEMTetraP1P1Cell.hh
NStokesStabilized.hh
NStokesStabilized.cxx
PassiveScalar.hh
PassiveScalar.cxx
HeatConduction.cxx
HeatConduction.hh
StandardKEpsilon.cxx
StandardKEpsilon.hh
SpalartAllmaras.cxx
SpalartAllmaras.hh
Hydrogen.cxx
Hydrogen.hh
)

LIST ( APPEND OPTIONAL_dirfiles 
  LES.cxx
  LES.hh
  PLaSCoupling.cxx
  PLaSCoupling.hh
)

IF ( CF_HAVE_PLAS )
  LIST ( APPEND UFEMTetraP1P1Cell_includedirs ${PLAS_INCLUDE_DIR} )
	LIST ( APPEND UFEMTetraP1P1Cell_files  PLaSCoupling.cxx PLaSCoupling.hh )
  LIST ( APPEND UFEMTetraP1P1Cell_cflibs PLaSTracking ${PLAS_LIBRARY} )
ENDIF()

# LES extension
FIND_PATH(LESMODELS_INCLUDE_DIR les_interface_v4.h ${LESMODELS_HOME}/interface NO_DEFAULT_PATH)
FIND_LIBRARY(LESMODELS_LIBRARIES vkiles ${LESMODELS_HOME}/vkiles/lib NO_DEFAULT_PATH)
IF ( LESMODELS_INCLUDE_DIR AND LESMODELS_LIBRARIES )
  LIST ( APPEND UFEMTetraP1P1Cell_includedirs ${LESMODELS_INCLUDE_DIR} )
  LIST ( APPEND UFEMTetraP1P1Cell_libs ${LESMODELS_LIBRARIES})
  LIST ( APPEND UFEMTetraP1P1Cell_files LES.cxx  LES.hh )
ENDIF()

LIST ( APPEND UFEMTetraP1P1Cell_cflibs UFEM )

CF_ADD_PLUGIN_LIBRARY ( UFEMTetraP1P1Cell )

CF_WARN_ORPHAN_FILES()


