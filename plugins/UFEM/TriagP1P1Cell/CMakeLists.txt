LIST ( APPEND UFEMTriagP1P1Cell_files
CellProps.hh
CellProps.cxx
UFEMTriagP1P1Cell.hh
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
Mitrem.hh
Mitrem.cxx
MitremBin.hh
MitremBin.cxx
)

LIST ( APPEND OPTIONAL_dirfiles 
  LES.cxx
  LES.hh
  PLaSCoupling.cxx
  PLaSCoupling.hh
)

IF ( CF_HAVE_PLAS )
  LIST ( APPEND UFEMTriagP1P1Cell_includedirs ${PLAS_INCLUDE_DIR} )
	LIST ( APPEND UFEMTriagP1P1Cell_files  PLaSCoupling.cxx PLaSCoupling.hh )
  LIST ( APPEND UFEMTriagP1P1Cell_cflibs PLaSTracking ${PLAS_LIBRARY} )
ENDIF()

# LES extension
FIND_PATH(LESMODELS_INCLUDE_DIR les_interface_v4.h ${LESMODELS_HOME}/interface NO_DEFAULT_PATH)
FIND_LIBRARY(LESMODELS_LIBRARIES vkiles ${LESMODELS_HOME}/vkiles/lib NO_DEFAULT_PATH)
IF ( LESMODELS_INCLUDE_DIR AND LESMODELS_LIBRARIES )
  LIST ( APPEND UFEMTriagP1P1Cell_includedirs ${LESMODELS_INCLUDE_DIR} )
  LIST ( APPEND UFEMTriagP1P1Cell_libs ${LESMODELS_LIBRARIES})
  LIST ( APPEND UFEMTriagP1P1Cell_files LES.cxx  LES.hh )
ENDIF()

LIST ( APPEND UFEMTriagP1P1Cell_cflibs UFEM )

CF_ADD_PLUGIN_LIBRARY ( UFEMTriagP1P1Cell )

CF_WARN_ORPHAN_FILES()


