LIST ( APPEND CFmeshExtruder_files
CFmeshExtruder.hh
Extruder2D.cxx
Extruder2D.hh
Extruder2DFVM.cxx
Extruder2DFVM.hh
Extruder2DDGM.cxx
Extruder2DDGM.hh
)

IF (CF_HAVE_MPI)
LIST ( APPEND CFmeshExtruder_files Extruder2DFVMMPI.cxx Extruder2DFVMMPI.hh )
ENDIF()

LIST ( APPEND CFmeshExtruder_cflibs Framework )

CF_ADD_PLUGIN_LIBRARY ( CFmeshExtruder )
CF_WARN_ORPHAN_FILES()
