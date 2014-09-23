LIST ( APPEND ATDModel_files
ATDModel.hh
ATDModelLibrary.hh
ATDModelLibrary.cxx
CombustionModelLibrary.hh
CombustionModelLibrary.cxx
)

LIST ( APPEND ATDModel_cflibs Framework )
CF_ADD_PLUGIN_LIBRARY ( ATDModel )

CF_WARN_ORPHAN_FILES()
