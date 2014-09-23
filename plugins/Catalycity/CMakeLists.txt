LIST ( APPEND Catalycity_files
Catalycity.hh
BasicCatalycityModel.cxx
BasicCatalycityModel.hh
)

LIST ( APPEND Catalycity_cflibs Framework )
CF_ADD_PLUGIN_LIBRARY ( Catalycity )

##################################################################
CF_WARN_ORPHAN_FILES()
