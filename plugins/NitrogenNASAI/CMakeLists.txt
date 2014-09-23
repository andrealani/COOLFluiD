IF( NOT CF_SKIP_FORTRAN )

LIST ( APPEND NitrogenNASAI_files
NitrogenNASA.hh
NitrogenNASALibrary.hh
NitrogenNASALibrary.cxx
)

LIST ( APPEND NitrogenNASAI_cflibs Framework NitrogenNASA )
LIST ( APPEND NitrogenNASAI_requires_mods NitrogenNASA )
CF_ADD_PLUGIN_LIBRARY ( NitrogenNASAI )

CF_WARN_ORPHAN_FILES()

ENDIF( NOT CF_SKIP_FORTRAN )
