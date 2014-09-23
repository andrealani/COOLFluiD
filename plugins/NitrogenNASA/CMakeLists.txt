IF( NOT CF_SKIP_FORTRAN )

LIST ( APPEND NitrogenNASA_files
n2initialize.f
n2thermo.f
n2arrhenius.f
n2lookup.f
commondata.cmn
)

CF_ADD_PLUGIN_LIBRARY ( NitrogenNASA )
CF_WARN_ORPHAN_FILES()

ENDIF( NOT CF_SKIP_FORTRAN )

