LIST ( APPEND BackwardEuler_files
BwdEulerData.hh
BackwardEuler.hh
StdUnSetup.hh
UpdateSol.hh
StdSetup.hh
BwdEuler.hh
UpdateSol.cxx
BwdEulerData.cxx
StdUnSetup.cxx
StdSetup.cxx
BwdEuler.cxx
)

LIST ( APPEND BackwardEuler_cflibs Framework )

CF_ADD_PLUGIN_LIBRARY ( BackwardEuler )

##################################################################

LIST ( APPEND BackwardEulerMHD_files
BackwardEulerMHD.hh
UpdateSolMHD.hh
UpdateSolMHD.cxx
UpdateSolOginoMHD.hh
UpdateSolOginoMHD.cxx
StdSetupOginoMHD.hh
StdSetupOginoMHD.cxx
)

LIST ( APPEND BackwardEulerMHD_requires_mods MHD )
LIST ( APPEND BackwardEulerMHD_cflibs BackwardEuler MHD )


CF_ADD_PLUGIN_LIBRARY ( BackwardEulerMHD )

##################################################################

CF_WARN_ORPHAN_FILES()
