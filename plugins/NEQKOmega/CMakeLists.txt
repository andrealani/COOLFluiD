LIST ( APPEND NEQKOmega_files
Euler2DNEQKOmegaRhoivtToCons.cxx
Euler2DNEQKOmegaRhoivtToCons.hh
EulerNEQKOmegaConsVarSet.cxx
EulerNEQKOmegaVarSet.cxx
NEQKOmega.hh
NavierStokesNEQKOmega.ci
NavierStokesNEQKOmega.cxx
NavierStokesNEQKOmega.hh
NavierStokesNEQKOmegaRhoivt.ci
NavierStokesNEQKOmegaRhoivt.cxx
NavierStokesNEQKOmegaRhoivt.hh
NavierStokesNEQKOmegaVarSet.ci
NavierStokesNEQKOmegaVarSet.hh
)

LIST ( APPEND NEQKOmega_cflibs NavierStokes NEQ KOmega )
LIST ( APPEND NEQKOmega_requires_mods NavierStokes NEQ KOmega )

CF_ADD_PLUGIN_LIBRARY ( NEQKOmega )
CF_WARN_ORPHAN_FILES()
ADD_SUBDIRECTORY ( testcases )
