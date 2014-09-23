LIST ( APPEND LESvki_files
SmagorinskyPhysicalModel.ci
SmagorinskyPhysicalModel.cxx
SmagorinskyPhysicalModel.hh
Smagorinsky2DVarSet.hh
Smagorinsky2DVarSet.cxx
Smagorinsky2DCons.hh
Smagorinsky2DCons.cxx
Smagorinsky3DVarSet.hh
Smagorinsky3DVarSet.cxx
Smagorinsky3DCons.hh
Smagorinsky3DCons.cxx
WALESPhysicalModel.ci
WALESPhysicalModel.cxx
WALESPhysicalModel.hh
WALES2DVarSet.hh
WALES2DVarSet.cxx
WALES2DCons.hh
WALES2DCons.cxx
WALES3DVarSet.hh
WALES3DVarSet.cxx
WALES3DCons.hh
WALES3DCons.cxx
ClarkPhysicalModel.ci
ClarkPhysicalModel.cxx
ClarkPhysicalModel.hh
Clark2DVarSet.hh
Clark2DVarSet.cxx
Clark2DCons.hh
Clark2DCons.cxx
GradientPhysicalModel.ci
GradientPhysicalModel.cxx
GradientPhysicalModel.hh
Gradient2DVarSet.hh
Gradient2DVarSet.cxx
Gradient2DCons.hh
Gradient2DCons.cxx
LESvki.hh
)
LIST ( APPEND LESvki_cflibs Framework NavierStokes )
LIST ( APPEND LESvki_requires_mods NavierStokes )
CF_ADD_PLUGIN_LIBRARY ( LESvki )
CF_WARN_ORPHAN_FILES()
