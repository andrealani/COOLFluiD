LIST ( APPEND FluctSplitNEQ_files
DiffTermCarbuncleFixTCNEQ.cxx
#DiffTermCarbuncleFixTCNEQ_Ec.ci
#DiffTermCarbuncleFixTCNEQ_Ec.hh
ExtraDissActiveSetup.cxx
ExtraDissActiveSetup.hh
ExtraDissActiveUnSetup.cxx
ExtraDissActiveUnSetup.hh
FluctSplitNEQ.hh
NSchemeCSysNEQ.cxx
NSchemeCSysNEQ.hh
SplitStrategyNEQ.hh
SplitStrategyNEQ.ci
SplitStrategyNEQ.cxx
StrongNoSlipWallIsothermalTCNEQ2DImpl.cxx
StrongNoSlipWallIsothermalTCNEQ2DImpl.hh
TCNEQDiffTerm.ci
TCNEQDiffTerm.cxx
TCNEQDiffTerm.hh
TCNEQAxiSourceTerm.hh
TCNEQAxiSourceTerm.ci
TCNEQSourceTerm.ci
TCNEQSourceTerm.cxx
TCNEQSourceTerm.hh
WeakSlipWallTCNEQ2DImpl.cxx
WeakSlipWallTCNEQ2DImpl.hh
WeakSlipWallTCNEQ2DImplCNEQ.cxx
WeakSlipWallTCNEQ2DImplCNEQ.hh
)

LIST ( APPEND FluctSplitNEQ_cflibs NavierStokes NEQ FluctSplit FluctSplitScalar FluctSplitSystem FluctSplitNavierStokes )
LIST ( APPEND FluctSplitNEQ_requires_mods FluctSplit NavierStokes NEQ )
     
CF_ADD_PLUGIN_LIBRARY ( FluctSplitNEQ )
CF_WARN_ORPHAN_FILES()
