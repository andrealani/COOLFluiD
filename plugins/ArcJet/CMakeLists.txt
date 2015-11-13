LIST ( APPEND ArcJet_files
ArcJet.hh
ArcJetConvVarSet.hh
ArcJetConvVarSet.ci
ArcJetConvVarSet.cxx
ArcJetDiffVarSet.hh
ArcJetDiffVarSet.ci
ArcJetDiffVarSet.cxx
ArcJetInductionConvVarSet.hh
ArcJetInductionConvVarSet.ci
ArcJetInductionConvVarSet.cxx
ArcJetInductionTerm.hh
ArcJetInductionTerm.ci
ArcJetInductionDiffVarSet.hh
ArcJetInductionDiffVarSet.ci
ArcJetInductionDiffVarSet.cxx
ArcJetLTEPhysicalModel.hh
ArcJetLTEPhysicalModel.ci
ArcJetLTEPhysicalModel.cxx
ArcJetPvtToCons.hh
ArcJetPvtToCons.ci
ArcJetPvtToCons.cxx
ArcJetLTEPvtToCons.hh
ArcJetLTEPvtToCons.ci
ArcJetLTEPvtToCons.cxx
ArcJetPhysicalModel.hh
ArcJetPhysicalModel.ci
ArcJetPhysicalModel.cxx
ArcJetTerm.hh
ArcJetTerm.ci
EulerArcJetPhysicalModel.hh
EulerArcJetPhysicalModel.ci
)

LIST ( APPEND ArcJet_requires_mods NavierStokes LTE )
LIST ( APPEND ArcJet_cflibs Common MathTools Framework NavierStokes LTE )

CF_ADD_PLUGIN_LIBRARY ( ArcJet )

LIST ( APPEND ArcJetSALTE_files
ArcJetSALTE.hh
ArcJetSALTEConvVarSet.cxx
ArcJetSALTEDiffVarSet.cxx
ArcJetSALTEPhysicalModel.cxx
ArcJetSALTEPvtToCons.cxx
)

LIST ( APPEND ArcJetSALTE_requires_mods NavierStokes LTE ArcJet SA SALTE )
LIST ( APPEND ArcJetSALTE_cflibs Common MathTools Framework NavierStokes LTE SA SALTE )

CF_ADD_PLUGIN_LIBRARY ( ArcJetSALTE )

CF_WARN_ORPHAN_FILES()
#ADD_SUBDIRECTORY ( testcases )
