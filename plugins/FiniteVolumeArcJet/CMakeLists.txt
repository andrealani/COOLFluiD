LIST ( APPEND FiniteVolumeArcJet_files
ArcJetST.cxx
ArcJetST.hh
ArcJetInsulatedWall.ci
ArcJetInsulatedWall.cxx
ArcJetInsulatedWall.hh
ArcJetPhiBCTurb.cxx
ArcJetPhiElectrode.ci
ArcJetPhiElectrode.cxx
ArcJetPhiElectrode.hh
ArcJetPhiInsulatedWall.ci
ArcJetPhiInsulatedWall.cxx
ArcJetPhiInsulatedWall.hh
ArcJetPhiOutlet.ci
ArcJetPhiOutlet.cxx
ArcJetPhiOutlet.hh
ArcJetPhiST.cxx
ArcJetPhiST.hh
CurrentProps.cxx
CurrentProps.hh
FiniteVolumeArcJet.hh
)

LIST ( APPEND FiniteVolumeArcJet_cflibs NavierStokes LTE SA KOmega SALTE ArcJet FiniteVolume FiniteVolumeNavierStokes FiniteVolumeLTE FiniteVolumeSA FiniteVolumeKOmega )
LIST ( APPEND FiniteVolumeArcJet_requires_mods NavierStokes LTE SA KOmega ArcJet FiniteVolume FiniteVolumeNavierStokes FiniteVolumeTurb )

CF_ADD_PLUGIN_LIBRARY ( FiniteVolumeArcJet )

CF_WARN_ORPHAN_FILES()
