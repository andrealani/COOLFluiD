LIST ( APPEND FiniteVolumeArcJet_files
ArcJetST.cxx
ArcJetST.hh
ArcJetInsulatedWall.ci
ArcJetInsulatedWall.cxx
ArcJetInsulatedWall.hh
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
Radiation.cxx
Radiation.hh
FiniteVolumeArcJet.hh
)

LIST ( APPEND FiniteVolumeArcJet_cflibs NavierStokes LTE ArcJet FiniteVolume FiniteVolumeNavierStokes )
LIST ( APPEND FiniteVolumeArcJet_requires_mods NavierStokes LTE ArcJet FiniteVolume FiniteVolumeNavierStokes )

CF_ADD_PLUGIN_LIBRARY ( FiniteVolumeArcJet )

CF_WARN_ORPHAN_FILES()
