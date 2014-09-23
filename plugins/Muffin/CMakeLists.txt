
LIST (APPEND Muffin_files
  MuffinModule.hh
  BC.cxx                     BC.hh
  BCElectroosmoticWall.cxx   BCElectroosmoticWall.hh
  BCFile.cxx                 BCFile.hh
  BCFixPressure.cxx          BCFixPressure.hh
  BCFixVelocity.cxx          BCFixVelocity.hh
  BCFunction.cxx             BCFunction.hh
  BCMassFlux.cxx             BCMassFlux.hh
  BCWall.cxx                 BCWall.hh
  CC.cxx                     CC.hh
  ComMultipleStates.cxx      ComMultipleStates.hh
  ComSave.cxx                ComSave.hh
  ComSingleStates.cxx        ComSingleStates.hh
  Element.cxx                Element.hh
  Loop.cxx                   Loop.hh
  MuffinMeshDataBuilder.cxx  MuffinMeshDataBuilder.hh
  Muffin.cxx                 Muffin.hh
  MuffinData.cxx             MuffinData.hh
  StdLink.cxx                StdLink.hh
  StdSetup.cxx               StdSetup.hh
  StdSolve.cxx               StdSolve.hh
  StdUnSetup.cxx             StdUnSetup.hh
  System.cxx                 System.hh
  SystemFlow.cxx             SystemFlow.hh
  SystemTemp.cxx             SystemTemp.hh
  SystemTurb.cxx             SystemTurb.hh
)
LIST (APPEND OPTIONAL_dirfiles
  BCBulk.cxx                 BCBulk.hh
  BCElectrode.cxx            BCElectrode.hh
  CCBuoyancyFromMITReM.cxx   CCBuoyancyFromMITReM.hh
  ComUpdateGasRate.cxx       ComUpdateGasRate.hh
  LoopBalanceCurrent.cxx     LoopBalanceCurrent.hh
  SystemMITReM.cxx           SystemMITReM.hh

  CCMomentumFromBubbles.cxx  CCMomentumFromBubbles.hh
  SystemPLaS.cxx             SystemPLaS.hh

  SystemMITReMPLaS.cxx       SystemMITReMPLaS.hh
)

LIST ( APPEND Muffin_cflibs Framework )
CF_ADD_PLUGIN_LIBRARY (Muffin)


IF (CF_HAVE_MITREM)
  LIST (APPEND MuffinMITReM_files
    BCBulk.cxx                BCBulk.hh
    BCElectrode.cxx           BCElectrode.hh
    CCBuoyancyFromMITReM.cxx  CCBuoyancyFromMITReM.hh
    ComUpdateGasRate.cxx      ComUpdateGasRate.hh
    LoopBalanceCurrent.cxx    LoopBalanceCurrent.hh
    SystemMITReM.cxx          SystemMITReM.hh
  )
  LIST (APPEND MuffinMITReM_libs        ${MITREM_LIBRARY})
  LIST (APPEND MuffinMITReM_includedirs ${MITREM_INCLUDE_DIR})
  CF_ADD_PLUGIN_LIBRARY (MuffinMITReM)
ENDIF (CF_HAVE_MITREM)


IF (CF_HAVE_PLAS)
  LIST (APPEND MuffinPLaS_files
    CCMomentumFromBubbles.cxx  CCMomentumFromBubbles.hh
    SystemPLaS.cxx             SystemPLaS.hh
  )
  LIST (APPEND MuffinPLaS_includedirs ${PLAS_INCLUDE_DIR})
  CF_ADD_PLUGIN_LIBRARY (MuffinPLaS)
ENDIF (CF_HAVE_PLAS)


IF (CF_HAVE_MITREM AND CF_HAVE_PLAS)
  SET (MuffinMITReMPLaS_files
    SystemMITReMPLaS.cxx  SystemMITReMPLaS.hh
  )
  LIST (APPEND MuffinMITReMPLaS_libs        ${MITREM_LIBRARY})
  LIST (APPEND MuffinMITReMPLaS_includedirs ${MITREM_INCLUDE_DIR} ${PLAS_INCLUDE_DIR})
  CF_ADD_PLUGIN_LIBRARY (MuffinMITReMPLaS)
ENDIF (CF_HAVE_MITREM AND CF_HAVE_PLAS)


CF_WARN_ORPHAN_FILES()


ADD_SUBDIRECTORY ( testcases )

