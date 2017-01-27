LIST ( APPEND RadiativeTransfer_files
RadiativeTransfer.hh
RadiativeTransferOutputFVMCC.cxx
RadiativeTransferOutputFVMCC.hh
RadiationLibrary/RadiationPhysics.hh
RadiationLibrary/RadiationPhysics.cxx
RadiationLibrary/RadiationPhysicsHandler.hh
RadiationLibrary/RadiationPhysicsHandler.cxx
RadiationLibrary/Radiator.hh
RadiationLibrary/Radiator.cxx
RadiationLibrary/Reflector.hh
RadiationLibrary/Reflector.cxx
RadiationLibrary/Models/Grey/GreyRadiator.hh
RadiationLibrary/Models/Grey/GreyRadiator.cxx
RadiationLibrary/Models/Null/NullRadiator.hh
RadiationLibrary/Models/Null/NullRadiator.cxx
RadiationLibrary/Models/Null/NullReflector.hh
RadiationLibrary/Models/Null/NullReflector.cxx
RadiationLibrary/Models/ArcJet/ArcJetRadiator.hh
RadiationLibrary/Models/ArcJet/ArcJetRadiator.cxx
RadiationLibrary/Models/PARADE/ParadeRadiator.hh
RadiationLibrary/Models/PARADE/ParadeRadiator.cxx
RadiationLibrary/Models/Reflection/DiffuseReflector.hh
RadiationLibrary/Models/Reflection/DiffuseReflector.cxx
RadiationLibrary/Models/Reflection/SpecularReflector.hh
RadiationLibrary/Models/Reflection/SpecularReflector.cxx
Solvers/MonteCarlo/RadiativeTransferMonteCarlo.hh
Solvers/MonteCarlo/RadiativeTransferMonteCarlo.cxx
Solvers/MonteCarlo/RandomNumberGenerator.hh
Solvers/MonteCarlo/RandomNumberGenerator.cxx
Solvers/FiniteVolumeDOM/RadiativeTransferFVDOM.cxx
Solvers/FiniteVolumeDOM/RadiativeTransferFVDOM.hh
PostProcess/PostProcessAverage.cxx
PostProcess/PostProcessAverage.hh
PostProcess/PostProcess.cxx
PostProcess/PostProcess.hh
PostProcess/PostProcessNull.cxx
PostProcess/PostProcessNull.hh
PostProcess/PostProcessModule.hh
)

IF (CF_HAVE_GSL)
   LIST ( APPEND RadiativeTransfer_includedirs ${GSL_INCLUDE_DIR} )
   LIST ( APPEND RadiativeTransfer_libs ${GSL_LIBRARY} ${BLAS_LIBRARY} )
   LIST ( APPEND RadiativeTransfer_files
     PostProcess/PostProcessLSF.cxx
     PostProcess/PostProcessLSF.hh
   )
ENDIF()

IF (CF_HAVE_CUDA)
  LIST ( APPEND RadiativeTransfer_includedirs ${MPI_INCLUDE_DIR} ${CUDA_INCLUDE_DIR} )
  LIST ( APPEND RadiativeTransfer_libs ${CUDA_LIBRARIES} )
  LIST ( APPEND RadiativeTransferCUDA_cflibs FiniteVolume RadiativeTransfer )
  LIST ( APPEND RadiativeTransfer_files
    RadiationLibrary/Models/PARADE/ParadeRadiatorCUDA.cu
    RadiationLibrary/Models/PARADE/ParadeRadiatorCUDA.hh
    Solvers/FiniteVolumeDOM/RadiativeTransferFVDOMCUDA.cu
    Solvers/FiniteVolumeDOM/RadiativeTransferFVDOMCUDA.hh
    )
ENDIF()

LIST ( APPEND RadiativeTransfer_cflibs Framework LagrangianSolver FiniteVolume )
CF_ADD_PLUGIN_LIBRARY ( RadiativeTransfer )
CF_WARN_ORPHAN_FILES()
ADD_SUBDIRECTORY ( testcases )
