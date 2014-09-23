IF ( CF_HAVE_LESMODELS )

  LIST ( APPEND LESDataProcessing_files
  LESDataProcessing.hh
  LESProcessingData.hh
  LESProcessingData.cxx
  LESProcessingMethod.hh
  LESProcessingMethod.cxx
  LESProcessingCom.hh
  TimeAveraging.hh
  TimeAveraging.cxx
  GradientComputer.hh
  GradientComputer.cxx
  GradientComputerFVMCC.hh
  GradientComputerFVMCC.cxx
  ComputeTurbulenceFunctions.hh
  ComputeTurbulenceFunctions.cxx
  TurbulenceFunction.hh
  TurbulenceFunction.cxx
  Vorticity2D.hh
  Vorticity2D.cxx
  Vorticity3D.hh
  Vorticity3D.cxx
  Qcriterion2D.hh
  Qcriterion2D.cxx
  Qcriterion3D.hh
  Qcriterion3D.cxx
  SGSViscosity.hh
  SGSViscosity.cxx
  SGSKineticEnergy.hh
  SGSKineticEnergy.cxx
  VelocityGradients2D.hh
  VelocityGradients2D.cxx
  VelocityGradients3D.hh
  VelocityGradients3D.cxx
  CFL.hh
  CFL.cxx
  LocalTimeStep.hh
  LocalTimeStep.cxx
  )

  LIST ( APPEND LESDataProcessing_cflibs Framework )
  LIST ( APPEND LESDataProcessing_cflibs LES)

  CF_ADD_PLUGIN_LIBRARY ( LESDataProcessing )

  CF_WARN_ORPHAN_FILES()

ENDIF ( CF_HAVE_LESMODELS )

