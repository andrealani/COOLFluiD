IF ( CF_HAVE_LESMODELS )
  LIST ( APPEND LES_files
  LES2DCons.cxx
  LES2DCons.hh
  LES2DPuvt.cxx
  LES2DPuvt.hh
  LES2DVarSet.cxx
  LES2DVarSet.hh
  LES3DCons.cxx
  LES3DCons.hh
  LES3DPvt.cxx
  LES3DPvt.hh
  LES3DVarSet.cxx
  LES3DVarSet.hh
  LESModule.hh
  LESPhysicalModel.ci
  LESPhysicalModel.cxx
  LESPhysicalModel.hh
  LESVarSet.cxx
  LESVarSet.hh
  )

  LIST ( APPEND LES_includedirs ${LESMODELS_INCLUDE_DIR} )
  LIST ( APPEND LES_libs ${LESMODELS_LIBRARIES})

  LIST ( APPEND LES_cflibs Framework NavierStokes)


  CF_ADD_PLUGIN_LIBRARY ( LES )
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/les_config.hh.in ${COOLFluiD_BINARY_DIR}/les_config.hh)

  CF_WARN_ORPHAN_FILES()

ENDIF ( CF_HAVE_LESMODELS )

