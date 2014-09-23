IF ( CF_HAVE_LIBTOOL ) # needs libtool to load the libs

  INCLUDE_DIRECTORIES ( ${LIBTOOL_INCLUDE_DIR} )

  # these files do not produce orphan warnings
  LIST ( APPEND OPTIONAL_dirfiles
          RunNschemeTVMET.cxx
          NschemeTVMET.hh
          NschemeTVMET.ci
          RunRoeFluxFlens.cxx
          RunRoeFluxFlens.hh
  )

  IF ( CF_HAVE_TVMET )
    LIST ( APPEND RunTools_includedirs ${TVMET_INCLUDE_DIR} )
    ADD_DEFINITIONS(-DTVMET_OPTIMIZE)
    LIST ( APPEND BenchmarkLoader_files
          RunNschemeTVMET.cxx
          NschemeTVMET.hh
          NschemeTVMET.ci
        )
  ENDIF ( CF_HAVE_TVMET )

  IF ( CF_HAVE_BLITZ )
    LIST ( APPEND RunTools_includedirs ${BLITZ_INCLUDE_DIR} )
  ENDIF ( CF_HAVE_BLITZ )

  IF ( CF_HAVE_FLENS )
    LIST ( APPEND RunTools_includedirs ${FLENS_INCLUDE_DIR} )
    LIST ( APPEND BenchmarkLoader_includedirs ${FLENS_INCLUDE_DIR} )
    LIST ( APPEND BenchmarkLoader_libs ${FLENS_LIBRARY} )
    LIST ( APPEND BenchmarkLoader_files RunRoeFluxFlens.cxx  RunRoeFluxFlens.hh )
  ENDIF ( CF_HAVE_FLENS )

  LIST ( APPEND RunTools_files RunTools.cxx )

  LIST ( APPEND BenchmarkLoader_files
        BenchmarkLoader.cxx
        CodeBuilders.hh
        TemplateClassLoader.hh
        convert_to_string.hh
        RunRoeFlux.hh
        RunCFMathTools.hh
        RunBlitz.hh
        RunFlens.hh
        RunTinyBlitz.hh
        RunRoeFluxTinyBlitz.hh
        RunRoeFluxCFLET.hh
        AutoTemplateLoader.hh
        RunRoeFluxBlitz.hh
        RunTVMET.hh
        BaseRunLib.hh
        RunRoeFluxTVMET.hh
        RunNscheme.hh
        RunNscheme.cxx
        RunNschemeD.cxx
        Nscheme.hh
        Nscheme.ci
        NschemeD.hh
        NschemeD.ci
        MatrixInverterTVMET.hh
      )

  LIST ( APPEND BenchmarkLoader_libs Environment MathTools ${LIBTOOL_LIBRARY} )
  LIST ( APPEND RunTools_libs        Environment MathTools ${LIBTOOL_LIBRARY} )

  CF_ADD_PLUGIN_APP ( BenchmarkLoader )
  CF_ADD_PLUGIN_APP ( RunTools )

  CF_WARN_ORPHAN_FILES()

ENDIF ( CF_HAVE_LIBTOOL )
