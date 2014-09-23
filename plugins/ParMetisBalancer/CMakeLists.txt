LIST ( APPEND ParMetisBalancer_files
DataStorage.hh
DataStorage.cxx
ParMetisBalancer.hh
ParMetisBalancer.cxx
ParMetisBalancerData.hh
ParMetisBalancerData.cxx
ParMetisBalancerModule.hh
StdSetup.hh
StdSetup.cxx
StdUnSetup.hh
StdUnSetup.cxx
)

LIST ( APPEND OPTIONAL_dirfiles StdRepart.hh StdRepart.cxx ) # avoids warning of orphan when parmeti is not present
IF(CF_HAVE_PARMETIS)
  LIST ( APPEND ParMetisBalancer_files StdRepart.hh StdRepart.cxx    )
  LIST ( APPEND ParMetisBalancer_includedirs ${PARMETIS_INCLUDE_DIR} )
  LIST ( APPEND ParMetisBalancer_libs        ${PARMETIS_LIBRARIES}   )
ENDIF(CF_HAVE_PARMETIS)

LIST ( APPEND ParMetisBalancer_cflibs Framework )

CF_ADD_PLUGIN_LIBRARY ( ParMetisBalancer )

CF_WARN_ORPHAN_FILES()
