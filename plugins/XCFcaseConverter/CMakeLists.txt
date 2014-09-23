LIST ( APPEND xcfcase_converter_files main.cxx )
LIST ( APPEND xcfcase_converter_libs ${CF_KERNEL_LIBS} ${CF_KERNEL_STATIC_LIBS}  ${CF_Boost_LIBRARIES}  )

CF_ADD_PLUGIN_APP ( xcfcase_converter )

CF_WARN_ORPHAN_FILES()
ADD_SUBDIRECTORY ( testcases )
