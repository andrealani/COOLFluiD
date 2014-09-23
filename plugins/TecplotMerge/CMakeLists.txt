LIST ( APPEND tecplot_merge_files TecplotMerge.cxx )
LIST ( APPEND tecplot_merge_kopacek_files TecplotMergeKopacek.cxx )
LIST ( APPEND coef_merge_files    CoefMerge.cxx )

CF_ADD_PLUGIN_APP ( tecplot_merge_kopacek )
CF_ADD_PLUGIN_APP ( tecplot_merge )
CF_ADD_PLUGIN_APP ( coef_merge )

CF_WARN_ORPHAN_FILES()
